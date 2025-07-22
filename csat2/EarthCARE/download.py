"""
EarthCARE download module for the csat2 library.

- Uses `locator.get_folder()` to determine local storage paths
- Reads credentials from ~/.csat2/earthcare_auth.json
"""

import os
import json
import re
import subprocess
from pathlib import Path
import time
import requests
import pickle
import zipfile
from csat2 import locator
import csat2.misc.time
from tqdm import tqdm
import logging
from csat2.EarthCARE.utils import DEFAULT_VERSION, get_product_level

log = logging.getLogger(__name__)

def load_earthcare_auth():
    """
    Loads EarthCARE credentials from ~/.csat2/earthcare_auth.json.
    Raises FileNotFoundError or json.JSONDecodeError if invalid.
    """
    path = os.path.expanduser("~/.csat2/earthcare_auth.json")
    if not os.path.exists(path):
        raise FileNotFoundError(f"To download data, place your ESA/EarthCARE credentials in: {path}")
    with open(path) as f:
        return json.load(f)


def download_file_locations(product,
                            year=None, doy=None, hour=None, minute=None,
                            dtime=None,
                            orbit=None, frame=None,
                            version=DEFAULT_VERSION):
    """
    List available ZIP filenames for an EarthCARE Level-2 product on a given date.
    """
    product_level = get_product_level(product)
    url = (f'https://eocat.esa.int/eo-catalogue/collections/EarthCAREL{product_level}Validated/items?'+
           f'&limit=2000&productType={product}')

    dataflag = False
    if year and doy:
        if dtime:
            raise ValueError('Use either year/doy or datetime')
        # Format date strings
        year, month, day = csat2.misc.time.doy_to_date(year, doy)
        tstr = f'{year}-{month:0>2}-{day:0>2}'
        if hour:
            if minute is None:
                raise ValueError('Minute must be specified if hour is')
            tstr = f'{year}-{month:0>2}-{day:0>2}T{hour:0>2}:{minute:0>2}:00.00Z'
        url += f'&datetime={tstr}'
        dataflag = True
    if dtime:
        url += f'&datetime={dtime.replace(microsecond=0).isoformat()}Z'
        dataflag = True
    if orbit:
        url += f'&orbitNumber={orbit}'
        dataflag = True
    if frame:
        url += f'&frame={frame}'
        dataflag = True
    if version:
        url += f'&productVersion={version}'
        dataflag = True
        
    if dataflag == False:
        raise ValueError('You need to select fewer files, try a year/doy or orbit')
    
    # This request doesn't accept authorisation
    response = requests.Session().get(url, stream=True)
    file_dict = json.loads(response.text)

    if file_dict['numberMatched'] == 0:
        raise ValueError('No EarthCare files for the given parameters')

    output_names = []
    for feature in file_dict['features']:
        if feature['collection'] == f'EarthCAREL{product_level}Validated':
            output_names.append(feature['id'])

    return sorted(output_names)


def get_auth_cookies(oads_hostname,
                     expiry_limit=3600):
    '''Getting new authentication cookies each time we make a request is slow.
    This allows them to persist on disk until the need to be renewed.'''

    cookie_file = os.path.expanduser(f"~/.csat2/{oads_hostname}.cjar")
    try:
        with open(cookie_file, 'rb') as f:
            auth_cookies = pickle.load(f)
            log.debug('Auth cookie found and loaded')
    except FileNotFoundError:
        log.debug('No auth cookie file')
        auth_cookies = refresh_auth_cookies(oads_hostname)

    expires = None
    for cookie in auth_cookies:
        if cookie.name == '_saml_idp':
            if (cookie.expires - time.time()) < expiry_limit:
                log.debug('Auth cookie expiring')
                auth_cookies = refresh_auth_cookies(oads_hostname)
            break
    else:
        log.debug('No valid auth cookie')
        auth_cookies = refresh_auth_cookies(oads_hostname)

    with open(cookie_file, 'wb') as f:
        pickle.dump(auth_cookies, f)
    log.debug('Stored new auth cookie')
    os.chmod(cookie_file, 0o600)
    log.debug('Fixed auth cookie file permissions')
    
    return auth_cookies


def refresh_auth_cookies(oads_hostname):
    # SAML logins require three requests
    # https://stackoverflow.com/questions/52618451/python-requests-saml-login-redirect
    # See also https://github.com/koenigleon/oads-download/ for EarthCARE specific details/servers
    
    # Request 1: The 'target' - in this case the login page, but it could be anything
    response1 = requests.get(f"https://{oads_hostname}/oads/access/login")

    response1_cookies = response1.cookies
    for r in response1.history:
        response1_cookies = requests.cookies.merge_cookies(response1_cookies, r.cookies)
    sessionDataKey = re.search(
        r"sessionDataKey=([\:\/\w-]*)",
        response1.content.decode('ascii')).groups(0)[0]
    log.info('SessionKey Identified')

    # Request 2: We were redirected to the login platform, so make a new request here
    # using the user data and the session key that we got from the first request
    earthcare_auth = load_earthcare_auth()
    request2_data = {
        "username": earthcare_auth['username'],
        "password": earthcare_auth['password'],
        "sessionDataKey": sessionDataKey,
        "tocommonauth": "true" # Apparently required for EarthCARE
    }
    
    response2 = requests.post("https://eoiam-idp.eo.esa.int/samlsso",
                              data=request2_data,
                              cookies=response1_cookies,
                              )

    # Extract the necessary form info (hidden parameters here)
    # Note that this regex form is fragile, consider upgrading in the future
    relayState = re.search(
        r"RelayState.*?value=\'([\w\:\/\.-]*?)\'",
        response2.content.decode('ascii')).groups(0)[0]
    SAMLResponse = re.search(
        r"SAMLResponse.*?value=\'([\w\+=]*)",
        response2.content.decode('ascii')).groups(0)[0]
    saml_redirect_url = re.search(
        r"samlsso-response-form.*action=\"([\w\.\"\:\/-]*)\"",
        response2.content.decode('ascii')).groups(0)[0]
    log.info('SAML authentication successful')
        

    # Request 3: Send the SAML response from the authentication platform
    # to the service - in this case the data server. Once this is validated,
    # we get a cookie we can use for downloading files.
    request3_data = {
        "RelayState": relayState,
        "SAMLResponse": SAMLResponse,
    }

    response3 = requests.post(saml_redirect_url,
                              data=request3_data
                              )

    # Lovely tasty authentication cookies for downloading data
    response3_cookies = response3.cookies
    for r in response3.history:
        response3_cookies = requests.cookies.merge_cookies(response3_cookies, r.cookies)
    log.info('SAML authentication validated by data server')

    return response3_cookies


def download(product, year=None, doy=None, orbit=None, frame=None,
             version=DEFAULT_VERSION, force_redownload=False, quiet=False):
    """
    Downloads EarthCARE files for a specified product and date/orbit/frame.

    Files are downloaded via HTTPS into a local folder defined by csat2's locator system.

    Args:
        product_type (str): EarthCARE product name (e.g. "CPR_CLD_2A").
        version (str): Version code (e.g. "AB").
        year (int): Year of data (e.g. 2025).
        doy (int): Day of year in NASA format (1st Jan is DOY1)

    Returns:
        list[str]: List of successfully downloaded filenames.
    """
    file_locations = download_file_locations(product, year, doy, orbit, frame, version)

    session = requests.Session()
    session.headers["user-agent"] = 'csat2'
    auth_cookies = get_auth_cookies('ec-pdgs-dissemination1.eo.esa.int')
    
    product_level = get_product_level(product)
    
    for file_location in file_locations:
        url = (f'https://ec-pdgs-dissemination1.eo.esa.int/oads/data/'+
               f'EarthCAREL{product_level}Validated/{file_location}.ZIP')

        # We could be downloading an orbit where the frames go into a different day
        # We need to check the folder for each new file
        timestr = file_location.split('_')[5]
        year = int(timestr[:4])
        _, doy = csat2.misc.time.date_to_doy(
            year, int(timestr[4:6]), int(timestr[6:8]))

        local_folder = locator.get_folder(
            "EarthCARE", product=product,
            year=year, doy=doy,
            version=version,
        )
        os.makedirs(local_folder, exist_ok=True)

        newfile = f"{local_folder}/{file_location}.ZIP"
        if (not os.path.exists(newfile.replace('.ZIP', '.h5'))) or force_redownload:
            if force_redownload:
                os.remove(newfile)
       
            response = session.get(
                url,
                stream=True,
                cookies=auth_cookies)

            if response.status_code != 200:
                raise ValueError(f"Download failed status:{response.status_code}")
            total = int(response.headers.get("content-length", 0))
            log.debug("Got response")

            with open(newfile, 'w+b') as out:
                with tqdm(
                        desc=os.path.basename(url),
                        total=total,
                        unit="iB",
                        unit_scale=True,
                        unit_divisor=1024,
                        disable=quiet,
                    ) as bar:
                        blocksize = max(4096, total // 100)
                        for data in response.iter_content(chunk_size=blocksize):
                            size = out.write(data)
                            bar.update(size)
            with zipfile.ZipFile(newfile, 'r') as zip_ref:
                zip_ref.extractall(local_folder)
            os.remove(newfile)
        else:
            log.info("Skipping {}".format(os.path.basename(url)))

def check(product,
          year, doy,
          orbit=None, frame=None,
          version=DEFAULT_VERSION):
    """
    Check if a specific EarthCARE file exists locally.

    Args:
        product_type (str): EarthCARE product name (e.g., "CPR_CLD_2A").
        version (str): Version code (e.g., "AB").
        year (int): Year of data.
        month (int): Month of data.
        day (int): Day of data.
        orbit (int): Orbit number.
        orbit_id (str): Orbit ID suffix (e.g., "H").

    Returns:
        bool: True if file is found locally, False otherwise.

    Example filename:
        ECA_EXAB_CPR_CLD_2A_20250320T195315Z_20250320T223619Z_04603H.ZIP

    Example usage:
        check(product_type='CPR_CLD_2A', version='AB',
              year=2025, month=3, day=20, orbit=4603, orbit_id='H')
    """
    filename = locator.search("EarthCARE", product,
                              version=version,
                              year=year, doy=doy,
                              orbit=orbit, frame=frame,
                              )
    if len(filename) == 1:
        return True
    else:
        return False


def list_local_files_for_day(product,
                             year, doy,
                             version=DEFAULT_VERSION):
    """
    List all EarthCARE files available locally for the given date.

    Returns:
        list[Path]: Paths to matching local files.
    """
    filenames = locator.search("EarthCARE", product,
                               year=year, doy=doy,
                               version=version,
                               orbit='*****', #Note that orbit wil not expand by default
                               frame='*',
                               )
    return sorted(filenames)


def test_connection():
    """
    Verifies EarthCARE connectivity and credentials by ensuring a valid SAML
    authentication response.
    """
    try:
        _ = refresh_auth_cookies('ec-pdgs-dissemination1.eo.esa.int')
        print("EarthCARE connection test succeeded.")
    except Exception as e:
        print(f"EarthCARE connection test failed: {e}")

if __name__ == "__main__":
    test_connection()

