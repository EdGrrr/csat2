"""
EarthCARE download module for the csat2 library.

- Uses `locator.get_folder()` to determine local storage paths
- Reads credentials from ~/.csat2/earthcare_auth.json
"""

import os
import json
import time
import requests
import zipfile
from csat2 import locator
import csat2.misc.time
from tqdm import tqdm
import logging
import fsspec
import xarray as xr
from csat2.EarthCARE.utils import DEFAULT_BASELINE, get_product_level

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
                            baseline=DEFAULT_BASELINE,
                            limit=200):
    """
    List available ZIP filenames for an EarthCARE product on a given date.
    """
    product_level = get_product_level(product)
    url = (f'https://catalog.maap.eo.esa.int/catalogue/collections/EarthCAREL{product_level}Validated_MAAP/items?'+
           f'&limit={limit}&productType={product}')

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
    if baseline:
        url += f'&productVersion={baseline}'
        dataflag = True
        
    if dataflag == False:
        raise ValueError('You need to select fewer files, try a year/doy or orbit')
    
    # This request doesn't accept authorisation
    token = esa_maap_token.refresh_token()

    response = requests.Session().get(url, stream=True, headers={"Authorization": f"Bearer {token}"})
    file_dict = json.loads(response.text)

    if file_dict['numberMatched'] == 0:
        raise ValueError('No EarthCare files for the given parameters')

    output_names = []
    for feature in file_dict['features']:
        if feature['collection'] == f'EarthCAREL{product_level}Validated_MAAP':
            output_names.append(
                {'id': feature['id'],
                 'maap_h5': feature['assets']['enclosure_h5']['href'],
                 'maap_zip': feature['assets']['product']['href'],
                 'maap_thumbnail': feature['assets']['thumbnail']['href']}
            )

    return sorted(output_names, key=lambda x: x["id"])


def get_maap_token():
    """Use OFFLINE_TOKEN to fetch a short-lived access token.

    From the MAAP/ESA example code at https://docs.maap-project.org/en/latest/science/EarthCARE/EarthCARE_access_and_visualize.html"""
    earthcare_auth = load_earthcare_auth()
    
    OFFLINE_TOKEN = earthcare_auth["OFFLINE_TOKEN"]
    CLIENT_ID = earthcare_auth["CLIENT_ID"]
    CLIENT_SECRET = earthcare_auth["CLIENT_SECRET"]

    if not all([OFFLINE_TOKEN, CLIENT_ID, CLIENT_SECRET]):
        raise ValueError("Missing OFFLINE_TOKEN, CLIENT_ID, or CLIENT_SECRET in credentials file ~/.csat2/earthcare_auth.json")

    url = "https://iam.maap.eo.esa.int/realms/esa-maap/protocol/openid-connect/token"
    data = {
        "client_id": CLIENT_ID,
        "client_secret": CLIENT_SECRET,
        "grant_type": "refresh_token",
        "refresh_token": OFFLINE_TOKEN,
        "scope": "offline_access openid"
    }

    response = requests.post(url, data=data)
    response.raise_for_status()

    response_json = response.json()
    access_token = response_json.get('access_token')

    if not access_token:
        raise RuntimeError("Failed to retrieve access token from IAM response")

    return access_token

class ESAMaapToken():
    def __init__(self, expiry_time=3600):
        self.token = None
        self.expiry_time = expiry_time

    def refresh_token(self):
        if (self.token is None):
            log.info('No cookie from this session')
            self.token = (get_maap_token(), time.time())
        elif ((time.time()-self.token[1]) >= self.expiry_time):
            log.info('Cookie is too old')
            self.token = (get_maap_token(), time.time())
        return self.token[0]

esa_maap_token = ESAMaapToken()

def download(product, year=None, doy=None, orbit=None, frame=None,
             baseline=DEFAULT_BASELINE, force_redownload=False, quiet=False, method='MAAP'):
    """
    Downloads EarthCARE files for a specified product and date/orbit/frame.

    Files are downloaded via HTTPS into a local folder defined by csat2's locator system.

    Args:
        product_type (str): EarthCARE product name (e.g. "CPR_CLD_2A").
        baseline (str): Baseline code (e.g. "AB").
        year (int): Year of data (e.g. 2025).
        doy (int): Day of year in NASA format (1st Jan is DOY1)

    Returns:
        list[str]: List of successfully downloaded filenames.
    """
    file_locations = download_file_locations(
        product, year=year, doy=doy, orbit=orbit, frame=frame, baseline=baseline)

    session = requests.Session()
    session.headers["user-agent"] = 'csat2'
    if method == 'OADS':
        auth_cookies = esa_cookie.get_auth_cookies('ec-pdgs-dissemination1.eo.esa.int')
    elif method == 'MAAP':
        token = esa_maap_token.refresh_token()
    else:
        raise ValueError("Download method must be 'OADS' or 'MAAP'")
    
    product_level = get_product_level(product)
    
    for file_location in file_locations:
        if method == 'OADS':
            url = (f'https://ec-pdgs-dissemination1.eo.esa.int/oads/data/'+
                   f"EarthCAREL{product_level}Validated/{file_location['id']}.ZIP")
        elif method == 'MAAP':
            url = f'{file_location["maap_zip"]}'
            
        # We could be downloading an orbit where the frames go into a different day
        # We need to check the folder for each new file
        timestr = file_location['id'].split('_')[5]
        year = int(timestr[:4])
        _, doy = csat2.misc.time.date_to_doy(
            year, int(timestr[4:6]), int(timestr[6:8]))

        local_folder = locator.get_folder(
            "EarthCARE", product=product,
            year=year, doy=doy,
            baseline=baseline,
        )
        os.makedirs(local_folder, exist_ok=True)

        newfile = f"{local_folder}/{file_location['id']}.ZIP"
        if (not os.path.exists(newfile.replace('.ZIP', '.h5'))) or force_redownload:
            try:
                os.remove(newfile)
            except FileNotFoundError:
                pass

            if method == 'OADS':
                response = session.get(
                    url,
                    stream=True,
                    cookies=auth_cookies)
            elif method == 'MAAP':
                response = session.get(
                    url,
                    stream=True,
                    headers={"Authorization": f"Bearer {token}"})

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
                # Remove any subfolders in zip extraction (OADS to MAAP update)
                for zip_info in zip_ref.infolist():
                    if zip_info.is_dir():
                        continue
                    zip_info.filename = os.path.basename(zip_info.filename)
                    log.debug(f'Extracting {zip_info.filename}')
                    zip_ref.extract(zip_info, local_folder)
            #os.remove(newfile)
        else:
            log.info("Skipping {}".format(os.path.basename(url)))

def open_maap_stream(product, orbit, frame, baseline=DEFAULT_BASELINE, fail_multiple=True):
    streams = download_file_locations(product, orbit=orbit, frame=frame, baseline=baseline)
    if len(streams) == 0:
        raise ValueError('No valid files for this granule')
    if (len(streams) > 1) and fail_multiple:
        raise ValueError('Multiple valid files for this granule')
    stream_location = streams[0]['maap_h5']

    token = esa_maap_token.refresh_token()
    
    fs = fsspec.filesystem("https", headers={"Authorization": f"Bearer {token}"})
    f = fs.open(stream_location, "rb")  
    ds = xr.open_dataset(f, engine="h5netcdf", group="ScienceData")
    return ds
            
def check(product,
          orbit=None, frame=None,
          baseline=DEFAULT_BASELINE):
    """
    Check if a specific EarthCARE file exists locally.

    Args:
        product (str): EarthCARE product name (e.g., "CPR_CLD_2A").
        baseline (str): Baseline code (e.g., "AB").
        orbit (int): Orbit number.
        frame (str): Frame (e.g., "H").

    Returns:
        bool: True if file is found locally, False otherwise.

    Example filename:
        ECA_EXAB_CPR_CLD_2A_20250320T195315Z_20250320T223619Z_04603H.ZIP

    Example usage:
        check(product='ATL_NOM_1B', baseline='AE',
              orbit=4603, frame='A')
    """
    raise NotImplementedError('Requires work with geometa')
    
    # filename = locator.search("EarthCARE", product,
    #                           baseline=baseline,
    #                           year=year, doy=doy,
    #                           orbit=orbit, frame=frame,
    #                           )
    # if len(filename) == 1:
    #     return True
    # else:
    #     return False


def list_local_files_for_day(product,
                             year, doy,
                             baseline=DEFAULT_BASELINE):
    """
    List all EarthCARE files available locally for the given date.

    Returns:
        list[Path]: Paths to matching local files.
    """
    filenames = locator.search("EarthCARE", product,
                               year=year, doy=doy,
                               baseline=baseline,
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
        _ = get_maap_token()
        print("EarthCARE token collected.")
    except Exception as e:
        print(f"EarthCARE connection test failed: {e}")

if __name__ == "__main__":
    test_connection()

