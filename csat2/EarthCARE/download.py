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


# helpers for extending download
def _extract_zip(zip_path, out_dir, force_clean=False):
    """
    Extract a downloaded ZIP into out_dir.

    If force_clean=True, only keep .h5 files and remove the ZIP and header file afterwards.
    """
    with zipfile.ZipFile(zip_path, 'r') as zip_ref:
        for zip_info in zip_ref.infolist():
            if zip_info.is_dir():
                continue
            name = os.path.basename(zip_info.filename)
            if force_clean and (not name.lower().endswith('.h5')):
                continue
            zip_info.filename = name
            zip_ref.extract(zip_info, out_dir)

    if force_clean:
        base = os.path.splitext(os.path.basename(zip_path))[0]
        hdr1 = os.path.join(out_dir, base + '.HDR')
        hdr2 = os.path.join(out_dir, base + '.hdr')
        if os.path.exists(hdr1):
            os.remove(hdr1)
        if os.path.exists(hdr2):
            os.remove(hdr2)
        if os.path.exists(zip_path):
            os.remove(zip_path)


def _file_date_from_id(file_id):
    """
    Extract (year, doy) from an EarthCARE file id.
    """
    timestr = file_id.split('_')[5]
    year, doy = csat2.misc.time.date_to_doy(
        int(timestr[:4]), int(timestr[4:6]), int(timestr[6:8]))
    return year, doy


def _drop_processing_time(file_id):
    """
    It removes the ProcessingTime from the file ID
    """
    parts = file_id.split('_')
    if len(parts) < 8:
        return file_id
    return '_'.join(parts[:6] + parts[7:])


def _processing_time(file_id):
    """
    It only extract the ProcessingTime from the file ID
    """
    parts = file_id.split('_')
    if len(parts) < 8:
        return ''
    return parts[6]


def _frame_token(file_id):
    return file_id.split('_')[-1][-1]


def _pick_latest_iterations(file_locations):
    """
    to avoid downloading one data file multiple time with different <ProcessingTime>
    Keep the file with the latest processing timestamp.
    If only the immediate version exists for that file -> use it.
    """
    latest = {}
    for file_location in file_locations:
        key = _drop_processing_time(file_location['id'])
        if (key not in latest) or (_processing_time(file_location['id']) > _processing_time(latest[key]['id'])):
            latest[key] = file_location
    return list(latest.values())


def _find_local_iterations(local_folder, file_id):
    """
    it checks if already a file is downloaded, it can be different ProcessingTime though.
    """
    if not os.path.exists(local_folder):
        return []

    wanted_key = _drop_processing_time(file_id)
    matches = []

    for name in os.listdir(local_folder):
        if not (name.endswith('.h5') or name.endswith('.ZIP')):
            continue
        stem = os.path.splitext(name)[0]
        if _drop_processing_time(stem) == wanted_key:
            matches.append(os.path.join(local_folder, name))

    return sorted(matches)


def _remove_local_iterations(local_folder, file_id):
    """
    remove the local file irrespective of ProcessingTime
    """
    for path in _find_local_iterations(local_folder, file_id):
        try:
            os.remove(path)
        except FileNotFoundError:
            pass


def download(product, 
             year=None, doy=None,
             orbit=None, frames=None,
             baselines=DEFAULT_BASELINE,
             max_files=None,
             force_clean=False,
             force_redownload=False,
             quiet=False,
             method='MAAP'):
    """
    Download EarthCARE files for a given product and date/orbit selection.

    Files are queried from the EarthCARE catalogue and downloaded via HTTPS into
    the local folder structure defined by csat2's locator system. The function
    supports both MAAP and OADS download methods.

    You can select files by:
    - year + doy
    - optionally orbit
    - optionally one or more frames

    If multiple baselines are provided, they are tried in order until matching
    files are found. If multiple processing-time iterations exist for the same
    logical granule, only the newest iteration is kept. Downloaded ZIP files are
    then extracted locally, and optional cleanup can remove temporary ZIP/HDR
    files while keeping the HDF5 data.

    Args:
    product (str):
        EarthCARE product name, e.g. "ATL_NOM_1B" or "CPR_CLD_2A".
    year (int, optional):
        Year of the requested data.
    doy (int, optional):
        Day of year (1 = 1 January).
    orbit (int or str, optional):
        Orbit number used to filter the search.
    frames (str or list[str], optional):
        Single frame or list of frames to download, e.g. "F" or ["F", "G"].
    baselines (str or list[str], optional):
        Single baseline or list of baselines to try in order.
    max_files (int, optional):
        Maximum number of new files to download.
    force_clean (bool, optional):
        If True, remove ZIP/HDR files after extraction and keep only h5.
    force_redownload (bool, optional):
        If True, remove existing local iterations and download again.
    quiet (bool, optional):
        If True, disable the progress bar output.
    method (str, optional):
        Download source, either "MAAP" or "OADS".

    Returns:
    list[str]:
        List of successfully downloaded file IDs.
    """
    
    if isinstance(baselines, str):
        baselines = [baselines]
    if isinstance(frames, str):
        frames = [frames]

    session = requests.Session()
    session.headers["user-agent"] = 'csat2'
    if method == 'OADS':
        auth_cookies = esa_cookie.get_auth_cookies('ec-pdgs-dissemination1.eo.esa.int')
    elif method == 'MAAP':
        token = esa_maap_token.refresh_token()
    else:
        raise ValueError("Download method must be 'OADS' or 'MAAP'")

    product_level = get_product_level(product)
    fetched = []

    if frames is None:
        remaining_frames = None
    else:
        remaining_frames = list(frames)

    query_limit = 300 if max_files is None else max(300, max_files)

    for baseline in baselines:
        if max_files is not None and len(fetched) >= max_files:
            break

        if remaining_frames is None:
            try:
                file_locations = download_file_locations(
                    product,
                    year=year, doy=doy,
                    orbit=orbit,
                    baseline=baseline,
                    limit=query_limit,
                )
            except Exception:
                continue
        else:
            file_locations = []
            for frame in remaining_frames:
                try:
                    one = download_file_locations(
                        product,
                        year=year, doy=doy,
                        orbit=orbit,
                        frame=frame,
                        baseline=baseline,
                        limit=query_limit,
                    )
                    file_locations.extend(one)
                except Exception:
                    pass

        if len(file_locations) == 0:
            continue

        file_locations = _pick_latest_iterations(file_locations)

        file_locations = sorted(file_locations, key=lambda x: x['id'])

        processed_frames = set()

        for file_location in file_locations:
            if max_files is not None and len(fetched) >= max_files:
                break

            if method == 'OADS':
                url = (f'https://ec-pdgs-dissemination1.eo.esa.int/oads/data/' +
                       f'EarthCAREL{product_level}Validated/{file_location["id"]}.ZIP')
            else:
                url = file_location['maap_zip']

            file_year, file_doy = _file_date_from_id(file_location["id"])

            local_folder = locator.get_folder(
                'EarthCARE', product=product,
                baseline=baseline,
                year=file_year, doy=file_doy,
            )

            os.makedirs(local_folder, exist_ok=True)

            existing_local = _find_local_iterations(local_folder, file_location['id'])

            if existing_local and not force_redownload:
                log.info('Skipping %s (an iteration already exists locally)', file_location['id'])
                processed_frames.add(_frame_token(file_location['id']).upper())
                continue

            if force_redownload and existing_local:
                _remove_local_iterations(local_folder, file_location['id'])

            newfile = f"{local_folder}/{file_location['id']}.ZIP"

            try:
                os.remove(newfile)
            except FileNotFoundError:
                pass

            if method == 'OADS':
                response = session.get(url, stream=True, cookies=auth_cookies)
            else:
                response = session.get(url, stream=True, headers={'Authorization': f'Bearer {token}'})

            if response.status_code != 200:
                raise ValueError(f'Download failed status:{response.status_code}')

            total = int(response.headers.get('content-length', 0))

            with open(newfile, 'w+b') as out:
                with tqdm(
                    desc=os.path.basename(url),
                    total=total,
                    unit='iB',
                    unit_scale=True,
                    unit_divisor=1024,
                    disable=quiet,
                ) as bar:
                    blocksize = max(4096, total // 100) if total > 0 else 4096
                    for data in response.iter_content(chunk_size=blocksize):
                        size = out.write(data)
                        bar.update(size)

            _extract_zip(newfile, local_folder, force_clean=force_clean)
            fetched.append(file_location['id'])
            processed_frames.add(_frame_token(file_location['id']).upper())

        if remaining_frames is not None:
            remaining_frames = [frame for frame in remaining_frames if frame.upper() not in processed_frames]
            if len(remaining_frames) == 0:
                break
        else:
            break

    return fetched
             


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

