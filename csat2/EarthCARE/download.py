"""
EarthCARE download module for the csat2 library.

- Uses `locator.get_folder()` to determine local storage paths
- Reads credentials from ~/.csat2/earthcare_auth.json
- Lists and downloads ZIP files via FTPS using lftp
"""

import os
import json
import subprocess
from pathlib import Path
from csat2 import locator
from .utils import DEFAULT_BASELINE, DEFAULT_PRODUCT_TYPE




# Load credentials
#CREDS_PATH = os.path.expandvars("${HOME}/.csat2/earthcare_auth.json")
#with open(CREDS_PATH) as cf:
#    earthcare_auth = json.load(cf)


def load_earthcare_auth():
    """
    Loads EarthCARE credentials from ~/.csat2/earthcare_auth.json.
    Raises FileNotFoundError or json.JSONDecodeError if invalid.
    """
    path = os.path.expanduser("~/.csat2/earthcare_auth.json")
    if not os.path.exists(path):
        raise FileNotFoundError(f"[ERROR] Credentials file not found: {path}")
    with open(path) as f:
        return json.load(f)



# Path to the lftp binary; can be overridden via the LFTP_BIN environment variable
LFTP_BIN = os.environ.get("LFTP_BIN", "lftp")

# EarthCARE FTPS endpoint for ESA data downloads
EARTHCARE_SERVER = "ftps://ec-pdgs-dissemination2.eo.esa.int:990"



def download_file_locations(product_type=DEFAULT_PRODUCT_TYPE,
                            baseline=DEFAULT_BASELINE,
                            year=None, month=None, day=None):
    """
    List available ZIP filenames for an EarthCARE Level-2 product on a given date.
    """
    
    if year is None or month is None or day is None:
        raise ValueError(
            "Missing required date inputs: year, month, and day must all be specified.\n"
            "Example usage:\n"
            "  download_file_locations(product_type='CPR_CLD_2A', baseline='AB', year=2025, month=3, day=20)"
        )
    
    
    # Format date strings
    y, m, d = f"{year:04d}", f"{month:02d}", f"{day:02d}"
    remote_dir = f"/EarthCARE/EarthCAREL2Validated/{product_type}/{baseline}/{y}/{m}/{d}"

    # Build lftp commands to list files
    cmd = f"""
        set ssl:verify-certificate no;
        set ftp:ssl-force true;
        set ftp:ssl-protect-data true;
        cd {remote_dir};
        cls -1 *.ZIP;
        bye
    """
    auth = load_earthcare_auth()
    proc = subprocess.run(
        [LFTP_BIN, "-u", f"{auth['username']},{auth['password']}",
         EARTHCARE_SERVER, "-e", cmd],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        text=True, check=True
    )
    return sorted(proc.stdout.splitlines())

def download(product_type=DEFAULT_PRODUCT_TYPE, baseline=DEFAULT_BASELINE,
             year=None, month=None, day=None, max_files=None):
    """
    Downloads EarthCARE ZIP files for a specified product and date.

    Files are downloaded via FTPS into a local folder defined by csat2's locator system.

    Args:
        product_type (str): EarthCARE product name (e.g. "CPR_CLD_2A").
        baseline (str): Baseline code (e.g. "AB").
        year (int): Year of data (e.g. 2025).
        month (int): Month of data (1–12).
        day (int): Day of data (1–31).
        max_files (int, optional): Maximum number of missing files to download.

    Returns:
        list[str]: List of successfully downloaded filenames.
    """
    if year is None or month is None or day is None:
        raise ValueError(
            "Missing required date inputs: year, month, and day must all be specified.\n"
            "Example: download(product_type='CPR_CLD_2A', baseline='AB', year=2025, month=3, day=20, max_files=2)"
        )


    zips = download_file_locations(product_type, baseline, year, month, day)
    if not zips:
        raise ValueError(f"No remote .ZIP files for {product_type} {baseline} on {year}-{month:02d}-{day:02d}")

    local_root = locator.get_folder(
        "EARTHCARE", product=product_type,
        baseline=baseline,
        year=year, month=month, day=day
    )
    os.makedirs(local_root, exist_ok=True)

    missing = [f for f in zips if not os.path.exists(os.path.join(local_root, f))]
    to_fetch = missing if max_files is None else missing[:max_files]
    fetched = []

    auth = load_earthcare_auth()
    for fname in to_fetch:
        zip_path = os.path.join(local_root, fname)
        h5_path = zip_path.replace(".ZIP", ".h5")

        # Skip if either .ZIP or .h5 exists
        if os.path.exists(zip_path) or os.path.exists(h5_path):
           print(f" Skipping {fname} (ZIP or H5 already exists)")
           continue

        # Download with lftp
        cmd = f"""
            set ssl:verify-certificate no;
            set ftp:ssl-force true;
            set ftp:ssl-protect-data true;
            cd /EarthCARE/EarthCAREL2Validated/{product_type}/{baseline}/{year}/{month:02d}/{day:02d};
            lcd {local_root};
            get {fname};
            bye
        """

        subprocess.run(
            [LFTP_BIN, "-u", f"{auth['username']},{auth['password']}",
             EARTHCARE_SERVER, "-e", cmd], check=True
        )
        print(f" Downloaded {fname}")
        fetched.append(fname)

    for f in set(zips) - set(fetched):
        print(f"Skipped {f} (already exists)")

    return fetched



def check(product_type=DEFAULT_PRODUCT_TYPE,
          baseline=DEFAULT_BASELINE,
          year=None, month=None, day=None,
          orbit=None, orbit_id=None):
    """
    Check if a specific EarthCARE file exists locally.

    Args:
        product_type (str): EarthCARE product name (e.g., "CPR_CLD_2A").
        baseline (str): Baseline code (e.g., "AB").
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
        check(product_type='CPR_CLD_2A', baseline='AB',
              year=2025, month=3, day=20, orbit=4603, orbit_id='H')
    """
    if None in (year, month, day, orbit, orbit_id):
        raise ValueError("Missing required inputs: year, month, day, orbit, and orbit_id.")

    matches = locator.search("EARTHCARE", product_type,
        baseline=baseline,
        year=year, month=month, day=day,
        orbit=orbit, orbit_id=orbit_id,
    )

    if matches:
        print(f"Found file: {os.path.basename(matches[0])}")
        return True
    else:
        print("No matching file found.")
        return False


def list_local_files_for_day(product_type=DEFAULT_PRODUCT_TYPE,
                             baseline=DEFAULT_BASELINE,
                             year=None, month=None, day=None):
    """
    List all EarthCARE files available locally for the given date.

    Returns:
        list[Path]: Paths to matching local files.
    """
    if None in (year, month, day):
        raise ValueError("Missing required inputs: year, month, and day must be specified.")

    try:
        folder = locator.get_folder(
            "EARTHCARE", product=product_type,
            baseline=baseline,
            year=year, month=month, day=day
        )
    except Exception as e:
        print(f"[ERROR] Could not determine folder path: {e}")
        return []

    if not os.path.exists(folder):
        print(f"Folder not found: {folder}")
        return []

    files = sorted(Path(folder).glob("*"))

    if not files:
        print(f"(No files found in {folder})")
        return []

    print(f"\nFound {len(files)} files in {folder}:\n")
    for f in files:
        print(" -", f.name)

    return files


def test_connection():
    """
    Verifies EarthCARE FTPS connectivity and credentials.
    """
    try:
        _ = download_file_locations("CPR_CLD_2A", "AB", 2025, 4, 17)
        print("EarthCARE connection test succeeded.")
    except Exception as e:
        print(f"EarthCARE connection test failed: {e}")

if __name__ == "__main__":
    test_connection()

