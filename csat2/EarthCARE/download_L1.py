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
EARTHCARE_SERVER = "ftps://ec-pdgs-dissemination1.eo.esa.int:990"

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
    remote_dir = f"/EarthCARE/EarthCAREL1Validated/{product_type}/{baseline}/{y}/{m}/{d}"

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
            cd /EarthCARE/EarthCAREL1Validated/{product_type}/{baseline}/{year}/{month:02d}/{day:02d};
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


