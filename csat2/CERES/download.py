import requests
from pathlib import Path
from csat2.download.earthdata import get_token, geturl
from bs4 import BeautifulSoup

# Base URL for CERES SYN1deg-1Hour data



BASE_URL = 'https://asdc.larc.nasa.gov/data/CERES/SYN1deg-1Hour/Terra-Aqua-NOAA20_Edition4B/'

def download_files(year: int, month: int, dom: int, local_path: Path):
    """
    Download CERES SYN1deg-1Hour data for a specific date using Earthdata username/password
    """
    remote_filename = f'CER_SYN1deg-1Hour_Terra-Aqua-NOAA20_Edition4B_407412.{year}{month:02d}{dom:02d}.hdf'
    file_url = BASE_URL + f'{year}/{month:02d}/{remote_filename}'

    local_path = Path(local_path)
    local_path.parent.mkdir(parents=True, exist_ok=True)

    print(f"Downloading: {remote_filename}")
    print(f"From: {file_url}")
    print(f"To: {local_path}")

    try:
        with open(local_path, "wb") as f:
            geturl(file_url, out=f, quiet=False)  # automatically uses SessionWithHeaderRedirection

        print(f"Successfully downloaded: {local_path}")
        print(f"Final file size: {local_path.stat().st_size / 1024**2:.1f} MB")

    except Exception as e:
        print(f"Download failed: {e}")
        if local_path.exists():
            local_path.unlink()
            print(f"Removed partial file: {local_path}")
        raise



