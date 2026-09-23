import requests
from pathlib import Path
from csat2.download.earthdata import get_token, geturl
from bs4 import BeautifulSoup
from urllib.parse import urljoin

# Base URL for CERES SYN1deg-1Hour data



BASE_URL = 'https://asdc.larc.nasa.gov/data/CERES/SYN1deg-1Hour/Terra-Aqua-NOAA20_Edition4B/'

CMR_BASE_URL = ("https://cmr.earthdata.nasa.gov/virtual-directory/collections/C3181056140-LARC_CLOUD/temporal/"
)


def _html_to_text(html):
    if isinstance(html, bytes):
        return html.decode("utf-8", errors="replace")
    return str(html)

def find_granule_url(year: int, month: int, dom: int):
    """
    Authenticated directory listing fetch + parse to find hourly granule for given date.
    Returns full URL to the .hdf file (including granule id).
    """

    date_str = f"{year}{month:02d}{dom:02d}"

    old_dir_url = f"{BASE_URL}{year}/{month:02d}/"

    # Use csat2 geturl to fetch the directory HTML authenticated.
    try:

        html = geturl(old_dir_url, out=None, quiet=True)
        text = _html_to_text(html)

    

        soup = BeautifulSoup(text, "html.parser")
        files = [a["href"] for a in soup.find_all("a") if a.get("href", "").endswith(".hdf")]


        # Match pattern: *_<granuleID>.YYYYMMDD.hdf (I can't work out how the granule IDs are assigned)
        for f in files:
            if f.endswith(f".{date_str}.hdf"):
                return urljoin(old_dir_url, f)
            
    except Exception as e:
        print('old file lokup failed, trying CMR method')

    cmr_day_url = f"{CMR_BASE_URL}{year}/{month:02d}/{dom:02d}"

    try:
        html = geturl(cmr_day_url, out=None, quiet=True)
        text = _html_to_text(html)

        soup = BeautifulSoup(text, "html.parser")

        candidates = []
        for a in soup.find_all("a"):
            href = a.get("href", "")
            label = a.get_text(strip=True)

            if date_str in href or date_str in label:
                candidates.append(href)

        if not candidates:
            preview = "\n".join(a.get("href", "") for a in soup.find_all("a")[:10])
            raise ValueError(
                f"No matching CMR granule link found for {date_str} in {cmr_day_url}\n"
                f"First links found:\n{preview}"
            )

        granule_url = urljoin(cmr_day_url + "/", candidates[0])

        # The CMR page often links to the granule landing/redirect URL.
        # geturl should follow redirects/auth when downloading.
        return granule_url

    except Exception as e:
        raise ValueError(
            f"No CERES hourly file found for {date_str}\n"
            f"Tried old URL:\n  {old_dir_url}\n"
            f"Tried CMR URL:\n  {cmr_day_url}\n"
            f"CMR error: {e}"
        )


def download_files(year: int, month: int, dom: int, local_path: Path):
    """
    Download CERES SYN1deg-1Hour data for a specific date.
    - Dynamically finds the granule id using authenticated directory listing.
    - Downloads using geturl() which handles Earthdata auth and progress bar.
    """
    local_path = Path(local_path)
    local_path.parent.mkdir(parents=True, exist_ok=True)

    try:
        file_url = find_granule_url(year, month, dom)
        filename = file_url.split("/")[-1]

        print(f"Found granule: {filename}")
        print(f"Downloading from: {file_url}")
        print(f"Saving to: {local_path}")

        # Download the file into the local path, show progress bar
        with open(local_path, "wb") as f:
            geturl(file_url, out=f, quiet=False)

        print(f"Downloaded: {local_path}")
        print(f"File size: {local_path.stat().st_size/1024**2:.1f} MB")

    except Exception as e:
        print(f"Download failed: {e}")
        if local_path.exists():
            try:
                local_path.unlink()
                print(f"Removed partial file: {local_path}")
            except Exception:
                pass
        raise


