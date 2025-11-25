import requests
from pathlib import Path
from csat2.download.earthdata import get_token, geturl
from bs4 import BeautifulSoup

# Base URL for CERES SYN1deg-1Hour data



BASE_URL = 'https://asdc.larc.nasa.gov/data/CERES/SYN1deg-1Hour/Terra-Aqua-NOAA20_Edition4B/'

def find_granule_url(year: int, month: int, dom: int):
    """
    Authenticated directory listing fetch + parse to find hourly granule for given date.
    Returns full URL to the .hdf file (including granule id).
    """
    dir_url = f"{BASE_URL}{year}/{month:02d}/"

    # Use csat2 geturl to fetch the directory HTML authenticated.

    html = geturl(dir_url, out=None, quiet=True)
    # geturl returns bytes or content; ensure we have text
    if isinstance(html, bytes):
        text = html.decode("utf-8", errors="replace")
    else:
        text = str(html)

    soup = BeautifulSoup(text, "html.parser")
    files = [a["href"] for a in soup.find_all("a") if a.get("href", "").endswith(".hdf")]

    date_str = f"{year}{month:02d}{dom:02d}"  # YYYYMMDD

    # Match pattern: *_<granuleID>.YYYYMMDD.hdf (I can't work out how the granule IDs are assigned)
    for f in files:
        if f.endswith(f".{date_str}.hdf"):
            return dir_url + f

    # If none matched, raise debug message listing available few files, this is highly unlikely to occur normally.
    preview = "\n".join(files[:5]) if files else "(no .hdf entries found)"
    raise ValueError(
        f"No CERES hourly file found for {date_str} in {dir_url}\n"
        f"Available .hdf files (first 5):\n{preview}"
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
        print(f"Saving to:       {local_path}")

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


