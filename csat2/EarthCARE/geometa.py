from pathlib import Path
from datetime import datetime, timedelta
from csat2 import locator
import csat2.misc.time
from csat2.EarthCARE.utils import get_orbit_date_approx
from csat2.EarthCARE.download import download_file_locations
from csat2.EarthCARE.utils import get_orbit_date_approx, DEFAULT_VERSION
import numpy as np


def create_geometa(year):
    # Geometa files currently have four columns
    # Orbit Year DOY Hour(decimal time)
    files = download_file_locations('ATL_NOM_1B', year=year, frame='A', limit=6000)
    orbdata = [(
        f.split('_')[-1][:5],
        csat2.misc.time.datetime_to_ydh(
            datetime.strptime(f.split('_')[5], '%Y%m%dT%H%M%SZ'))
    ) for f in files]
    
    path = Path(locator.format_filename("EarthCARE", "GEOMETA", year=year))
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, 'w') as f:
        f.write('# Orbit Year DOY Hour\n')
        for data in orbdata:
            f.write(f'{data[0]} {data[1][0]} {data[1][1]} {data[1][2]:.2f}\n')
    return


def read_geometa(year):
    filename = locator.format_filename("EarthCARE", "GEOMETA", year=year)
    return np.genfromtxt(filename, names=True, dtype=['i', 'i', 'i', 'f'])


def get_or_create_earthcare_geometa(
        year,
        orbit=None,
        version=DEFAULT_VERSION
):

    """
    Ensures the EarthCARE GEOMETA file for a given year exists.
    If orbit_number is provided, ensures that orbit is listed.

    Defaults:
    - year = 2025
    - product = CPR_CLD_2A
    - version = AB
    """

    # Step 1: Estimate orbit date (if orbit provided)
    if orbit_number is not None:
        year_est, doy_est, date_est, _ = get_orbit_date_approx(orbit_number)
        year = date_est.year  # override year to match actual orbit

    # Step 2: Construct GEOMETA path based on final year
    path = Path(locator.format_filename("EarthCARE", "GEOMETA", year=year))

    # Early exit if file exists and no orbit check is needed
    if path.exists() and orbit_number is None:
        return path

    # Ensure folder exists
    path.parent.mkdir(parents=True, exist_ok=True)

    # ==== DELAYED IMPORT here to avoid circular issue
    from csat2.EarthCARE.download import download_file_locations
    orbit_str = f"{orbit_number:05d}"
    entries = []

    # Scan ±1 day around estimated date
    for offset in [-1, 0, 1]:
        year_offset, doy_offset = csat2.misc.time.doy_step(
            year_est, doy_est, offset)
        try:
            files = download_file_locations(product=product,
                                            year=year_offset, doy=doy_offset,
                                            version=version)
        except Exception as e:
            print(f"[WARNING] Could not list files for {d.date()}: {e}")
            continue

        for f in files:
            if f"_{orbit_str}" in f:
                entries.append(Path(f).name)

    if not entries:
        raise FileNotFoundError(
            f"No files found for orbit {orbit_number} near {year_est}, {doy_est}")

    # Write or append to GEOMETA file
    entries = sorted(set(entries))
    if not path.exists():
        with open(path, "w") as f:
            f.write("\n".join(entries) + "\n")
    else:
        with open(path, "r+") as f:
            existing = set(f.read().splitlines())
            new_entries = [e for e in entries if e not in existing]
            if new_entries:
                f.write("\n" + "\n".join(new_entries) + "\n")

    return path

