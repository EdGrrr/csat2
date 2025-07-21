"""
EarthCARE readfiles helpers
===========================

Utility functions for listing, parsing, and reading EarthCARE Level-2 product
files (HDF5 or ZIP).

Typical filename
----------------
ECA_EXAB_CPR_CLD_2A_20250311T000554Z_20250311T005723Z_04451A.ZIP

Time units inside datasets are usually::

    b'seconds since 2000-01-01 00:00:00.0'
"""

import os
from pathlib import Path
from datetime import datetime, timedelta
import h5py
from csat2 import locator
from csat2.misc.hdf import read_hdf5   
from csat2.EarthCARE.utils import DEFAULT_BASELINE, DEFAULT_PRODUCT_TYPE
from csat2.misc.time import doy_to_date
from csat2.EarthCARE.geometa import get_or_create_earthcare_geometa

from csat2.EarthCARE.utils import get_orbit_date_approx, get_orbit_number_approx
from csat2.EarthCARE.geometa import get_or_create_earthcare_geometa



def files_available(product_type=DEFAULT_PRODUCT_TYPE,
                    baseline=DEFAULT_BASELINE,
                    year=None, month=None, day=None):
    """
    Return a sorted list of local EarthCARE filenames for one day.

    Example: files_available(year=2025, month=3, day=20)
    """
    if None in (year, month, day):
        raise ValueError("year, month, and day must all be provided")

    filenames = locator.search(
        "EARTHCARE", product_type,
        baseline=baseline,
        year=year, month=month, day=day,
        orbit="*****", orbit_id="*",
        exit_first=False,
    )

    if not filenames:
        print(f"No local EarthCARE files for {year:04d}-{month:02d}-{day:02d}")

    return sorted(filenames)

# ------------------------------------------------------------------
def variables_available(product_type=DEFAULT_PRODUCT_TYPE,
                        baseline=DEFAULT_BASELINE,
                        year=None, month=None, day=None):
    """
    Return a dict of all datasets in the first local EarthCARE HDF5 file
    found for the given date.

    Example
    -------
    >>> vars = variables_available(year=2025, month=3, day=20)
    >>> list(vars)
    ['ScienceData/latitude', 'ScienceData/longitude', ...]
    """
    if None in (year, month, day):
        raise ValueError("year, month, and day must all be provided")

    files = locator.search(
        "EARTHCARE", product_type,
        baseline=baseline,
        year=year, month=month, day=day,
        orbit="*****", orbit_id="*",
        exit_first=False,
    )

    h5files = [f for f in files if f.endswith(".h5")]
    if not h5files:
        raise IOError("No .h5 files found for this date (did you unzip the .ZIP archives?).")

    try:
        return read_hdf5(h5files[0])
    except Exception as e:
        raise IOError(f"Failed to read {h5files[0]}: {e}")


def available_orbits(product_type=DEFAULT_PRODUCT_TYPE,
                     baseline=DEFAULT_BASELINE,
                     year=None, month=None, day=None):
    """
    List the orbit identifiers present on disk for one day.

    Example
    -------
    >>> available_orbits(year=2025, month=3, day=20)
    ['04603H', '04604A', ...]
    """
    # Require complete date
    if None in (year, month, day):
        raise ValueError("year, month, and day must all be provided")

    # Look for every matching file (ZIP or H5) on this date
    files = locator.search(
        "EARTHCARE", product_type,
        baseline=baseline,
        year=year, month=month, day=day,
        orbit="*****", orbit_id="*",
        exit_first=False,
    )

    if not files:
        raise IOError("No EarthCARE files found locally for that date.")

    files = sorted(files)
    orbits = [
        os.path.basename(f).split("_")[-1].split(".")[0]
        for f in files
    ]
    return sorted(set(orbits))

# ------------------------------------------------------------------
def get_orbit_filename(orbit_number=None,
                       product_type=DEFAULT_PRODUCT_TYPE,
                       baseline=DEFAULT_BASELINE):
    """
    Return a list of EarthCARE filenames for a given orbit number.
    
    Uses local geometa file if available (or creates it if missing).
    Falls back to remote lookup if needed.
    """
    if orbit_number is None:
        raise ValueError("Must specify an orbit_number.")

    orbit_str = f"{int(orbit_number):05d}"
    
    # Estimate year from orbit number
    (year, _), _, _ = get_orbit_date_approx(orbit_number)

    # Step 1: Ensure GEOMETA file exists and includes this orbit
    geometa_path = get_or_create_earthcare_geometa(
        orbit_number=orbit_number,
        year=year,
        product=product_type,
        baseline=baseline
    )

    # Step 2: Search matching lines in geometa
    with open(geometa_path) as f:
        matches = [line.strip() for line in f if f"_{orbit_str}" in line]

    if not matches:
        raise FileNotFoundError(f"No file found for orbit {orbit_number} in geometa: {geometa_path}")

    return matches
# ------------------------------------------------------------------
def get_orbit_date(orbit_number=None):
    """
    Return a sorted list of unique (year, DOY) for a given EarthCARE orbit number.
    Uses GEOMETA file, populates it if needed.
    """
    from csat2.EarthCARE.geometa import get_or_create_earthcare_geometa

    if orbit_number is None:
        raise ValueError("orbit_number is required")

    (est_year, _), _, _ = get_orbit_date_approx(orbit_number)
    orbit_str = f"{int(orbit_number):05d}"
    doy_set = set()

    for year in [est_year, est_year + 1]:
        try:
            path = get_or_create_earthcare_geometa(orbit_number=orbit_number, year=year)
            with open(path) as f:
                for line in f:
                    if f"_{orbit_str}" in line:
                        try:
                            timestamp = line.split("_")[5]
                            dt = datetime.strptime(timestamp, "%Y%m%dT%H%M%SZ")
                            doy_set.add((dt.year, dt.timetuple().tm_yday))
                        except Exception as e:
                            print(f"[WARNING] Skipping line due to parse error: {line.strip()} - {e}")
        except (FileNotFoundError, ValueError):
            continue

    if not doy_set:
        raise FileNotFoundError(f"No entry found for orbit {orbit_number} after updating GEOMETA.")

    return sorted(doy_set)





def filename_to_datetime(fname):
    """
    Extracts start datetime from EarthCARE filename.
    Works with .ZIP or .h5/.H5 files.
    """
    fname = fname.upper().replace(".ZIP", "").replace(".H5", "")  # normalize
    parts = fname.split("_")
    start_str = parts[5]  # should be like 20250320T231010Z
    return datetime.strptime(start_str, "%Y%m%dT%H%M%SZ")





# def get_orbit_datetimes(orbit_number=None,
#                         product_type=DEFAULT_PRODUCT_TYPE,
#                         baseline=DEFAULT_BASELINE):
#     """
#     Return a sorted list of unique start datetimes (as datetime objects)
#     for all EarthCARE files corresponding to the given orbit.
#     """
#     if orbit_number is None:
#         raise ValueError("Must specify an orbit_number.")

#     filenames = get_orbit_filename(
#         orbit_number=orbit_number,
#         product_type=product_type,
#         baseline=baseline
#     )

#     # Extract and deduplicate start datetimes
#     datetimes = {filename_to_datetime(f) for f in filenames}
#     return sorted(datetimes)



def get_orbit_datetimes(orbit_number,
                        orbit_id,
                        product_type=DEFAULT_PRODUCT_TYPE,
                        baseline=DEFAULT_BASELINE):
    """
    Return the start datetime (as a datetime object) for a specific
    EarthCARE granule identified by orbit number and orbit ID.
    """
    orbit_str = f"{int(orbit_number):05d}{orbit_id}"

    filenames = get_orbit_filename(
        orbit_number=orbit_number,
        product_type=product_type,
        baseline=baseline
    )

    matching = [f for f in filenames if f.endswith(f"{orbit_str}.ZIP") or f.endswith(f"{orbit_str}.h5")]

    if not matching:
        raise FileNotFoundError(f"No matching file found for orbit {orbit_str}")

    if len(matching) > 1:
        raise RuntimeError(f"Multiple files found for orbit {orbit_str}, expected only one.")

    return filename_to_datetime(matching[0])

# ------------------------------------------------------------------

def get_orbit_by_time(dtime: datetime,
                      product_type: str = DEFAULT_PRODUCT_TYPE,
                      baseline: str     = DEFAULT_BASELINE) -> str:
    """
    Given a datetime, return the EarthCARE orbit ID whose start time
    is closest. E.g.
    """
    # 1) estimate candidate orbits (±1 orbit)
    (lo, hi), _ = get_orbit_number_approx(dtime)

    best = {"delta": float("inf"), "filename": None}

    # 2) scan each candidate orbit
    for orbit in range(lo, hi + 1):
        geometa_path = get_or_create_earthcare_geometa(
            orbit_number=orbit,
            product=product_type,
            baseline=baseline
        )
        orbit_str = f"{orbit:05d}"
        with open(geometa_path, "r") as fh:
            files = [
                line.strip()
                for line in fh
                if line.strip() and f"_{orbit_str}" in line
            ]

        # 3) find closest‐in‐time granule
        for fname in files:
            start_dt = datetime.strptime(
                fname.split("_")[5].rstrip("Z"),
                "%Y%m%dT%H%M%S"
            )
            delta = abs((dtime - start_dt).total_seconds())
            if delta < best["delta"]:
                best.update(delta=delta, filename=fname)

    if best["filename"] is None:
        raise FileNotFoundError(f"No local EarthCARE files found near {dtime}.")

    # 4) extract and return the orbit ID (e.g. "04606A")
    orbit_id = best["filename"].rsplit("_", 1)[-1].replace(".ZIP", "")
    return orbit_id

# ------------------------------------------------------------------
def readin_earthcare_curtain(product_type=DEFAULT_PRODUCT_TYPE,
                              baseline=DEFAULT_BASELINE,
                              orbit_number=None,
                              orbit_id=None,
                              sds=None):
    """
    Reads one EarthCARE curtain HDF5 file by orbit number and orbit ID.
    """
    if orbit_number is None or orbit_id is None:
        raise ValueError("Both orbit_number and orbit_id must be specified.")

    orbit_str = f"{int(orbit_number):05d}"
    date_tuples = get_orbit_date(orbit_number)  # Could return 2 days

    found_file = None
    for year, doy in date_tuples:
        _, month, day = doy_to_date(year, doy)

        files = locator.search(
            "EARTHCARE", product_type,
            baseline=baseline,
            year=year, month=month, day=day,
            orbit=orbit_str,
            orbit_id=orbit_id,
            exit_first=False,
        )

        h5files = [f for f in files if f.endswith(".h5")]
        zipfiles = [f for f in files if f.lower().endswith(".zip")]

        if h5files:
            found_file = h5files[0]
            break  # we found it, no need to check more

        if zipfiles and not found_file:
            found_file = f"[ZIP FOUND] {zipfiles[0]}"  # Save ZIP message in case no .h5 is found


    if found_file is None:
        raise IOError("No .h5 files found for this orbit on any nearby day.")

    if found_file.startswith("[ZIP FOUND]"):
        raise IOError(f".h5 file not found, but .ZIP archive exists. Please unzip:\n{found_file[11:]}")

    print(f"[INFO] Reading file: {found_file}")

    if isinstance(sds, str):
        sds = [sds]

    return read_hdf5(found_file, sds=sds)

