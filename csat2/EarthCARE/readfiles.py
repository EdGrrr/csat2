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
import netCDF4
import logging
import xarray as xr
import numpy as np
from csat2 import locator
import csat2.misc.hdf
from csat2.misc.time import doy_to_date
from csat2.EarthCARE.utils import DEFAULT_VERSION, get_orbit_date_approx, get_orbit_number_approx
from csat2.EarthCARE.geometa import get_or_create_earthcare_geometa
from csat2.EarthCARE.download import download_file_locations

log = logging.getLogger(__name__)


def files_available(product,
                    year, doy,
                    version=DEFAULT_VERSION):
    """
    Return a sorted list of local EarthCARE filenames for one day.
    """
    filenames = locator.search(
        "EarthCARE", product,
        year=year, doy=doy,
        version=version,
        orbit="*****", frame="*",
        exit_first=False,
    )

    if not filenames:
        print(f"No local EarthCARE files for {year:04d}-{month:02d}-{day:02d}")

    return sorted(filenames)


def variables_available(product, year, doy, orbit='*****', frame='*',
                        version=DEFAULT_VERSION):
    """
    Return a dict of all datasets in the first local EarthCARE HDF5 file
    found for the given date.

    Example
    -------
    >>> vars = variables_available(year=2025, doy=100)
    >>> list(vars)
    ['ScienceData/latitude', 'ScienceData/longitude', ...]
    """
    files = locator.search(
        "EarthCARE", product,
        year=year, doy=doy,
        version=version,
        orbit=orbit, frame=frame,
        exit_first=False,
    )
    return read_hdf5(files[0])


def available_orbits(product,
                     year, doy,
                     version=DEFAULT_VERSION):
    """
    List the orbit identifiers present on disk for one day.

    Example
    -------
    >>> available_orbits('ATL_NOM_1B', 2025, 100)
    ['04603H', '04604A', ...]
    """
    # Look for every matching file (ZIP or H5) on this date
    files = locator.search(
        "EarthCARE", product,
        year=year, doy=doy,
        version=version,
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


def get_orbit_filename(product, orbit,
                       version=DEFAULT_VERSION):
    """
    Return a list of EarthCARE filenames for a given orbit number.
    
    Uses local geometa file if available (or creates it if missing).
    Falls back to remote lookup if needed.
    """
    orbit_str = f"{int(orbit):05d}"
    
    # Estimate year from orbit number
    (year, _), _, _ = get_orbit_date_approx(orbit)

    # Step 1: Ensure GEOMETA file exists and includes this orbit
    geometa_path = get_or_create_earthcare_geometa(
        orbit=orbit,
        year=year,
        product=product,
        version=version
    )

    # Step 2: Search matching lines in geometa
    with open(geometa_path) as f:
        matches = [line.strip() for line in f if f"_{orbit_str}" in line]

    if not matches:
        raise FileNotFoundError(f"No file found for orbit {orbit_number} in geometa: {geometa_path}")

    return matches


def get_orbit_date(orbit):
    """
    Return a sorted list of unique (year, DOY) for a given EarthCARE orbit number.
    Uses GEOMETA file, populates it if needed.
    """
    from csat2.EarthCARE.geometa import get_or_create_earthcare_geometa

    if orbit_number is None:
        raise ValueError("orbit_number is required")

    (est_year, _), _, _ = get_orbit_date_approx(orbit)
    orbit_str = f"{int(orbit):05d}"
    doy_set = set()

    for year in [est_year, est_year + 1]:
        try:
            path = get_or_create_earthcare_geometa(year=year)
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
        raise FileNotFoundError(f"No entry found for orbit {orbit} after updating GEOMETA.")

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


def get_orbit_datetimes(orbit,
                        frame,
                        version=DEFAULT_VERSION):
    """
    Return the start datetime (as a datetime object) for a specific orbit/frame combination
    """
    orbit_str = f"{int(orbit):05d}{frame}"

    filenames = get_orbit_filename(
        product,
        orbit=orbit,
        version=version
    )

    matching = [f for f in filenames if f.endswith(f"{orbit_str}.ZIP") or f.endswith(f"{orbit_str}.h5")]

    if not matching:
        raise FileNotFoundError(f"No matching file found for orbit {orbit_str}")

    if len(matching) > 1:
        raise RuntimeError(f"Multiple files found for orbit {orbit_str}, expected only one.")

    return filename_to_datetime(matching[0])


def get_orbit_by_time(dtime: datetime,
                      product: str,
                      version: str     = DEFAULT_VERSION) -> str:
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
            version=version
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


def readin_earthcare_curtain(product,
                             orbit,
                             frame,
                             version=DEFAULT_VERSION,
                             sds=None):
    """
    Reads one EarthCARE curtain HDF5 file by orbit number and orbit ID.
    """
    orbit_str = f"{int(orbit):05d}"
    date_tuples = get_orbit_date(orbit)  # Could return 2 days

    found_file = None
    for year, doy in date_tuples:
        _, month, day = doy_to_date(year, doy)

        files = locator.search(
            "EarthCARE", product,
            year=year, doy=doy,
            orbit=orbit,
            frame=frame,
            version=version,
            exit_first=False,
        )

        h5files = [f for f in files if f.endswith(".h5")]

        if h5files:
            found_file = h5files[0]
            break  # we found it, no need to check more

    if found_file is None:
        raise IOError("No .h5 files found for this orbit on any nearby day.")

    if isinstance(sds, str):
        sds = [sds]

    return csat2.misc.hdf.read_hdf5(found_file, sds=sds)


def readin_earthcare_curtain_filename(filename,
                                      sds=None):
    """
    Reads one EarthCARE curtain HDF5 file by filename
    """
    if isinstance(sds, str):
        sds = [sds]

    log.debug(filename)
    with netCDF4.Dataset(filename) as ncdf:
        ds = xr.Dataset()
        tdims = []
        for name in sds:
            var = ncdf['ScienceData/'+name]
            var.set_auto_mask(False)
            dims = var.dimensions
            indata = var[:]
            try:
                indata = np.where(indata == var._Fillvalue, np.nan, indata)
            except AttributeError:  # No fill value
                pass
            ds[name] = xr.DataArray(indata, dims=dims)
            try:
                ds[name].attrs["units"] = var.units
            except AttributeError:
                pass
            tdims.extend(dims)
        tdims = set(tdims)
        for tdim in tdims:
            try:
                ds[tdim] = xr.DataArray(ncdf['ScienceData/'+tdim][:], dims=(tdim,))
            except (KeyError, IndexError):
                pass
        return ds
