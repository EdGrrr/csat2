from pathlib import Path
from datetime import timedelta, date
from csat2 import locator
from csat2.EarthCARE.utils import (
    get_orbit_date_approx,
    DEFAULT_BASELINE,
    DEFAULT_PRODUCT_TYPE,
    DEFAULT_YEAR,
)


# NOTE:
# For now, GEOMETA is stored as a plain text file containing one full
# EarthCARE filename per line.
#
# This format is kept for consistency with other EarthCARE modules,
# which already search GEOMETA as text and parse timestamps directly
# from full filename.
#
# A table-based GEOMETA format (for example orbit / year / DOY / hour)
# can be added later. For now, consistency with the rest of the codebase
# is the priority.


def _geometa_path(year):
    """Return the GEOMETA path for a given year and ensure its folder exists."""

    path = Path(locator.format_filename("EarthCARE", "GEOMETA", year=year))
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def read_geometa(year):
    """
    Return all non-empty GEOMETA entries for a given year.

    GEOMETA currently stores one full EarthCARE file ID per line.
    """
    path = _geometa_path(year)

    if not path.exists():
        return []

    with open(path, "r") as f:
        return [line.strip() for line in f if line.strip()]


def _append_unique_lines(path, entries):
    """
    Append unique non-empty lines to a GEOMETA file.

    The file stores one full EarthCARE file ID per line.
    """
    entries = sorted({entry.strip() for entry in entries if entry and entry.strip()})

    if not entries:
        return path

    if path.exists():
        with open(path, "r") as f:
            existing = {line.strip() for line in f if line.strip()}
    else:
        existing = set()

    new_entries = [entry for entry in entries if entry not in existing]

    if new_entries:
        mode = "a" if path.exists() else "w"
        with open(path, mode) as f:
            for entry in new_entries:
                f.write(entry + "\n")
    
    return path


def _list_remote_ids_for_day(year, doy,
                             product, baseline, limit= 1000,):
    """
    Return remote EarthCARE file IDs for one day as strings ending in '.ZIP'.
    """

    from csat2.EarthCARE.download import download_file_locations

    items = download_file_locations(
        product=product,
        year=year,
        doy=doy,
        baseline=baseline,
        limit=limit,
    )

    return sorted(
        f"{item['id']}.ZIP"
        for item in items
        if item.get("id")
    )



def get_or_create_earthcare_geometa(
    orbit,
    year = DEFAULT_YEAR,
    product = DEFAULT_PRODUCT_TYPE,
    baseline = DEFAULT_BASELINE,
):
    """
    Ensure the GEOMETA file for the requested year exists and contains entries
    for the requested orbit.
    """
    (_, _), date_est, _ = get_orbit_date_approx(orbit)

    # Keep the caller-provided year if given.
    # Only fall back to the estimated orbit year when year is None.
    if year is None:
        year = date_est.year

    path = _geometa_path(year)

    orbit_str = f"{orbit:05d}"
    entries = []

    for offset in (-1, 0, 1):
        d = date_est + timedelta(days=offset)

        # Only collect entries that belong to the GEOMETA file year being built.
        if d.year != year:
            continue

        try:
            ids = _list_remote_ids_for_day(
                year=d.year,
                doy=d.timetuple().tm_yday,
                product=product,
                baseline=baseline,
            )

        except Exception as e:
            print(f"[WARNING] Could not list files for {d}: {e}")
            continue

        for file_id in ids:
            if f"_{orbit_str}" in file_id:
                entries.append(file_id)

    if not entries:
        raise FileNotFoundError(f"No files found for orbit {orbit} in year {year} near {date_est}")

    return _append_unique_lines(path, entries)



def fill_geometa_for_day(year, doy,
                         product = DEFAULT_PRODUCT_TYPE,
                         baseline = DEFAULT_BASELINE):
    """
    Ensure the GEOMETA file for the given year contains all EarthCARE file IDs
    found for the requested day.
    """
    path = _geometa_path(year)

    try:
        ids = _list_remote_ids_for_day(
            year=year,
            doy=doy,
            product=product,
            baseline=baseline,
        )
    except Exception as e:
        print(f"[ERROR] Could not list files for {year}-{doy:03d}: {e}")
        return None

    if not ids:
        return path

    return _append_unique_lines(path, ids)


def create_geometa(
    year,
    product = DEFAULT_PRODUCT_TYPE,
    baseline = DEFAULT_BASELINE,
):
    """
    Populate the yearly GEOMETA file by scanning day by day.

    For the current year, only scan up to today's date.
    For future years, return the GEOMETA path without scanning.
    """
    today = date.today()

    start = date(year, 1, 1)
    end = date(year, 12, 31)

    if year > today.year:
        return _geometa_path(year)

    if year == today.year:
        end = today

    current = start
    while current <= end:
        fill_geometa_for_day(
            year=current.year,
            doy=current.timetuple().tm_yday,
            product=product,
            baseline=baseline,
        )
        current += timedelta(days=1)

    return _geometa_path(year)


