from pathlib import Path
from datetime import timedelta
from csat2 import locator
from csat2.EarthCARE.utils import get_orbit_date_approx

#from csat2.EarthCARE.download import download_file_locations
from csat2.EarthCARE.utils import get_orbit_date_approx, DEFAULT_BASELINE, DEFAULT_PRODUCT_TYPE, DEFAULT_YEAR

def get_or_create_earthcare_geometa(
    orbit_number=None,
    year=DEFAULT_YEAR,
    product=DEFAULT_PRODUCT_TYPE,
    baseline=DEFAULT_BASELINE
):

    """
    Ensures the EarthCARE GEOMETA file for a given year exists.
    If orbit_number is provided, ensures that orbit is listed.

    Defaults:
    - year = 2025
    - product = CPR_CLD_2A
    - baseline = AB
    """

    # Step 1: Estimate orbit date (if orbit provided)
    if orbit_number is not None:
        (_, _), date_est, _ = get_orbit_date_approx(orbit_number)
        year = date_est.year  # override year to match actual orbit

    # Step 2: Construct GEOMETA path based on final year
    path = Path(locator.format_filename("EARTHCARE", "GEOMETA", year=year))

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
        d = date_est + timedelta(days=offset)
        try:
            files = download_file_locations(product_type=product, baseline=baseline,
                                            year=d.year, month=d.month, day=d.day)
        except Exception as e:
            print(f"[WARNING] Could not list files for {d.date()}: {e}")
            continue

        for f in files:
            if f"_{orbit_str}" in f:
                entries.append(Path(f).name)

    if not entries:
        raise FileNotFoundError(f"No files found for orbit {orbit_number} near {date_est}")

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

