from datetime import datetime, timedelta

# Default EarthCare product baseline. This is much less consistent than MODIS,
# so we might remove this and force users to specify it.
DEFAULT_BASELINE = "AE"

# same constants as before
REFERENCE_ORBIT   = 4451
REFERENCE_DATE    = datetime(2025, 3, 11, 0, 52, 6)     # orbit 4451 timestamp
ORBIT_DAY_PORTION = 0.0642650462962963                 # days per orbit ≈ 1/14.57


frame_names = ['A', 'B', 'C', 'D',
               'E', 'F', 'G', 'H']

# Names for the longitude/latitude variables in each file
lonlat_vars = {'ATL_NOM_1B': ["ellipsoid_longitude",
                              "ellipsoid_latitude"],
               'CPR_NOM_1B': ["longitude", "latitude"],
               'MSI_RGR_1C': ["longitude", "latitude"]}



def get_product_level(product):
    return product.split('_')[-1][0]

def get_orbit_date_approx(orbit_number=None):
    """
    Approximate an EarthCARE orbit’s calendar date.

    Reference
    ---------
    • Orbit 4451 was recorded at 2025-03-11 00:52:06 UTC
      (taken from filename “04451E”).
    • Average orbital rate ≈ 14.57 orbits day-¹
      (1 orbit ≈ 0.064265 day).
    """
    if orbit_number is None:
        raise ValueError("orbit_number is required")


    # Mean fraction of a day per orbit
    orbits_per_day    = 1 / ORBIT_DAY_PORTION

    # Offset in orbits → offset in days
    delta_orbits  = orbit_number - REFERENCE_ORBIT
    estimated_dt  = REFERENCE_DATE  + timedelta(days=delta_orbits / orbits_per_day)
    estimated_date = estimated_dt.date()

    year = estimated_date.year
    doy  = estimated_dt.timetuple().tm_yday

    note = ("Returns ((year, DOY), date). "
            "Accuracy is roughly ±1 day due to orbital timing drift.")

    return (year, doy), estimated_date, note


def get_orbit_number_approx(dt):
    """
    Given a datetime (UTC), approximate which EarthCARE orbit you're in.
    Returns (float_estimate, (min_orbit, max_orbit), note).
    """
    # 1) delta in days
    delta_days = (dt - REFERENCE_DATE).total_seconds() / 86400
    
    # 2) delta in orbits (can be fractional)
    delta_orbits = delta_days / ORBIT_DAY_PORTION
    
    # 3) absolute orbit number (float)
    est_orbit = REFERENCE_ORBIT + delta_orbits
    
    # 4) because each orbit spans 8 granules across ±~2 hours,
    #    and your average-period model drifts ~1 orbit/day,
    #    we’ll give ±1 orbit as a conservative uncertainty:
    lo = int(est_orbit) - 1
    hi = int(est_orbit) + 1
    
    note = ("Estimate is float; integer orbits will be around "
            f"{lo}–{hi} (±1) due to model drift.")
    return (lo, hi), note

