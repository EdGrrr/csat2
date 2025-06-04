import calendar


def check_valid_args( year: int, doy: int, time: int,product: str, version: str):
    '''Check that ISCCP granule arguments are valid. Raises ValueError if invalid.'''

    valid_products = ['isccp', 'isccp-basic']
    valid_versions = ['hgg', 'hgh', 'hgm', 'hgx']
    valid_times = range(0, 24, 3)

    max_doy = 366 if calendar.isleap(year) else 365

    min_year, max_year = 1983, 2017 ## TODO: make this more dynamic incase more data is added

    if not (min_year <= year <= max_year):
        raise ValueError(f"Year must be between {min_year} and {max_year}, got {year}.")  

    if not (1 <= doy <= max_doy):
        raise ValueError(f"Day-of-year (doy) must be between 1 and {max_doy} for year {year}, got {doy}.")
    
    if time not in valid_times:
        raise ValueError("Time must be in 3-hour UTC intervals (0, 3, 6, ..., 21).")

    if product not in valid_products:
        raise ValueError(f"Invalid product: '{product}'. Must be one of {valid_products}.")

    if version not in valid_versions:
        raise ValueError(f"Invalid version: '{version}'. Must be one of {valid_versions}.")

    if product == 'isccp-basic' and version == 'hgx':
        raise ValueError("Version 'hgx' is unavailable for ISCCP Basic.")


    

