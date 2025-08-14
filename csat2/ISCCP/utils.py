import calendar
import numpy as np
import xarray as xr


def check_valid_gran_args( year: int, doy: int, time: int):
    '''Check that ISCCP granule arguments are valid. Raises ValueError if invalid.'''


    valid_times = range(0, 24, 3)

    max_doy = 366 if calendar.isleap(year) else 365

    min_year, max_year = 1983, 2017 ## TODO: make this more dynamic incase more data is added

    if not (min_year <= year <= max_year):
        raise ValueError(f"Year must be between {min_year} and {max_year}, got {year}.")  

    if not (1 <= doy <= max_doy):
        raise ValueError(f"Day-of-year (doy) must be between 1 and {max_doy} for year {year}, got {doy}.")
    
    if time not in valid_times:
        raise ValueError("Time must be in 3-hour UTC intervals (0, 3, 6, ..., 21).")

def geolocate(self, collection, product, varname, target_lons, target_lats, method='nearest'):
    '''Find colocated variable values for arbitrary longitude/latitude coordinates.
    
    Parameters:
    -----------
    collection : str
        Collection name (e.g., 'isccp-basic', 'isccp')
    product : str
        Product name (e.g., 'hgg', 'hgh', 'hgm')
    varname : str
        Variable name to extract
    target_lons : array-like
        Target longitude coordinates (can be single value or array)
    target_lats : array-like
        Target latitude coordinates (can be single value or array)
    method : str, optional
        Interpolation method ('nearest' or 'linear'), default is 'nearest'
    
    Returns:
    --------
    xarray.DataArray
        Colocated variable values at target coordinates
    '''
    
    # Ensure inputs are numpy arrays
    target_lons = np.atleast_1d(np.array(target_lons))
    target_lats = np.atleast_1d(np.array(target_lats))
    
    # Check that lon/lat arrays have same length
    if len(target_lons) != len(target_lats):
        raise ValueError("target_lons and target_lats must have the same length")
    
    # Normalize longitudes to [0, 360] to match ISCCP convention
    target_lons = target_lons % 360
    
    # Get the variable data
    var_data = self.get_variable(collection, product, varname)
    
    # Get grid coordinates
    grid_lons, grid_lats = self.get_lonlat(collection, product)
    
    # Create coordinate arrays for xarray interpolation
    if method == 'nearest':
        # Use xarray's sel with nearest method for exact grid point matching
        colocated_data = []
        
        for i, (target_lon, target_lat) in enumerate(zip(target_lons, target_lats)):
            try:
                # Find nearest grid point
                point_data = var_data.sel(
                    lon=target_lon, 
                    lat=target_lat, 
                    method='nearest'
                )
                colocated_data.append(point_data)
            except KeyError as e:
                # Handle case where coordinates are outside grid bounds
                print(f"Warning: Point ({target_lon:.2f}, {target_lat:.2f}) outside grid bounds")
                # Create a masked/NaN value with same structure as valid data
                nan_data = point_data * np.nan if 'point_data' in locals() else np.nan
                colocated_data.append(nan_data)
        
        # Combine results into single DataArray
        if len(colocated_data) > 1:
            # Create new coordinate for the points
            point_coord = np.arange(len(colocated_data))
            result = xr.concat(colocated_data, dim='point')
            result = result.assign_coords(point=point_coord)
            
            # Add coordinate metadata
            result = result.assign_coords(
                target_lon=('point', target_lons),
                target_lat=('point', target_lats)
            )
        else:
            result = colocated_data[0]
            # Add coordinate metadata for single point
            result = result.assign_coords(
                target_lon=target_lons[0],
                target_lat=target_lats[0]
            )
    
    elif method == 'linear':
        # Use xarray's interp for linear interpolation
        target_points = xr.Dataset({
            'lon': ('point', target_lons),
            'lat': ('point', target_lats)
        })
        
        result = var_data.interp(
            lon=target_points.lon,
            lat=target_points.lat,
            method='linear'
        )
        
        # Add coordinate metadata
        result = result.assign_coords(
            target_lon=('point', target_lons),
            target_lat=('point', target_lats)
        )
        
    else:
        raise ValueError(f"Method '{method}' not supported. Use 'nearest' or 'linear'")
    
    # Add metadata about the geolocation
    result.attrs.update({
        'geolocation_method': method,
        'target_coordinates': f'{len(target_lons)} point(s)',
        'source_granule': f'{self.year}-{self.doy:03d}T{self.time:02d}:00'
    })
    
    return result     


def check_valid_collection_product(collection:str,product:str):
    '''check that the products and collections are valid, some other product collection combinations haven't been implemented'''
        
    valid_collections = ['isccp', 'isccp-basic']

    valid_products = ['hgg', 'hgh', 'hgm', 'hgx']
            

    if collection not in valid_collections:
        raise ValueError(f"Invalid collection: '{collection}'. Must be one of {valid_collections}.")

    if product not in valid_products:
       raise ValueError(f"Invalid product: '{product}'. Must be one of {valid_products}.")

    if collection == 'isccp-basic' and product == 'hgx':
        raise ValueError("product 'hgx' is unavailable for ISCCP Basic.")


    

