from datetime import datetime, timedelta
from csat2 import locator
from csat2 import misc
import xarray as xr
from csat2.ISCCP.download import download_files
import calendar
import os
from .utils import check_valid_gran_args,check_valid_collection_product
import numpy as np
import netCDF4


class Granule:
    '''An object for accessing and managing a specific ISCCP granule.
    Note that the lat is [-90, 90] and lon is [0, 360] in the granule files.'''

    def __init__(self, year: int, doy: int, time: int):
        
        '''Initialize a Granule object for a specific year, day-of-year, and time.
        The granule object is independant of hte product and type'''

        check_valid_gran_args(year, doy, time)
        month = misc.time.doy_to_date(year, doy)[1]  # Get month from day of year

        self.year = year
        self.month = month ## month is an integer from 1 to 12
        self.doy = doy
        self.time = time
        self.lonlat = None # this is okay to cache since the grid is not changing

    
    def get_fileloc(self,collection,product):
        '''Return the local directory and filename (location) for the granule, given a certain product e.g. (hgg,hgh,hgs,hgx) and collection (this referes to either isccp-basic or isccp, isccp-basic is just a subset of isccp containing useful variables).
         Note that this does not require the file to be downloaded'''

        check_valid_collection_product(collection,product)
        ## converts the product name as it appears differently inthe filename and directory names
        product_uppercase = product.upper()

        # Convert collection to uppercase as it appears uppercase in the local filename
        collection_uppercase = {'isccp-basic': 'ISCCP-Basic', 'isccp': 'ISCCP'} 
        
        time_str = f"{self.time:02d}" + '00'
        if product == 'hgg':
            fileloc = locator.format_filename(collection_uppercase[collection], product_uppercase, year=self.year, doy=self.doy, time=time_str)
        elif product == 'hgh':
            fileloc = locator.format_filename(collection_uppercase[collection], product_uppercase, year=self.year, month=self.month, time=time_str)
        elif product == 'hgm':
            fileloc = locator.format_filename(collection_uppercase[collection], product_uppercase, year=self.year, month = self.month)


        return fileloc

    

    def check(self,collection,product):
        '''Check if a specific product/ collection for the granule exists locally, returns True if the file exists'''


        check_valid_collection_product(collection,product)
        fileloc = self.get_fileloc(collection,product)

        if os.path.exists(fileloc):
            return True
        else:
            return False


    def download(self,collection,product,force_redownload = False):
        '''method for downloading a granule for a given product and collection'''

        if (not self.check(collection,product)) or force_redownload: # if the file does not already exist, or if force redownload is True, continue to download the file

            download_files(
                self.year, self.doy, self.time,
                collection=collection,
                product=product,
                force_redownload=True
            )
            return f'file downloaded to {self.get_fileloc(collection,product)}"'
        else:
            return f"file already exists at: {self.get_fileloc(collection, product)}"

    
    def get_variable(self,collection,product, varname: str):
        '''Read a specific variable from the granule file using netCDF4 and convert to xarray.
        netCDF4 handles scale/offset automatically: raw * scale + offset.
        '''
        self.check(collection,product) ## check if the file exists and download it if not

        with netCDF4.Dataset(self.get_fileloc(collection,product), mode='r') as nc:
            if varname not in nc.variables:
                raise ValueError(f"Variable '{varname}' not found in dataset.")
            
            var = nc.variables[varname]
            data = var[:] # this is a masked array object with fill values attribute .fill_value
            
            # Handle fill value
            fill_value = data.fill_value
            if fill_value is not None:
                #print(f"Replacing fill value {fill_value} with NaN in variable '{varname}'")
                data = np.ma.masked_equal(data, fill_value)


            dims = var.dimensions
            attrs = {attr: getattr(var, attr) for attr in var.ncattrs()}

            # Get coordinate variables if present
            coords = {}
            for dim in dims:
                if dim in nc.variables:
                    coords[dim] = nc.variables[dim][:]

            # Cache lon/lat if present and not yet set
            if self.lonlat is None and 'lon' in nc.variables and 'lat' in nc.variables:
                self.lonlat = (nc.variables['lon'][:], nc.variables['lat'][:])

        # Convert to xarray DataArray (no fill_value argument needed)
        da = xr.DataArray(data, dims=dims, coords=coords, attrs=attrs, name=varname)
        return da
            
        
    def get_lonlat(self,collection,product, grid=False):
        '''Return the longitude and latitude arrays from the granule file.
        This doesn't change with each timestep, so it is cached after the first call.
        If grid is True, return the lon/lat as 2D meshgrid arrays (lon_grid, lat_grid).'''

        # Return cached result if available
        if self.lonlat is not None:
            return self.lonlat

        # Load lon/lat from file
        self.check(collection,product)
        with xr.open_dataset(self.get_fileloc(collection,product)) as ds:
            if 'lon' not in ds or 'lat' not in ds:
                raise ValueError("Longitude and latitude coordinates not found in dataset.")

            lon = ds['lon'].values
            lat = ds['lat'].values

            if grid:
                lon, lat = np.meshgrid(lon, lat)

            self.lonlat = (lon, lat)

        return self.lonlat

    def get_metadata(self,collection,product):
        '''Load and return subset of metadata (adapt based off what is of interest).'''

        target_metadata = [
            'summary',
            'project',
            'time_coverage_start',
            'time_coverage_end',
            'geospatial_bounds',
            'isccp_number_of_satellites_contributing',
            'platform',
            'instrument',
            'isccp_percent_full_cells',
            'geospatial_lat_resolution',
            'geospatial_lon_resolution'
        ]

        self.check(collection,product)
        ds = xr.open_dataset(self.get_fileloc(collection,product))

        # Filter attributes to only include target_metadata
        metadata = {
            'attributes': {k: v for k, v in ds.attrs.items() if k in target_metadata},
            'coordinates': {coord  for coord in ds.coords}
        }
        ds.close()
        return metadata
        
    def geolocate(self, collection, product, varname: str, lon, lat, method='gridded'):
        '''For a specific granule and set of lon and lats, returns the corresponding ISCCP data as an xarray object
        Note, that is only accepts longitudes in the range [0, 360] and assumes that the ISCCP data is on a regular 1 x 1 degree grid.
        i.e. lons [0.5, 1.5, ..., 359.5] and lats [-89.5, -88.5, ..., 89.5]

        Parameters:
        -----------
        lon : array like or scalar, can convert between different longitude conventions
        lat : array like, this needs to be between [-90,90]
        collection : str
            Collection name
        product : str
            Product name
        varname : Variable of interest, note it returns an xarray object with the same structure as the target
        method : gridded, this is the only method currently implemented, this will fail for some of the ISCCP data that is on a differnt grid
        Returns:
        --------
        data : xarray object
            ISCCP data at the specified locations, maintaining input array shape
        '''
        
        if method == 'gridded':
            # Check that ISCCP data has expected grid otherwise the lon/ lat to index conversion will be invalid
            ISCCP_lon, ISCCP_lat = self.get_lonlat(collection, product)
            
            expected_lon = np.linspace(0.5, 359.5, 360)
            expected_lat = np.linspace(-89.5, 89.5, 180)
            
            if not (np.allclose(ISCCP_lon, expected_lon) and np.allclose(ISCCP_lat, expected_lat)):
                raise ValueError("ISCCP longitude/latitude grids are not in the expected format.")
            
            # Convert inputs to numpy arrays
            lon = np.asarray(lon)
            lat = np.asarray(lat)
            
            # Check that lon and lat have the same shape
            if lon.shape != lat.shape:
                raise ValueError("Longitude and latitude arrays must have the same shape.")
            
            # Store original shape for output
            original_shape = lon.shape
            
            # Flatten arrays for processing
            lon_flat = lon.flatten()
            lat_flat = lat.flatten()

            ## ensure that the lons are in the [0,360) range, if not wrap them around such that they are
            lon_flat = np.asarray(lon_flat, dtype=float) % 360  ## modulo operator

            # Validate bounds, this is slighly superfluous for lons, but is nessecary for the lats
            if np.any((lon_flat < 0) | (lon_flat > 360)):
                raise ValueError("Longitude must be in [0, 360].")
            if np.any((lat_flat < -90) | (lat_flat > 90)):
                raise ValueError("Latitude must be in [-90, 90].")
            
            # Calculate indices
            lon_indices = np.floor(lon_flat).astype(int) ## any lon between [0,360) will map to [0,359] index
            lat_indices = np.floor(lat_flat + 90).astype(int) 
            
            # ensure the indices fall within  valid range, no indices should fall outside of this range however it is abit of a safety check to stop things breaking
            lon_indices = np.clip(lon_indices, 0, 359)
            lat_indices = np.clip(lat_indices, 0, 179)
            
            # Get the variable data
            var_data = self.get_variable(collection, product, varname)
            

            
            # Find lat and lon dimension positions
            dims = var_data.dims
            try:
                lat_dim_idx = dims.index('lat')
                lon_dim_idx = dims.index('lon')
            except ValueError:
                raise ValueError(f"Could not find 'lat' and 'lon' dimensions in variable '{varname}'")
            
            print(f"Lat dimension at index {lat_dim_idx}, Lon dimension at index {lon_dim_idx}")
            
            # Convert to numpy array for indexingk
            if hasattr(var_data, 'values'):
                data_array = var_data.values
            else:
                data_array = np.array(var_data)
            
            # Create advanced indexing list
            # We need to handle the case where lat/lon are not the last dimensions
            indexing_list = [slice(None)] * data_array.ndim ## generates a list of slices for all dimensions
            
            # Extract data at specified locations
            data_list = []
            for lat_idx, lon_idx in zip(lat_indices, lon_indices):
                # Create indexing list for this specific lat/lon pair
                current_indexing = indexing_list.copy()
                current_indexing[lat_dim_idx] = lat_idx # work out which are the lat and lon dimensions
                current_indexing[lon_dim_idx] = lon_idx
                
                # Extract data - this will preserve all other dimensions
                extracted_data = data_array[tuple(current_indexing)]
                data_list.append(extracted_data)
            
            # Convert to numpy array - this will have shape (n_points, *other_dims)
            if len(data_list) > 0:
                data_flat = np.array(data_list)
            else:
                data_flat = np.array([])
            
            # Reshape to original input shape
            if original_shape == ():
                # Handle scalar case
                if data_flat.ndim > 1:
                    # If there are other dimensions (like cloud type), keep them
                    data_values = data_flat[0]  # Take the single point
                else:
                    data_values = data_flat.item() if data_flat.size == 1 else data_flat
                
                # Create scalar DataArray with coordinates
                data = xr.DataArray(
                    data_values,
                    coords={'lon': lon.item(), 'lat': lat.item()},
                    attrs=var_data.attrs,
                    name=varname
                )
            else:
                # For array inputs e.g. a grid reshape the output to match the original structure
                # for flat data we have (n_points, *other_dims)
                # We want (other_dims..., *original_shape)
                
                if data_flat.ndim > 1:
                    # Move the points dimension to the end and reshape it to original_shape
                    # data_flat is (n_points, dim1, dim2, ...) -> (dim1, dim2, ..., *original_shape)
                    other_dims_shape = data_flat.shape[1:]  # Everything except the first (points) dimension
                    new_shape = other_dims_shape + original_shape
                    
                    # Transpose to move points dimension to the end, then reshape
                    axes = list(range(1, data_flat.ndim)) + [0]  # Move first dim to end
                    data_transposed = np.transpose(data_flat, axes)
                    data_values = data_transposed.reshape(new_shape)
                else:
                    data_values = data_flat.reshape(original_shape)
                
                # Create coordinate arrays in the same shape as the spatial dimensions
                lon_coords = np.broadcast_to(lon, original_shape)
                lat_coords = np.broadcast_to(lat, original_shape)
                
                # Get original non-spatial dimensions
                original_dims = list(var_data.dims)
                original_dims.remove('lat')
                original_dims.remove('lon')
                
                # Determine final dimension names
                if len(original_shape) == 1:
                    spatial_dims = ['points']
                    coords = {
                        'lon': ('points', lon_coords),
                        'lat': ('points', lat_coords)
                    }
                elif len(original_shape) == 2:
                    spatial_dims = ['y', 'x']
                    coords = {
                        'lon': (['y', 'x'], lon_coords),
                        'lat': (['y', 'x'], lat_coords)
                    }
                else:
                    # For higher dimensions use generic names
                    spatial_dims = [f'spatial_dim_{i}' for i in range(len(original_shape))]
                    coords = {
                        'lon': (spatial_dims, lon_coords),
                        'lat': (spatial_dims, lat_coords)
                    }
                
                final_dims = original_dims + spatial_dims
                
                # Create the data araray
                data = xr.DataArray(
                    data_values,
                    dims=final_dims,
                    coords=coords,
                    attrs=var_data.attrs,
                    name=varname
                )
                
        else:
            raise NotImplementedError(f"Method '{method}' is not implemented.")
        
        return data
        

    def next(self):
        '''Return a new Granule object for the next 3-hour interval, note that this is specific to hgg data at the minute.'''
        next_time = self.time + 3
        next_doy = self.doy
        next_year = self.year

        if next_time >= 24:
            next_time = 0
            next_doy += 1
            if next_doy > (366 if calendar.isleap(next_year) else 365):
                next_doy = 1
                next_year += 1

        return Granule(next_year, next_doy, next_time)
        


    

