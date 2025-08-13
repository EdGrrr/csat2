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

            self.local_path = download_files(
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
        with xr.open_dataset(self.local_path) as ds:
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
        ds = xr.open_dataset(self.local_path)

        # Filter attributes to only include target_metadata
        metadata = {
            'attributes': {k: v for k, v in ds.attrs.items() if k in target_metadata},
            'coordinates': {coord  for coord in ds.coords}
        }
        ds.close()
        return metadata
    
    def geolocate(self):
        ''''''
    
        
    def next(self):
        '''Return a new Granule object for the next 3-hour interval.'''
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
    


    

