from datetime import datetime, timedelta
from csat2 import locator
import xarray as xr
from csat2.ISCCP.download import download_files
import calendar
import os
from .utils import check_valid_args
import numpy as np


class Granule:
    '''An object for accessing and managing a specific ISCCP granule.
    Note that the lat is [-90, 90] and lon is [0, 360] in the granule files.'''

    def __init__(self, year: int, doy: int, time: int, product='isccp-basic', version='hgg'):
        
        '''Initialize a Granule object for a specific year, day-of-year, and time.'''

        check_valid_args(year, doy, time, product, version)
        self.year = year
        self.doy = doy
        self.time = time
        self.product = product
        self.version = version
        self.lonlat = None
        self.local_path = None

    
    def get_fileloc(self):
        '''Return the local directory and filename (location) for the granule.'''
        product_uppercase = {'isccp-basic': 'ISCCP-Basic', 'isccp': 'ISCCP'} ## converts the product name as it appears differently inthe filename and directory names
        version_uppercase = self.version.upper()  # Convert version to uppercase as it appears uppercase in the local filename
        #local_dir = locator.get_folder(product_uppercase[self.product],version_uppercase,year=self.year,doy = self.doy)

        time_str = f"{self.time:02d}" + '00'

        fileloc = locator.format_filename(product_uppercase[self.product], version_uppercase, year=self.year, doy=self.doy, time=time_str)

        return fileloc

    def check(self, download=True):
        '''Ensure the granule exists locally; if not, download it.'''

        fileloc = self.get_fileloc()
        self.local_path = fileloc  #cache it

        if not os.path.exists(fileloc):
            if download:
                self.local_path = download_files(
                    self.year, self.doy, self.time,
                    product=self.product,
                    version=self.version,
                    force_redownload=False
                )
            else:
                raise FileNotFoundError(f"Granule file not found: {fileloc}")

        return self.local_path

    def get_variable(self, varname: str):
        '''Read in a specific variable from the granule file.'''
        self.check()
        with xr.open_dataset(self.local_path) as ds:
            if varname not in ds:
                raise ValueError(f"Variable '{varname}' not found in dataset.")
            
            if self.lonlat is None and 'lat' in ds and 'lon' in ds:
                self.lonlat = (ds['lon'].values, ds['lat'].values)

            return ds[varname].load()
        
        
    def get_lonlat(self, grid=False):
        '''Return the longitude and latitude arrays from the granule file.
        This doesn't change with each timestep, so it is cached after the first call.
        If grid is True, return the lon/lat as 2D meshgrid arrays (lon_grid, lat_grid).'''

        # Return cached result if available
        if self.lonlat is not None:
            return self.lonlat

        # Load lon/lat from file
        self.check()
        with xr.open_dataset(self.local_path) as ds:
            if 'lon' not in ds or 'lat' not in ds:
                raise ValueError("Longitude and latitude coordinates not found in dataset.")

            lon = ds['lon'].values
            lat = ds['lat'].values

            if grid:
                lon, lat = np.meshgrid(lon, lat)

            self.lonlat = (lon, lat)

        return self.lonlat

    def get_metadata(self):
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

        self.check()
        ds = xr.open_dataset(self.local_path)

        # Filter attributes to only include target_metadata
        metadata = {
            'attributes': {k: v for k, v in ds.attrs.items() if k in target_metadata},
            'coordinates': {coord  for coord in ds.coords}
        }
        ds.close()
        return metadata    
        
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

        return Granule(next_year, next_doy, next_time, product=self.product, version=self.version)
    

    # def(self,lon,lat, varname=None, method='nearest'):
    #     '''function to read out variable and coloate given an arbitrary number of lon, lats'''
    #     )
   
    

