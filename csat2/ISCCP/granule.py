from datetime import datetime, timedelta
from csat2 import locator
import xarray as xr
from csat2.ISCCP.download import download_files
import calendar
import os
from .utils import check_valid_args


class Granule:
    '''An object for accessing and managing a specific ISCCP granule.'''

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

    def check(self,download = True):
        '''Ensure the granule exists locally; if not, download it.'''
        product_uppercase = {'isccp-basic': 'ISCCP-Basic', 'isccp': 'ISCCP'} ## converts the product name as it appears in two different ways, uppercase version in the remote filename and the local file naming convention
        version_uppercase = self.version.upper()  # Convert version to uppercase as it appears uppercase in the remote filename
        local_dir = locator.get_folder(product_uppercase[self.product],version_uppercase,year=self.year,doy = self.doy)

        time_str = f"{self.time:02d}" + '00'
        local_filename = locator.format_filename(product_uppercase[self.product], version_uppercase, year=self.year, doy=self.doy, time=time_str)
        self.local_path = os.path.join(local_dir, local_filename)

        # Download file if not found and download is True
        if not os.path.exists(self.local_path):
            if download:
                self.local_path = download_files(
                    self.year, self.doy, self.time,
                    product=self.product,
                    version=self.version,
                    force_redownload=False
                )
            else:
                raise FileNotFoundError(f"Granule file not found: {self.local_path}")

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

    def get_metadata(self):
        '''Load and return metadata (global attributes and coordinate info).'''
        self.check()
        ds = xr.open_dataset(self.local_path)

        metadata = {
            'attributes': dict(ds.attrs),
            'coordinates': {coord: ds[coord].values for coord in ds.coords}
        }

        # Cache lon/lat if available
        if self.lonlat is None and 'lat' in ds and 'lon' in ds:
            self.lonlat = (ds['lon'].values, ds['lat'].values)

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
    

