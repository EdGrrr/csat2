import os
from csat2 import locator
from csat2 import misc
from csat2.CERES.download import download_files
import netCDF4 as nc
import xarray as xr
import numpy as np



class Granule:
    '''An object for accessing and managing a specific CERES granule
    Note that this is only currently set up to deal with the CERES SYN - hourly product, and may have to be substancially editted if
    you want the data at different temporal frequencies.'''
    def __init__(self,year:int, doy :int, time:int):
        self.year = year
        _,month,dom = misc.time.doy_to_date(year, doy)  # Get month and day of month (dom) from year and day of year
        self.month = month
        self.dom = dom
        self.doy = doy
        self.time = time ## this is the UTC time in hours (0-23)
        

        
    def get_fileloc(self):
        ''' Return the local directory, note that this does not require that the file actually exists there.'''
        
        fileloc = locator.format_filename('CERES','SYN-hourly',year = self.year, month = self.month, dom  = self.dom)

        return fileloc
  

    def check(self):
        ''' Check if the file exists locally.'''
        fileloc = self.get_fileloc()
        if os.path.exists(fileloc):
            return True
        else:
            return False
        
    def download(self, force_redownload: bool = False):
            '''Download the file locally'''
            fileloc = self.get_fileloc()
            if not self.check() or force_redownload:
                download_files(
                    year=self.year,
                    month=self.month,  # Now passes month instead of doy
                    dom=self.dom,
                    local_path=fileloc
                )
                return f'Downloaded file to {fileloc}'
            else:
                return f'File already exists at {fileloc}, not redownloading.'
        


    def get_variable(self, varname: str):
        """
        Get a specific variable from the granule using xarray.
        Will download the file if it does not already exist locally.
        """
        # Ensure the file is downloaded
        self.download()
        
        fileloc = self.get_fileloc()
        
        # Check that the file exists
        if not os.path.exists(fileloc):
            raise FileNotFoundError(f"Granule file does not exist at {fileloc}")
        
        # Open the dataset with xarray (lazy-loading by default)
        ds = xr.open_dataset(fileloc, engine='netcdf4')  # or engine='h5netcdf' if needed
        
        if varname not in ds:
            raise KeyError(f"Variable '{varname}' not found in file {fileloc}")
        
        variable_data = ds[varname]  # This is an xarray.DataArray
        
        return variable_data
    
    def list_variables(self):
        """
        List all available variables in the granule.
        Will download the file if it does not already exist locally.
        """
        # Ensure the file is downloaded
        self.download()
        
        fileloc = self.get_fileloc()
        
        # Check that the file exists
        if not os.path.exists(fileloc):
            raise FileNotFoundError(f"Granule file does not exist at {fileloc}")
        
        # Open the NetCDF/HDF file and list the variables
        with nc.Dataset(fileloc, 'r') as dataset:
            variable_names = list(dataset.variables.keys())
        
        return variable_names
    
    def get_lonlat(self,grid:bool = False,lon_360:bool = False):
        ''' Return the longitude and latitude arrays from the granule, optional arguments to return as a grid or have longitudes in [0,360] range'''
        lon = self.get_variable('longitude').values
        lat = self.get_variable('latitude').values

        if grid:
            lon,lat = np.meshgrid(lon,lat)

        if lon_360:
            lon = np.where(lon <0, lon + 360, lon)
        return lon, lat
    
    def next(self):
        ''' Return a Granule object for the next hour'''
        new_time = self.time + 1
        new_doy = self.doy
        new_year = self.year
        
        if new_time >= 24:
            new_time = 0
            new_doy += 1
            # Check for year rollover
            if new_doy > misc.time.days_in_year(new_year):
                new_doy = 1
                new_year += 1
                
        return Granule(year=new_year, doy=new_doy, time=new_time)
    
    def geolocate(self, varname: str, target_lons, target_lats, method: str = "linear"):
            """
            Given target longitude and latitude, return the value of `varname`
            at the nearest CERES grid cell.

            Parameters
            ----------
            varname : str
                Name of CERES variable.
            target_lons, target_lats : float or array-like
                Coordinates to sample at.
            method : 'linear' (this is defaul and will be fastest)
                    'xarray' (this is jsut incase there are issues with linear, but will be slower, and not really needed for CERES SYN data, but included incase the grid happens to be non regular in future data products)

            Returns
            -------
            np.ndarray or float
                Sampled values at requested coordinates.
            """

            # Make arrays
            target_lons = np.asarray(target_lons)
            target_lats = np.asarray(target_lats)

            

            if method == "linear":
                
                return self._geolocate_linear(varname, target_lons, target_lats)
            

            elif method == "xarray":
                return self._geolocate_xarray(varname, target_lons, target_lats)

            else:
                raise NotImplementedError("method must be 'linear' or 'xarray'")
            
    def _geolocate_linear(self, varname, target_lons, target_lats):
            """
            Nearest-neighbor lookup using linear index mapping.
            Works because CERES SYN is on a regular lat/lon grid.
            The lon and lat have the following conventions:
            -lon: [-179.5, -178.5, ..., 178.5, 179.5]
            -lat: [89.5, 88.5, ..., -88.5, -89.5]
            - IMPORTANT: the data is stored as [time,lat, lon] in the file.
            ALSO IMPORTANT: This works best if the data is 2D (i.e., single time slice), but can handle higher dimensions
            VERY IMPORTANT: THIS ASSUMES THAT THE TARGET LONS ARE IN THE SAME CONVENTION AS THE CERES DATA (I.E., -179.5 TO 179.5).
            IF YOU ARE UNSURE USE THE XARRAY METHOD AS THIS IS PROBABLY MORE VERSATILE.
            """

            # Load grid
            lon,lat = self.get_lonlat() ## this test slows things down slightly but is important to ensure the method is valid
            assert lon[0] == -179.5 and lon[-1] == 179.5, "Longitude grid does not match expected CERES SYN grid. Use the 'xarray' method instead."
            assert lat[0] == 89.5 and lat[-1] == -89.5, "Latitude grid does not match expected CERES SYN grid. Use the 'xarray' method instead."
                
            target_lons_180 = ((target_lons + 180) % 360) - 180  # Convert to [-180, 180] range    

            lon_ind = np.round (target_lons_180 + 179.5).astype(int)
            lat_ind = np.round (89.5 - target_lats).astype(int)

            lon_ind = np.mod(lon_ind, lon.size)  # Wrap around longitudes
            lat_ind = np.clip(lat_ind, 0, lat.size - 1)  # Clip latitudes to valid range, no wrapping at the poles

            # Load variable data
            var = self.get_variable(varname).values
            assert var.ndim >=2, "Variable data must be at least 2D (lat, lon)."
            assert var.shape[-2] == lat.shape[0] and var.shape[-1] == lon.shape[0], \
                    "Variable data shape does not match expected CERES SYN grid shape."

            # Return values
            if var.ndim == 3:
                return np.array([var[:, li, lj] for li, lj in zip(lat_ind, lon_ind)])
            else:
                return var[lat_ind, lon_ind]

    def _geolocate_xarray(self, varname, target_lons, target_lats):
            """
            Nearest-neighbor lookup using xarray .sel(method='nearest').
            """

            ds = xr.open_dataset(self.get_fileloc(), engine="netcdf4")
            da = ds[varname]

            # Fix target longitude convention
            if da.lon.max() > 180:
                target_lons = np.mod(target_lons, 360)
            else:
                target_lons = ((target_lons + 180) % 360) - 180

            # Vectorized selection using DataArray indexers
            return da.sel(
                lon=xr.DataArray(target_lons, dims="points"),
                lat=xr.DataArray(target_lats, dims="points"),
                method="nearest"
            ).values
            
        