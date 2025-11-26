import os
from csat2 import locator
from csat2 import misc
from csat2.CERES.download import download_files
import netCDF4 as nc
import xarray as xr
import numpy as np

#### SOME NOTES ON THE CERES SYN1deg-1Hour PRODUCT
# The CERES SYN1deg-1Hour product is a gridded product on a 1 degree latitude by 1 degree longitude grid, with hourly temporal resolution.
# The data files are in HDF format, but can be read using netCDF4 or xarray in Python.

# Ncld is the cloud layer index and takes the values 1 - 5 (or 0 to 4 in python/ numpy indexing), with 1: high, 2: Upper mid, 3: lower mid, 4: low, 5: Total
# 1 = High (50-300 mb), 2 = UpperMid (300-500 mb), 3 = LowerMid (500-700 mb), 4 = Low (700 mb-Surface), 5 = Total (50 mb - Surface)
# Nhour is the hour of day index (0-23 UTC) the first (or zeroth) index is hour 0 and the final index is hour 23
# Nlat has length 180 Index #1 is defined at 89.5°N and #180 is at 89.5°S (-89.5) 
# Nlon has length 360 Index #1 is defined at 179.5°W (-179.5) and #360 is at 179.5°E
# Common variables useful for liquid water cloud analysis are:
# [sza] - Solar zenith angle in degrees (0-90) shape (Nhour,Nlat,Nlon)
# [obs_cld_amount] - Observed cloud amount (or CF) expressed as a percentage (0-100) it has shape (Ncld,Nhour,Nlat,Nlon)
# [obs_cld_od] - Observed cloud optical depth ( 3.7 um retrieval - unitless, range: 0 - 400) shape (Ncld,Nhour,Nlat,Nlon)
# [obs_cld_lwp] - Observed cloud liquid water path (in g/m^2) shape (Ncld,Nhour,Nlat,Nlon)
# [obs_cld_liq_radius] - Observed cloud liquid water radius (in micrometers, range 0-40) shape (Ncld,Nhour,Nlat,Nlon)






class Granule:
    '''An object for accessing and managing a specific CERES granule
    Note that this is only currently set up to deal with the CERES SYN - hourly product, and may have to be substancially editted if
    you want the data at different temporal frequencies.'''
    def __init__(self,year:int, doy :int, time:int= None):
        self.year = year
        _,month,dom = misc.time.doy_to_date(year, doy)  # Get month and day of month (dom) from year and day of year
        self.month = month
        self.dom = dom
        self.doy = doy
        self.time = time ## this is the UTC time in hours (0-23), if no time is given we will just use the full daily file
        

        
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
        


    def get_variable(self, varname: str, daily = False):
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
        
        if self.time is None or daily:
            variable_data = ds[varname]  # Return the full variable if no time is specified
        else:
        
            variable_data = ds[varname].isel(gmt_hr_index=self.time) # Extract the variable at the specified time
        
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
        lon = self.get_variable('longitude', daily=True).values ## since a single time slice is not really defined for lon and lat
        lat = self.get_variable('latitude', daily=True).values

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

            lon_ind = np.floor(target_lons_180 + 180).astype(int) % lon.size
            lat_ind = np.floor(89.5 - target_lats).astype(int)
            lat_ind = np.clip(lat_ind, 0, lat.size - 1)

            # Load variable data
            var = self.get_variable(varname)

            var_data = var.values
            assert var_data.shape[-2:] == (lat.size, lon.size), "Variable grid mismatch."
            assert var_data.ndim >=2, "Variable data must be at least 2D (lat, lon)."

            original_shape = target_lons.shape ## save the oroginal shape for reshaping later
            assert target_lats.shape == original_shape, "Target lon and lat must have the same shape."


            # Flatten the lon and lat
            lon_flat = target_lons.ravel()
            lat_flat = target_lats.ravel()

            # Convert lon convention
            lon_flat_180 = ((lon_flat + 180) % 360) - 180

            # Compute exact grid indices
            lon_ind = np.round(lon_flat_180 + 179.5).astype(int)
            lat_ind = np.round(89.5 - lat_flat).astype(int)

            # Wrap and clip
            lon_ind = np.mod(lon_ind, lon.size)
            lat_ind = np.clip(lat_ind, 0, lat.size - 1)

            # Load CERES variable
            var = self.get_variable(varname)
            var_data = var.values

            # indexing based on variable dimensions, remeber that the var data is stored as [Plevel,Time,lat, lon] in the file

            if var_data.ndim == 2:      # (lat,lon)
                out = var_data[lat_ind, lon_ind]

                new_shape = original_shape
                new_dims = ["dim_" + str(i) for i in range(len(original_shape))]

            elif var_data.ndim == 3:    # (time,lat,lon)
                out = var_data[:, lat_ind, lon_ind]

                new_shape = (var_data.shape[0],) + original_shape
                new_dims = [var.dims[0]] + ["dim_" + str(i) for i in range(len(original_shape))]

            elif var_data.ndim == 4:    # (Ncld,Nhour,lat,lon)
                out = var_data[:, :, lat_ind, lon_ind]

                new_shape = (var_data.shape[0], var_data.shape[1]) + original_shape
                new_dims = [var.dims[0], var.dims[1]] + \
                        ["dim_" + str(i) for i in range(len(original_shape))]

            else:
                raise ValueError("Unsupported variable dimension.")

            # Reshape back to target shape
            out = out.reshape(new_shape)

            # construct the new coordinates
            coords = {}

            # CERES dims
            for d in var.dims:
                if d not in ("latitude", "longitude"):
                    coords[d] = var.coords[d]

            # New dims (corresponding to target shape)
            for i, nd in enumerate(new_dims[-len(original_shape):]):
                coords[nd] = np.arange(original_shape[i])

            # Target coordinate fields
            coords["target_lon"] = (new_dims[-len(original_shape):], target_lons)
            coords["target_lat"] = (new_dims[-len(original_shape):], target_lats)

            return xr.DataArray(out, dims=new_dims, coords=coords)

    def _geolocate_xarray(self, varname, target_lons, target_lats):
        """
        Xarray version of nearest-neighbour geolocation.
        Works for arbitrary-shaped input lon/lat arrays.
        """

        ds = xr.open_dataset(self.get_fileloc(), engine="netcdf4")
        da = ds[varname]

        # Save original shape
        original_shape = target_lons.shape
        assert target_lats.shape == original_shape

        # Flatten
        lon_flat = target_lons.ravel()
        lat_flat = target_lats.ravel()

        # Fix longitude convention
        if da.longitude.max() > 180:
            lon_flat = np.mod(lon_flat, 360)
        else:
            lon_flat = ((lon_flat + 180) % 360) - 180

        # Run vectorised selection
        out = da.sel(
            longitude=xr.DataArray(lon_flat, dims="points"),
            latitude=xr.DataArray(lat_flat, dims="points"),
            method="nearest",
        )

        # Now reshape output dims
        # Example out dims: (time, points) or (points)
        extra_dims = ["dim_" + str(i) for i in range(len(original_shape))]

        new_dims = list(out.dims[:-1]) + extra_dims
        new_shape = out.shape[:-1] + original_shape

        out = out.values.reshape(new_shape)

        # Build coords
        coords = {}

        # Existing CERES dims
        for d in da.dims:
            if d not in ("latitude", "longitude"):
                coords[d] = da.coords[d]

        # Extra dims
        for i, nd in enumerate(extra_dims):
            coords[nd] = np.arange(original_shape[i])

        coords["target_lon"] = (extra_dims, target_lons)
        coords["target_lat"] = (extra_dims, target_lats)

        return xr.DataArray(out, dims=new_dims, coords=coords)

    



        