'''ECMWF.py
Utility files for using ECMWF data

Edward Gryspeerdt, Imperial College London, 2020
'''
from .. import locator
from csat2.misc.time import ydh_to_datetime
import numpy as np
from netCDF4 import Dataset, num2date
from csat2 import misc
import datetime
import xarray as xr
import logging


def readin(product, *args, **kwargs):
    if product in ['ERAInterim', 'ERA40', 'ERA5']:
        return readin_ERA(product=product, *args, **kwargs)
    # if product=='MACC-anth':
    #     return readin_MACCanth(*args, **kwargs)
    # if product=='MACC':
    #     return readin_MACC_aod(*args, **kwargs)
    else:
        print('Not a suitable product')


########################
# ERA readin functions #
########################
def readin_ERA(year, doy, variable, level=None, product='ERAInterim', time='LST',
               resolution='1grid', varname=None, *args, **kwargs):
    if variable == 'LTS':
        return readin_ERA_LTS(
            year, doy, product=product, time=time,
            resolution=resolution, *args, **kwargs)
    if variable == 'EIS':
        return readin_ERA_EIS(
            year, doy, product=product, time=time,
            resolution=resolution, *args, **kwargs)

    if not level:
        raise ValueError('level is required if variable not LTS or EIS')
    filename = locator.search('ECMWF', product,
                              year=year,
                              doy=doy,
                              variable=variable,
                              level=level,
                              resolution=resolution,
                              time=time)
    if (resolution is '1grid'):  # and (returndata['time'].dtype == np.float):
        # These should match the MODIS files, so will build the xarray directly
        # xarray is very slow at transposing arrays (for some reason)
        with Dataset(filename[0]) as ncdf:
            returndata = ncdf.variables[variable_names(ncdf)][:]
            try:
                rtime = num2date(
                    ncdf.variables['time'][:],
                    units=ncdf.variables['time'].units)
            except AttributeError:  # Times in files is in hour strings
                rtime = np.array([
                    ydh_to_datetime(year, doy,
                                    (t//100)+(t % 100)/60)
                    for t in ncdf.variables['time'][:]])
            rlon = ncdf.variables['lon'][:]
            if time == 'LST':
                return xr.DataArray(np.swapaxes(returndata, 1, 2),
                                    name=variable_names(ncdf),
                                    coords=[('time', rtime),
                                            ('lon', rlon),
                                            ('lat', ncdf.variables['lat'][:])])
            cycle = list(range(180, 360))+list(range(180))
            rlon = rlon[cycle]
            rlon = np.where(rlon > 180, rlon-360, rlon)
            return xr.DataArray(np.swapaxes(returndata, 1, 2)[:, cycle][:, :, ::-1],
                                name=variable_names(ncdf),
                                coords=[('time', rtime),
                                        ('lon', rlon),
                                        ('lat', ncdf.variables['lat'][::-1])])
    else:
        ds = xr.open_dataset(filename[0])
        varname = list(ds.data_vars.keys())[0]
        returndata = ds[varname][:]
        ds.close()
    return returndata


def variable_names(ncdf):
    ''' Returns the valid dataset variables in the file.'''
    pvars = [f for f in ncdf.variables.keys()
             if f not in ['lat', 'lon', 'time',
                          'latitude', 'longitude',
                          'level', 'plev',
                          'lat_bnds', 'lon_bnds']]
    if len(pvars) == 1:
        varname = pvars[0]
        return varname
    else:
        raise(KeyError,
              'Unable to infer variable name, {}'.format(pvars))


def readin_ERA_LTS(year, doy, product='ERAInterim', time='LST', resolution='1grid'):
    return _calc_lts(
        readin_ERA(year, doy, 'Temperature', '700hPa',
                   product=product, time=time,
                   resolution=resolution),
        readin_ERA(year, doy, 'Temperature', '1000hPa',
                   product=product, time=time,
                   resolution=resolution)).rename('lts')


def _es(temp):
    '''Returns SVP in Pa if temp is in K
    from Alduchov and Eskridge (1996)'''
    A, B, C = 6.1094, 17.625, 243.04
    t = temp-273.15
    return A*np.exp((B*t)/(C+t))*100


def _td(temp, rh):
    A, B, C = 6.1094, 17.625, 243.04
    temp = temp-273.15
    return C*(np.log(rh/100)+(B*temp)/(C+temp))/(B-np.log(rh/100)-B*temp/(C+temp))+273.15


def _sat_adiabatic_potential_gradient(T, p, rh):
    g = 9.81
    Lv = (2.501 - 0.00237*(T-273.15)) * 10**6  # From Bolton

    ep = 0.622
    es = 6.112*np.exp((17.67*(T-273.15))/((T-273.15)+243.5))*100  # From Bolton
    qs = ep * es/(p-es)

    Rsd = 287
    Rsw = 461.5
    cpd = 1005.7  # Bolton
    A = (Lv*qs/(Rsd*T))
    B = ((Lv*Lv*qs)/(cpd*Rsw*T*T))
    return (g/cpd) * (1 - (1+A)/(1+B))


def _calc_lcl(t1000, rh1000):
    # FAA calculation, assumes dry adiabat below LCL
    tlcl = 1/(1/(t1000-55) - np.log(rh1000/100)/2840) + 55
    return (t1000 - tlcl)/9.8 * 1000


def _calc_lts(t1000, t700):
    return 1.1082*t700 - t1000


def _calc_eis(t700, t1000, rh1000=80):
    # An averge value gives a slightly better match to WH06
    z700 = (287*0.5*(t700+t1000)/9.81)*np.log(1000/700)
    # z700 = (287*0.5*(t700+t1000)/9.81)*np.log(1010/700) # To match WH06
    g850 = _sat_adiabatic_potential_gradient(0.5*(t700+t1000), 85000, rh1000)
    # td = _td(t1000, rh1000)
    lcl = _calc_lcl(t1000, rh1000)
    lts = _calc_lts(t1000, t700)
    # lts = 1.11075*t700 - t1000 # To match WH06
    eis = lts - g850*z700+g850*lcl
    return lts, eis, g850*1000


def readin_ERA_EIS(year, doy, product='ERAInterim', time='LST', resolution='1grid', lcl='FAA', with_rh=True):
    t700 = readin_ERA(year, doy, 'Temperature', '700hPa',
                      product=product, time=time,
                      resolution=resolution)
    t1000 = readin_ERA(year, doy, 'Temperature', '1000hPa',
                       product=product, time=time,
                       resolution=resolution)
    if with_rh:
        rh1000 = readin_ERA(year, doy, 'Relative_humidity', '1000hPa',
                            product=product, time=time,
                            resolution=resolution)
        return _calc_eis(t700, t1000, rh1000)[1].rename('eis')
    else:
        return _calc_eis(t700, t1000)[1].rename('eis')


class ERA5Data():
    '''A class for using ERA5 Data. Holds the ncdf files open and updates them
    as necessary. Does not do any explict scaling, but netCDF4 will account
    for that.'''

    def __init__(self, variable, level, res='0.25grid', scaling=1, linear_interp=False):
        self.nc = None
        self.time = None
        self.lon = None
        self.lat = None
        self.year = None
        self.doy = None
        self.scaling = scaling
        self.variable = variable
        self.varkey = None
        self.level = level
        self.res = res
        # Define linear interp here, as if it is True or 'time'
        # we need to cache some extra data
        self.linear_interp = linear_interp

    def get_data_time(self, time):
        '''Returns the wind fields at a given time'''
        year, doy = time.year, time.timetuple().tm_yday
        if (doy != self.doy) or (year != self.year):
            self.update_files(year, doy)

        if self.linear_interp in ['time', 'both']:
            time_ind = np.argmin(time >= self.time) - 1
            if time_ind <= len(self.time)-1:
                # The weight for data at the first index
                weight = 1-((time-self.time[time_ind]) /
                            (self.time[time_ind+1]-self.time[time_ind]))
                u_ret = (
                    (weight*self.nc.variables[
                        self.varkey][time_ind, :, :]) +
                    ((1-weight)*self.nc.variables[
                        self.varkey][time_ind+1, :, :]))
            else:
                u_ret = self.nc.variables[self.varkey][-1, :, :]

        else:
            # No temporal interpolation, take the nearest index
            time_ind = np.argmin(np.abs(self.time - time))
            u_ret = self.nc.variables[self.varkey][time_ind, :, :]

        if (self.res == '1grid'):  # Arrange the data to match the MODIS files
            rlon = self.lon
            cycle = list(range(180, 360))+list(range(180))
            rlon = rlon[cycle]
            rlon = np.where(rlon > 180, rlon-360, rlon)
            return xr.DataArray(np.swapaxes(u_ret, 0, 1)[cycle][:, ::-1],
                                coords=[('lon', rlon),
                                        ('lat', self.lat[::-1])])
        else:
            return xr.DataArray(
                u_ret,
                dims=['lat', 'lon'],
                coords={'lon': self.lon, 'lat': self.lat})

    def get_data(self, lon, lat, time, simple=False):
        year, doy = time.year, time.timetuple().tm_yday
        if (doy != self.doy) or (year != self.year):
            self.update_files(year, doy)

        # Take the nearest index - note that the self.lon and self.lat
        # arrays store the gridbox centres, not the edges

        s_interp = False
        if self.linear_interp in ['space', 'both']:
            lat_ind = np.digitize(lat, self.lat)-1
            lat_weight = (np.mod(lat-self.lat[lat_ind], 1)/self.lat_inc)[..., None]

            lat_ind = np.repeat(lat_ind[..., None], 4, axis=-1)
            lat_ind[..., 2:] += 1
            lat_ind = np.clip(lat_ind, 0, len(self.lat)-1)

            lon_ind = np.digitize(np.mod(lon, 360), self.lon)-1
            lon_weight = (np.mod(lon - self.lon[lon_ind], 1)/self.lon_inc)[..., None]

            lon_ind = np.repeat(lon_ind[..., None], 4, axis=-1)
            lon_ind[..., [1, 3]] += 1
            lon_ind = np.mod(lon_ind, len(self.lon))
            
            weights = np.ones(lon_ind.shape)
            weights[..., [0, 2]] *= (1-lon_weight)
            weights[..., [1, 3]] *= lon_weight
            weights[..., :2] *= (1-lat_weight)
            weights[..., 2:] *= lat_weight
            
            s_interp = True
        else:
            lat_ind = np.digitize(lat, self.lat-0.5*self.lat_inc)-1
            lon_ind = np.digitize(np.mod(lon, 360), self.lon -
                                  0.5*self.lon_inc)-1

        t_interp = False
        if self.linear_interp in ['time', 'both']:
            time_ind = np.argmin(time >= self.time)-1
            if time_ind <= len(self.time)-1:
                # The weight for data at the first index
                tweight = 1-((time-self.time[time_ind]) /
                             (self.time[time_ind+1]-self.time[time_ind]))
                t_interp = True
        else:
            time_ind = np.argmin(np.abs(self.time - time))
            
        if simple:
            if len(self.nc.variables[self.varkey].shape) > 3:
                u_ret = self.scaling*self.nc.variables[self.varkey][
                    time_ind, 0, :, :][lat_ind, lon_ind]
            else:
                u_ret = self.scaling*self.nc.variables[self.varkey][
                    time_ind, :, :][lat_ind, lon_ind]
            return u_ret

        else:
            # Calculate the slice for the data to read in
            # - Done as the netCDF4 library is too smart
            # - It wont allow lists of indicies
            lat_range = lat_ind.min(), lat_ind.max()+1
            lon_range = lon_ind.min(), lon_ind.max()+1

            if len(self.nc.variables[self.varkey].shape) > 3:
                u_ret = self.scaling*self.nc.variables[self.varkey][
                    time_ind, 0, lat_range[0]:lat_range[1], lon_range[0]:lon_range[1]][
                        lat_ind-lat_range[0],
                        lon_ind-lon_range[0]]
                if t_interp:
                    u_retp = self.scaling*self.nc.variables[self.varkey][
                        time_ind+1, 0, lat_range[0]:lat_range[1], lon_range[0]:lon_range[1]][
                            lat_ind-lat_range[0],
                            lon_ind-lon_range[0]]
                    u_ret = tweight*u_ret + (1-tweight)*u_retp
            else:
                u_ret = self.scaling*self.nc.variables[self.varkey][
                    time_ind, lat_range[0]:lat_range[1], lon_range[0]:lon_range[1]][
                        lat_ind-lat_range[0],
                        lon_ind-lon_range[0]]
                if t_interp:
                    u_retp = self.scaling*self.nc.variables[self.varkey][
                        time_ind+1, lat_range[0]:lat_range[1], lon_range[0]:lon_range[1]][
                            lat_ind-lat_range[0],
                            lon_ind-lon_range[0]]
                    u_ret = tweight*u_ret + (1-tweight)*u_retp

            if s_interp:
                # This is where the averaging is done
                return np.sum(u_ret*weights, axis=-1)
            return u_ret

    def update_files(self, year, doy):
        # Need to load new data!
        nyear, ndoy = misc.time.doy_step(year, doy, 1)
        # u wind nc
        if self.nc:
            self.nc.close()
        self.nc = Dataset(locator.search(
            'ECMWF', 'ERA5',
            year=year, doy=doy,
            variable=self.variable,
            level=self.level,
            resolution=self.res,
            time='timed')[0])
        self.nc.set_auto_mask(False)
        self.varkey = variable_names(self.nc)
        # Time and date
        self.year = year
        self.doy = doy
        try:
            self.lon = self.nc.variables['lon'][:]
            self.lat = self.nc.variables['lat'][:]
        except KeyError:
            self.lon = self.nc.variables['longitude'][:]
            self.lat = self.nc.variables['latitude'][:]
        self.time = num2date(self.nc.variables['time'][:],
                             units=self.nc.variables['time'].units)
        if (len(set(np.diff(self.lon))) != 1) or (len(set(np.diff(self.lat))) != 1):
            raise ValueError('Must be an evenly girdded file')
        self.lon_inc = self.lon[1]-self.lon[0]
        self.lat_inc = self.lat[1]-self.lat[0]
        

class ERA5WindData():
    def __init__(self, level='850hPa', res='0.25grid', wind_scaling=(1, 1), *args, **kwargs):
        self.udata = ERA5Data('U-wind-component', level=level, res=res,
                              scaling=wind_scaling[0],
                              *args, **kwargs)
        self.vdata = ERA5Data('V-wind-component', level=level, res=res,
                              scaling=wind_scaling[1],
                              *args, **kwargs)

    def get_data_time(self, time):
        '''Returns the wind fields at a given time'''
        return (self.udata.get_data_time(time),
                self.vdata.get_data_time(time))

    def get_data(self, lon, lat, time, simple=False):
        return (
            self.udata.get_data(
                lon, lat, time, simple=simple),
            self.vdata.get_data(
                lon, lat, time, simple=simple))


##################
# MACC functions #
##################

def readin_MACCanth(year, doy):
    filename = locator.search('ECMWF', 'MACC-anth',
                              year=year,
                              doy=doy)
    try:
        fname = filename[0]
        with Dataset(fname, 'r') as ncdf:
            outdata = {}
            cycle = list(range(180, 360))+list(range(180))
            outdata['aod'] = ncdf.variables['aod'][:].squeeze().transpose()[
                cycle][:, ::-1]
            outdata['aaod'] = ncdf.variables['anthaod'][:].squeeze().transpose()[
                cycle][:, ::-1]
            return outdata
    except IndexError:
        return None


def readin_MACC_aod(year, doy, sds, time='LST'):
    if time == 'LST':
        time = '_'+time
    filename = locator.search('ECMWF', 'MACC-surf',
                              time=time, year=year, doy=doy)
    try:
        fname = filename[0]
        with Dataset(fname, 'r') as ncdf:
            outdata = {}
            cycle = list(range(180, 360))+list(range(180))
            for name in sds:
                outdata[name] = ncdf.variables[name][:].squeeze().transpose()[
                    cycle][:, ::-1]
            return outdata
    except IndexError:
        return None
