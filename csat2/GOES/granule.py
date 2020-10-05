'''Provides simple functions for using GOES data'''
from .. import locator
import os
import os.path
from google.cloud import storage
from fnmatch import fnmatch
from netCDF4 import Dataset
import numpy as np
import scipy
import datetime
import re
import misc
from .util import get_resolution as util_get_resolution
import logging
import xarray as xr


def readin_radiances_filename(filename):
    with xr.open_dataset(filename) as ds:
        return ds['Rad'][:]


def available_files(product, year, doy, hour, minute, sat, channel, area='RadF', mode=6):
    if minute == '*':
        minute = '**'
    if hour == '*':
        hour = '**'
    filenames = locator.search(
        'GOES', product,
        year=year, doy=doy, hour=hour, minute=minute,
        channel=channel, area=area, mode=mode, sat=sat)
    return filenames


class Granule():
    '''A class for handling GOES images. Will increment through images
    and should be able to cope with most common location transformation
    tasks (including parallax correction. Will also download files,
    assuming you place an appropriate API key in your csat2 config
    directory.

    Note, the different scan areas are handled separately, it is
    assumed you wont be converting between them.'''
    inc_minutes = {'RadM1': 1, 'RadM2': 1, 'RadC': 5, 'RadF': 10}

    def __init__(self, sat, area, year, doy, hour, minute, locator=None):
        self.sat = sat
        self.area = area
        self.year = year
        self.doy = doy
        self.hour = hour
        # Minute is stored as the start time for that block
        #  Use the inc_minutes value to shift to the start of the block
        self.minute = (minute-(minute % self.inc_minutes[area]))
        self.locator = locator

    @classmethod
    def fromtext(cls, gran_text):
        m = re.search(
            'G(?P<sat>..).(?P<year>....)(?P<doy>...).(?P<hour>..)(?P<minute>..).(?P<area>....)', gran_text)
        return cls(sat='G'+m.group('sat'),
                   area=m.group('area'),
                   year=int(m.group('year')),
                   doy=int(m.group('doy')),
                   hour=int(m.group('hour')),
                   minute=int(m.group('minute')))

    @classmethod
    def fromfilename(cls, filename):
        filename = os.path.basename(filename)
        m = re.search('OR_ABI-(.*?)-(?P<area>.*?)-M(?P<mode>.)C(?P<channel>..)' +
                      '_G(?P<sat>..)_s(?P<year>....)(?P<doy>...)(?P<hour>..)(?P<minute>..)', filename)
        return cls(sat='G'+m.group('sat'),
                   area=m.group('area'),
                   year=int(m.group('year')),
                   doy=int(m.group('doy')),
                   hour=int(m.group('hour')),
                   minute=int(m.group('minute')))

    def astext(self):
        return '{}.{}{:0>3}.{:0>2}{:0>2}.{}'.format(
            self.sat, self.year, self.doy,
            self.hour, self.minute, self.area)

    def __repr__(self):
        return self.astext()

    def datetime(self):
        return (datetime.datetime(
            year=self.year, month=1, day=1,
            hour=self.hour, minute=self.minute) +
                datetime.timedelta(days=self.doy-1))

    def get_lonlat(self, channel, product='L1b', mode='*'):
        '''Proes the lon-lat arrays for a specific channel.
        Note that as the channels have different spatial resolutions,
        ahannel must be specified here.'''
        ipc = self.get_index(channel, product, mode)
        llp = self.geolocate(np.array(ipc).reshape(
            2, -1).transpose(), channel, product, mode)
        return llp[:, 0].reshape(ipc[0].shape), llp[:, 1].reshape(ipc[0].shape)

    def get_resolution(self, channel, product='L1b', mode='*', fromfile=False):
        if fromfile:
            with Dataset(self.get_filename(channel, product=product, mode=mode)) as ncdf:
                return ncdf.variables['Rad'].resolution.split(' ')[1]
        else:
            return util_get_resolution(channel)*0.000028

    def geolocate(self, coords, channel=None, product='L1b', mode='*', interp=False):
        '''Returns the lon/lat of the gridbox indicies passed as coords'''
        if not self.locator:
            self.locator = GOESLocator(self.get_filename(
                channel=channel, product=product, mode=mode))
        return self.locator.geolocate(coords, interp)

    def locate(self, coords, alt=None, channel=None, product='L1b', mode='*', interp=False):
        
        if not self.locator:
            self.locator = GOESLocator(self.get_filename(
                channel=channel, product=product, mode=mode))
        return self.locator.locate(coords, alt, interp)

    def points_in_radius(self):
        pass

    def next(self, number=1, only_downloaded=False, only_exisiting=False):
        '''Increment image name'''
        dt = self.datetime()
        dt += datetime.timedelta(minutes=number*self.inc_minutes[self.area])
        year, doy = misc.time.date_to_doy(dt.year, dt.month, dt.day)
        return Granule(self.sat, self.area,
                       year, doy, dt.hour, dt.minute,
                       # Reasonable assumption that coordinates remain the same
                       locator=self.locator)

    def get_filename(self, channel, product='L1b', mode='*'):
        '''This is more complicated that the simple locator, as we
        have to account for the GOES image not being taken at an exact 
        time. Even if we can trust the timing, changes in the scan pattern
        can cause a change in the image time.

        Selects an image with a mane that puts it within the increment 
        timestep.'''
        filenames = locator.search('GOES', product,
                                   year=self.year, doy=self.doy,
                                   hour=self.hour, minute=str(
                                       int(self.minute)//10)+'*',
                                   area=self.area, sat=self.sat,
                                   channel=channel, mode=mode)
        minutes = np.array([int(f.split('_')[-3][-5:-3]) for f in filenames])
        ind = np.where(
            np.logical_and(
                ((minutes-self.minute) >= 0),
                ((minutes-self.minute) < self.inc_minutes[self.area])))
        if len(ind[0]) == 0:
            raise IndexError('No matching file')
        if len(ind[0]) > 1:
            raise IndexError('Non-unique filename')
        return filenames[ind[0][0]]

    def get_band_radiance(self, channel, product='L1b', mode='*', refl=False):
        with xr.open_dataset(
                self.get_filename(channel, product=product, mode=mode)) as ds:
            data = ds['Rad'][:]
            if refl:
                data *= ds['kappa0']
            return data

    def get_band_bt(self, channel, product='L1b', mode='*'):
        '''Conversion of radiances to BT using constants in files'''
        with xr.open_dataset(
                self.get_filename(channel, product=product, mode=mode)) as ds:
            data = ds['Rad'][:]
            fk1, fk2 = (ds['planck_fk1'],
                        ds['planck_fk2'])
            bc1, bc2 = (ds['planck_bc1'],
                        ds['planck_bc2'])
            logging.info(fk1, fk2, bc1, bc2)
            bt = fk2/(np.log((fk1/(data))+1))
            bt = (bt-bc1)/bc2
            return bt

    def get_shape(self, channel, product='L1b', mode='*'):
        with xr.open_dataset(
                self.get_filename(channel, product=product, mode=mode)) as ds:
            return ds['Rad'].shape

    def get_index(self, channel, product='L1b', mode='*'):
        granshape = self.get_shape(channel, product, mode)
        return np.array(np.meshgrid(np.arange(0, granshape[1]),
                                    np.arange(0, granshape[0])))

    def get_viewangles(self, channel, product='L1b', mode='*'):
        with xr.open_dataset(
                self.get_filename(channel, product=product, mode=mode)) as ds:
            return np.array(np.meshgrid(ds['x'][:],
                                        ds['y'][:]))

    def get_llcoord(self, channel, product='L1b', mode='*'):        
        with xr.open_dataset(
                self.get_filename(channel, product=product, mode=mode)) as ds:
            return ds['x'][0], ds['y'][-1]

    def download(self, channel, product='L1b', mode='*'):
        local_folder = locator.get_folder('GOES', product,
                                          year=self.year, doy=self.doy,
                                          hour=self.hour, minute=self.minute,
                                          channel=channel, area=self.area,
                                          mode=mode, sat=self.sat)

        # Create the local folder (if required)
        try:
            os.makedirs(local_folder)
        except FileExistsError:
            pass

        dl_minute = str(int(self.minute)//10)+'*'

        bucket = 'gcp-public-data-goes-{sat}'.format(sat=self.sat[1:])
        pattern = ('ABI-{product}-{area:.4}/{year}/{doy:0>3}/' +
                   '{hour:0>2}/OR_ABI-{product}-{area}-M{mode}C{channel:0>2}_' +
                   'G{sat}_s{year}{doy:0>3}{hour:0>2}{minute:0>2}*').format(
                       product=product,
                       year=self.year, doy=self.doy,
                       hour=self.hour, minute=dl_minute,
                       channel=channel, area=self.area,
                       mode=mode, sat=self.sat[1:])
        prefix = os.path.dirname(pattern)
        logging.debug(prefix)

        # You will need to create a credential file for a blank project
        # There are no permissions here, but otherwise this fails as you
        # have no authentication file. It may be possible to do this with
        # a blank file?
        storage_client = storage.Client.from_service_account_json(
            os.path.expandvars('${HOME}/.csat2/goes-test.json'))
        logging.debug(storage_client)
        files = storage_client.get_bucket(bucket).list_blobs(prefix=prefix)
        
        for fl in files:
            logging.debug(fl.name)
            if fnmatch(fl.name, pattern):
                re_minutes = int(fl.name.split('_')[-3][-5:-3])
                if (((re_minutes-self.minute) >= 0) and ((re_minutes-self.minute) < self.inc_minutes[self.area])):
                    logging.debug(os.path.basename(fl.name))
                    fl.download_to_filename(
                        os.path.join(local_folder, os.path.basename(fl.name)))


class GOESLocator():
    '''For converting between geodetic and satellite 
    coordinates, following
    https://www.goes-r.gov/users/docs/PUG-L1b-vol3.pdf'''

    def __init__(self, filename):
        m = re.search('OR_ABI-(.*?)-(?P<area>.*?)-M(?P<mode>.)C(?P<channel>..)' +
                      '_G(?P<sat>..)_s(?P<year>....)(?P<doy>...)(?P<hour>..)(?P<minute>..)', filename)
        self.res = util_get_resolution(int(m.group('channel')))
        with Dataset(filename) as ncdf:
            ncdf.set_auto_mask(False)
            gvar = ncdf.variables['goes_imager_projection']
            self.req = gvar.semi_major_axis
            self.invf = gvar.inverse_flattening
            self.rpol = gvar.semi_minor_axis
            self.e = 0.0818191910435
            self.pph = gvar.perspective_point_height
            self.H = self.pph + self.req
            self.lam0 = np.deg2rad(gvar.longitude_of_projection_origin)
            self.x = ncdf.variables['x'][:].astype('float64')
            self.y = ncdf.variables['y'][:].astype('float64')
            self.xint = None
            self.corrdx = None

    def _create_interp(self):
        self.xint = scipy.interpolate.interp1d(
            range(len(self.x)), self.x, fill_value='extrapolate')
        self.yint = scipy.interpolate.interp1d(
            range(len(self.y)), self.y, fill_value='extrapolate')
        self.rxint = scipy.interpolate.interp1d(
            self.x, range(len(self.x)), fill_value='extrapolate')
        self.ryint = scipy.interpolate.interp1d(
            self.y, range(len(self.y)), fill_value='extrapolate')

    def _sat_to_geodetic(self, x, y):
        '''x, y - the longitude and latitude in the satellite
        frame (in radians)'''
        a = np.sin(x)**2 + np.cos(x)**2 * (np.cos(y)**2 +
                                           (self.req**2/self.rpol**2)*np.sin(y)**2)
        b = -2*self.H*np.cos(x)*np.cos(y)
        c = self.H**2 - self.req**2
        rs = (-b - np.sqrt(b**2 - 4 * a * c))/(2*a)
        sx = rs*np.cos(x)*np.cos(y)
        sy = -rs*np.sin(x)
        sz = rs*np.cos(x)*np.sin(y)

        lon = self.lam0-np.arctan(sy/(self.H-sx))
        lat = np.arctan((self.req**2*sz)/(self.rpol**2 *
                                          np.sqrt((self.H-sx)**2+sy**2)))
        return lon, lat

    def _geodetic_to_sat(self, lon, lat):
        phic = np.arctan((self.rpol/self.req)**2*np.tan(lat))
        rc = self.rpol/np.sqrt(1-(self.e*np.cos(phic))**2)
        sx = self.H - rc*np.cos(phic)*np.cos(lon-self.lam0)
        sy = -rc*np.cos(phic)*np.sin(lon-self.lam0)
        sz = rc*np.sin(phic)

        x = np.arcsin(-sy/np.sqrt(sx**2+sy**2+sz**2))
        y = np.arctan(sz/sx)
        # Mask invisible points
        mask = np.where((self.H*(self.H-sx)) <
                        (sy**2+(self.req*sz/self.rpol)**2))
        x[mask] = np.nan
        y[mask] = np.nan
        return x, y

    def _rebin(self, locs, bins):
        binned = np.digitize(locs, bins)-1
        binned[locs < min(bins)] = -999
        binned[locs > max(bins)] = -999
        return binned

    def _create_alt_correction(self):
        # Angles as in logbook
        y, x = np.array(np.meshgrid(self.y, self.x, indexing='ij'))
        ipc = np.array(np.meshgrid(
            np.arange(0, len(self.x)),
            np.arange(0, len(self.y))))
        llp = self.geolocate(
            np.array(ipc).reshape(2, -1).transpose(),
            interp=False)
        lon = llp[:, 0].reshape(ipc[0].shape)
        theta = np.deg2rad(lon) - self.lam0
        phi = np.deg2rad(llp[:, 1].reshape(ipc[0].shape))

        # Correction for the satellite parallax
        # (km x-y per km height)
        self.corrdx = (np.sin(x)*np.cos(theta) +
                       np.cos(x)*np.sin(theta))*np.cos(phi)
        self.corrdy = (np.cos(x)*np.sin(y)*np.cos(theta)*np.cos(phi) -
                       np.sin(x)*np.sin(y)*np.sin(theta)*np.cos(phi) -
                       np.cos(y)*np.sin(phi))

    def geolocate(self, coords, interp=False):
        '''Returns the lon/lat of the gridbox indicies passed as coords'''
        if interp:
            if not self.xint:
                self._create_interp()
            x = self.xint(
                np.clip(coords[:, 0].astype('int'), 0, len(self.x)-1))
            y = self.yint(
                np.clip(coords[:, 1].astype('int'), 0, len(self.y)-1))
        else:
            x = self.x[np.clip(coords[:, 0].astype('int'), 0, len(self.x)-1)]
            y = self.y[np.clip(coords[:, 1].astype('int'), 0, len(self.y)-1)]
        lon, lat = self._sat_to_geodetic(x, y)
        return np.concatenate((np.rad2deg(lon)[:, None],
                               np.rad2deg(lat)[:, None]), axis=1)

    def locate(self, coords, alt=None, interp=False):
        '''Get the satellite coordinates for N lon-lat coords, shape - [N, 2]
        Performs a simple altitude parallax adjustment, passing altitude as km'''
        x, y = self._geodetic_to_sat(np.deg2rad(coords[:, 0]),
                                     np.deg2rad(coords[:, 1]))
        if interp:
            if not self.xint:
                self._create_interp()
            output = np.concatenate((self.rxint(x)[:, None],
                                     self.ryint(y)[:, None]), axis=1)
        else:
            output = np.concatenate((self._rebin(x, self.x)[:, None],
                                     self._rebin(y, self.y)[:, None]), axis=1)
        if alt is not None:
            # Correct for altitude parallax
            # Following method laid out in notebook
            # Vector dot product
            if self.corrdx is None:
                self._create_alt_correction()
            # store locations outside of grid
            ctlpos = output.astype('float')
            ctlpos[(ctlpos < 0).sum(axis=-1) > 0] = np.nan
            ctlmask = np.isnan(ctlpos)

            ctl_ind = output.astype('int')
            ctl_ind[ctlmask] = 0

            ctlpos_corrdx = (alt/self.res)*self.corrdx[
                ctl_ind[:, 1].astype('int'),
                ctl_ind[:, 0].astype('int')]
            ctlpos_corrdx[ctlmask[:, 0]] = np.nan
            ctlpos_corrdy = (alt/self.res)*self.corrdy[
                ctl_ind[:, 1].astype('int'),
                ctl_ind[:, 0].astype('int')]
            ctlpos_corrdy[ctlmask[:, 1]] = np.nan
            ctlpos_corr = np.concatenate(
                (ctlpos_corrdx[:, None],
                 ctlpos_corrdy[:, None]), axis=-1)
            output = ctlpos + ctlpos_corr
        return output
