from datetime import datetime, timedelta
import misc.time
from .. import locator
from .readfiles import (readin, readin_MODIS_L2, readin_metadata, field_interpolate, DEFAULT_COLLECTION)
from .util import band_centres, bowtie_correct
from .download import download, check
import numpy as np
from sklearn.neighbors import BallTree
import scipy.constants
from netCDF4 import Dataset


class Granule(object):
    """An object for plotting data for a MODIS granule"""
    def __init__(self, year, doy, time, sat, col=DEFAULT_COLLECTION):
        '''year, doy, time - all int
        sat = 'A' or 'T' '''
        super(Granule, self).__init__()
        self.year = year
        self.doy = doy
        self.time = time
        if sat not in ['A', 'T']:
            raise ValueError('Sat must be A or T') 
        self.sat = sat
        self.locator = None
        self.rect_locator = None
        self.lonlat = None
        self._orbit = None
        self._daynight = None
        # I don't really want the collection to be here, as it is not independent
        # of the data itself. However, I think it is the most sensible place to
        # have it. I can't think of many situations were you would use data from
        # different collections together (except perhaps geolocation data)
        self.col = col

    @classmethod
    def fromtext(cls, gran_text, col=DEFAULT_COLLECTION):
        return cls(int(gran_text[:4]),
                   int(gran_text[4:7]),
                   int(gran_text[8:12]),
                   gran_text[-1], col)

    @classmethod
    def fromfilename(cls, filename):
        basename = filename.split('/')[-1]
        sat = {'O': 'T', 'Y': 'A'}[basename[1]]
        year = int(basename.split('.')[1][1:5])
        doy = int(basename.split('.')[1][5:8])
        time = int(basename.split('.')[2])
        col = str(int(basename.split('.')[3]))
        return cls(year, doy, time, sat, col)

    def orbit(self):
        '''Return the satellite orbit number'''
        if not self._orbit:
            dat = readin('GEOMETA', self.year, self.doy,
                         {'A': 'aqua',
                          'T': 'terra'}[self.sat])
            self._orbit = dat['Orbit'][
                np.where(
                    dat['StartTime'] == self.datetime().strftime(
                        '%Y-%m-%d %H:%M'))[0]][0]
        return self._orbit

    def daynight(self):
        '''Is this granule day or night?'''
        if not self._daynight:
            dat = readin('GEOMETA', self.year, self.doy,
                         {'A': 'aqua',
                          'T': 'terra'}[self.sat])
            self._daynight = dat['DayNight'][
                np.where(
                    dat['StartTime'] == self.datetime().strftime(
                        '%Y-%m-%d %H:%M'))[0]][0]
        return self._daynight
    
    def timestr(self):
        '''Formatted time string'''
        return '{:0>4}'.format(self.time)

    def astext(self):
        return '{}{:0>3}.{:0>4}{}'.format(self.year, self.doy, self.time, self.sat)

    def datetime(self):
        '''The datetime of the Granule start'''
        _, mon, day = misc.time.doy_to_date(self.year, self.doy)
        return datetime(self.year, mon, day, int(self.time)//100, int(self.time) % 100)

    # def locate(self, locs, rectified=False, **kwargs):
    #     '''Returns the locations in the granule of lon,lat pairs,
    #     is an instance of MODISlocator'''
    #     if not self.locator:
    #         self.locator = _MODISlocator(self, rectified=rectified, col=self.col)
    #     return self.locator.locate(locs, **kwargs)

    def locate(self, locs, rectified=False, product=None, col=None, factor=2, **kwargs):
        '''Returns the locations in the granule of lon,lat pairs,
        is an instance of MODISlocator'''
        if rectified:
            #print('Rectified')
            if not self.rect_locator:
                #print('Rectified create')
                self.rect_locator = _MODISlocator(self, rectified=rectified, product=product, col=col, factor=factor)
            return self.rect_locator.locate(locs, **kwargs)
        else:
            if not self.locator:
                self.locator = _MODISlocator(self, product=product, col=col, factor=factor)
            return self.locator.locate(locs, **kwargs)

    def points_in_radius(self, loc, dist=20):
        if not self.locator:
            self.locator = _MODISlocator(self, col=self.col)
        return self.locator.points_in_radius(loc, dist)

    def _read_lonlat(self, product=None, col=None, dateline=False):
        if not col:
            col = self.col
        if not product:
            product = '021KM'
        data = readin_MODIS_L2(
            self.product_expand(product),
            self.year,
            self.doy,
            sds=[],
            times=[self.timestr()], col=col)
        if dateline:
            # Ignores the poles
            if (data['Longitude'].max()>160) and (data['Longitude'].min()<-160):
                data['Longitude'] = np.where(data['Longitude']<0,
                                             data['Longitude']+360,
                                             data['Longitude'])
        self.lonlat = [field_interpolate(data['Longitude']),
                       field_interpolate(data['Latitude'])]

    def get_lonlat(self, product=None, col=None, dateline=False):
        '''Get lon lat data - can specify the product or the collection
        used to the geolocation data from'''
        if not self.lonlat:
            self._read_lonlat(product=product, col=col, dateline=dateline)
        return self.lonlat

    def mloc_to_lonlat(self, mlocs):
        '''The the lon/lat positions for a given set of MODIS coordinates'''
        if not self.lonlat:
            self._read_lonlat()
        mlocslon = np.array(mlocs)[:, 0].clip(0, self.lonlat[0].shape[0]-1)
        mlocslat = np.array(mlocs)[:, 1].clip(0, self.lonlat[0].shape[1]-1)
        return list(zip(self.lonlat[0][mlocslon, mlocslat],
                        self.lonlat[1][mlocslon, mlocslat]))

    def get_angles(self, mst=False, product=None, col=None, dateline=False):
        '''Get the solar and scattering angles for a MODIS/VIIRS image
        Returns solar zenith, azimuth, satellite zenith, azimuth

        mst - use the MODIS science team values from MOD03'''
        if not col:
            col = self.col
        if mst:
            print('MST data')
            angles = readin_MODIS_L2(
                self.product_expand('03'),
                self.year, self.doy,
                sds = ['SolarZenith', 'SolarAzimuth',
                       'SensorZenith', 'SensorAzimuth'],
                times=[self.timestr()],
                col=col)
            sunz = angles['SolarZenith'][:, :-4]/100
            suna = angles['SolarAzimuth'][:, :-4]/100
            satz = angles['SensorZenith'][:, :-4]/100
            sata = angles['SensorAzimuth'][:, :-4]/100
        else:
            import pyorbital.astronomy, pyorbital.orbital
            
            lon, lat = self.get_lonlat(product=product, col=col, dateline=dateline)
            # Shift to time at granule centre
            # This gives an error of around 0.1 degrees compared to MODIS science team zenith
            # TODO: should this be more accurate
            time = self.datetime()+timedelta(seconds=300)
            sunz = pyorbital.astronomy.sun_zenith_angle(time, lon, lat)[:, :-5]
            # The solar azimuth error is larger when assuming an incorrect time (a few degrees)
            # This may need to be changed
            suna = pyorbital.astronomy.get_alt_az(time, lon, lat)[1]
            suna = np.rad2deg(suna)[:, :-5]

            sat_alt = {'A': 703.,
                       'T': 709.,
                       'N': 833.,
                       'J': 826.}[self.sat] # Altitude in km
            sat_lon = lon[:, lon.shape[1]//2]
            sat_lat = lat[:, lat.shape[1]//2]

            sata, satel = pyorbital.orbital.get_observer_look(
                sat_lon[:, None],
                sat_lat[:, None],
                sat_alt,
                self.datetime()+timedelta(seconds=300),
                lon, lat, 0)
            # Sensor zenith is accurate to less than a degree
            # A 10km height error gives a 0.2 error in satz towards swath edges
            # Note - these are almost constant for each row in the granule
            satz = 90 - satel[:, :-5]
            # Satellite azimuth is definitely not constant by granule. It will change
            # as you get near the pole (for example)
            sata = sata[:, :-5]
            sata = np.where(np.isfinite(sata), sata, 0)
            satz = np.where(np.isfinite(satz), satz, 0)
            
        # Do satellite-solar azimuth difference correctly (to within 0-180)
        # This is solar - satellite.... - not it isn't, note the 180 in the last line
        azidiff = suna-sata
        azidiff = np.mod(azidiff, 360)
        # Note this is 180-, which gives satellite-solar
        azidiff = 180-np.where(azidiff>180, 360-azidiff, azidiff)

        return sunz, suna, satz, sata, azidiff


    def get_cdnc(self, band=21, best=False, pcl=False, bowtie_corr=True):
        '''Calculate the CDNC for a MODIS granule for a 
        band in [16, 21, 37].
        best - Use the sampling from Grosvenor e 92018)
        pcl - use partially cloud pixels
        bowtie_corr - correct for bowtie effect (default=True)

        Currently does not perform temperature correction'''
        band_suffix = {21: '', 37: '_37', 16: '_16'}[band]
        sds = ['Cloud_Effective_Radius'+band_suffix,
               'Cloud_Optical_Thickness'+band_suffix,
               'Cloud_Phase_Optical_Properties']
        if best:
            sds += ['Cloud_Fraction', 'Cloud_Mask_SPI',
                    'Solar_Zenith_Day', 'Sensor_Zenith_Day']
        if pcl:
            sds += ['Cloud_Effective_Radius'+band_suffix+'_PCL',
                    'Cloud_Optical_Thickness'+band_suffix+'_PCL',
                    'Cloud_Phase_Optical_Properties']
        data = readin_MODIS_L2(
            self.product_expand('06_L2'),
            self.year,
            self.doy,
            sds=sds,
            times=[self.timestr()], col=self.col)
        if pcl:
            data['Cloud_Optical_Thickness'+band_suffix] = np.where(
                np.isfinite(data['Cloud_Optical_Thickness'+band_suffix]),
                data['Cloud_Optical_Thickness'+band_suffix],
                data['Cloud_Optical_Thickness'+band_suffix+'_PCL'])
            data['Cloud_Effective_Radius'+band_suffix] = np.where(
                np.isfinite(data['Cloud_Effective_Radius'+band_suffix]),
                data['Cloud_Effective_Radius'+band_suffix],
                data['Cloud_Effective_Radius'+band_suffix+'_PCL'])        
        cod = np.where(data['Cloud_Phase_Optical_Properties'] == 2,
                       data['Cloud_Optical_Thickness'+band_suffix], np.nan)
        re = 1e-6 * data['Cloud_Effective_Radius'+band_suffix]
        if best:
            mask = ((cod[:, :-4] > 4) *
                    (re[:, :-4] > 4e-6) *
                    (field_interpolate(data['Cloud_Fraction']) > 0.9) *
                    (data['Cloud_Mask_SPI'][:, :-4, 0] < 30) *
                    (field_interpolate(data['Solar_Zenith_Day']) < 65))
            if bowtie_corr:
                return (bowtie_correct((1.37e-11 * cod**0.5 * re**-2.5)),
                        bowtie_correct(mask))
            else:
                return ((1.37e-11 * cod**0.5 * re**-2.5),
                        mask)
        if bowtie_corr:
            return bowtie_correct((1.37e-11 * cod**0.5 * re**-2.5))
        else:
            return (1.37e-11 * cod**0.5 * re**-2.5)

    def get_variable(self, product, sds):
        '''Get data from a specific product e.g. '06_L2' '''
        data = readin_MODIS_L2(
            self.product_expand(product),
            self.year,
            self.doy,
            sds=sds,
            times=[self.timestr()], col=self.col)
        return data
        
    def get_band_radiance(self, band, col=None, refl=False, bowtie_corr=True):
        ''' 
        --Bowtie Corrected-- => bowtie_corr=True (default)

        Return the radiance for a specific band from a specific granule.
        Set 'refl' to True to get the reflectance instead (if supported)'''
        if not col:
            col = self.col
        var, ind, metadata = self.get_metadata_band(band, col=col)
        mod_data = readin(
            self.product_expand('021KM'),
            self.year,
            self.doy,
            [var],
            times=[self.timestr()], varind=ind, col=col)[var]
        if bowtie_corr:
            if refl:
                return bowtie_correct((
                    (mod_data.astype('float')-metadata['reflectance_offsets']) *
                    metadata['reflectance_scales'])[:, :-4])
            return bowtie_correct((
                (mod_data.astype('float')-metadata['radiance_offsets']) *
                metadata['radiance_scales'])[:, :-4])
        else:
            if refl:
                return (
                    (mod_data.astype('float')-metadata['reflectance_offsets']) *
                    metadata['reflectance_scales'])[:, :-4]
            return (
                (mod_data.astype('float')-metadata['radiance_offsets']) *
                metadata['radiance_scales'])[:, :-4]

        raise ValueError('Not a valid band')

    def get_band_bt(self, band, col=None):
        if not col:
            col = self.col
        return self._rad_to_bt(
            self.get_band_radiance(band, col=col),
            band_centres(band))

    def _rad_to_bt(self, radiance, wavelength):
        wvl = wavelength*1e-6  # Convert to m
        B = radiance*1e6  # Convert to m-1
        temp1 = ((wvl**5 * B)/(2*scipy.constants.h*scipy.constants.c**2))
        temp2 = np.log((1/(temp1) + 1))
        return (scipy.constants.h*scipy.constants.c/
                (scipy.constants.k*wvl*temp2))

    def product_expand(self, product):
        if not product.startswith('M'):
            # Add on the correct satellite prefix
            product = 'M{}D{}'.format(
                {'T': 'O', 'A': 'Y'}[self.sat],
                product)
        return product

    def download(self, product, col=None, force_redownload=False):
        if not col:
            col = self.col
        if (not self.check(product, col=col)) or force_redownload:
            download(self.product_expand(product),
                     self.year, self.doy,
                     times=[self.timestr()], col=col, force_redownload=force_redownload)

    def check(self, product, col=None):
        return check(self.product_expand(product),
                     self.year, self.doy,
                     time=self.timestr(), col=col)

    def get_metadata_band(self, band, col=None):
        if not col:
            col = self.col
        metadata = readin_metadata(
            self.product_expand('021KM'),
            self.year,
            self.doy,
            self.timestr(),
            col,
            ['EV_1KM_RefSB',
             'EV_1KM_Emissive',
             'EV_250_Aggr1km_RefSB',
             'EV_500_Aggr1km_RefSB'])
        for var in metadata.keys():
            try:
                ind = metadata[var]['band_names'].split(',').index(str(band))
                output = metadata[var]
                output['radiance_offsets'] = output['radiance_offsets'][ind]
                output['radiance_scales'] = output['radiance_scales'][ind]
                try:
                    output['reflectance_offsets'] = output['reflectance_offsets'][ind]
                    output['reflectance_scales'] = output['reflectance_scales'][ind]
                    output['corrected_counts_offsets'] = output['corrected_counts_offsets'][ind]
                    output['corrected_counts_scales'] = output['corrected_counts_scales'][ind]
                except IndexError:
                    pass
                return var, ind, output
            except:
                pass
        raise ValueError('Not a valid band')

    def __repr__(self):
        return 'Granule: '+self.astext()

    def increment(self, number=1):
        '''Step forward 'number' granules'''
        dt = self.datetime()
        dt += timedelta(minutes=5*number)
        _, doy = misc.time.date_to_doy(dt.year, dt.month, dt.day)
        ntime = '{:0>2}{:0>2}'.format(dt.hour, dt.minute)
        return Granule(dt.year, doy, ntime, self.sat, col=self.col)


class _MODISlocator():
    ''' Converts lat-lon to MODIS grid locations for a specified accuracy (factor)'''
    def __init__(self, granule, factor=2, rectified=False, col=None, product=None, dateline=False):
        if not col:
            col = granule.col
        if not product:
            product = '021KM'    
        self.factor = factor

        # Store the granule this locator references
        self.granule_name = granule.astext()  
        data = readin_MODIS_L2(
                granule.product_expand(product),
                granule.year,
                granule.doy,
                sds=[],
                col=col,
                times=[granule.time])
            
        lat = field_interpolate(data['Latitude'])[:, :-5]
        if dateline:
            # Ignores the poles
            if (data['Longitude'].max()>160) and (data['Longitude'].min()<-160):
                data['Longitude'] = np.where(data['Longitude']<0,
                                             data['Longitude']+360,
                                             data['Longitude'])
        lon = field_interpolate(data['Longitude'])[:, :-5]

        if rectified:
            file_correct = locator.search(
                'MODIS', 'ctrect', res='1km')[0]
            with Dataset(file_correct) as ncdf:
                self.rect_conv = ncdf.variables['remap_locate'][:]
            self.rectified = True
        else:
            self.rectified = False
        
        # FUTURE: Transform coordinates then use cKDTree for speed
        self.locator_tree = BallTree(np.array(list(zip(
            np.deg2rad(lon[::factor, ::factor]).ravel(),
            np.deg2rad(lat[::factor, ::factor]).ravel()))), metric='haversine')
        self.shape = lat[::factor, ::factor].shape

    def unraveler(self, loc_ind):
        return [(self.factor * (ind[0] // self.shape[1]),
                 self.factor * (ind[0] % self.shape[1]))
                for ind in loc_ind]

    def locate(self, locs, filter_outside=10, remove_outside=False):
        '''Get the MODIS pixel locations.  Locations further than
        filter_outside km from the nearest pixel return np.nan

        Currently returns nans if passed nan'''
        if not np.all(np.isfinite(locs)):
            return [(np.nan, np.nan)]*len(locs)
        loc_ind = self.locator_tree.query(np.deg2rad(locs))
        if np.min(loc_ind[0] * 6378) > filter_outside:
            return [(np.nan, np.nan)]
        output = np.array(self.unraveler(loc_ind[1]))
        if remove_outside:
            os_points = np.where(loc_ind[0]*6378 > filter_outside)
            # Cannot use nan here as index array
            output[os_points[0], :] = -1  
        if self.rectified:
            output[:, 1] = self.rect_conv[output[:, 1]]
        return output

    def points_in_radius(self, loc, dist):
        '''Return the indexs of all the points within dist km of loc'''
        loc_ind = self.locator_tree.query_radius(np.deg2rad(loc), r=dist/6378)
        return np.array(self.unraveler(loc_ind[0].reshape(-1, 1)))

    def coordinate_shift(self, lon, lat, centre_lon, centre_lat,
                         travel_lon, travel_lat):
        '''Shift a set of coordinates such that the equator passes through the middle of the
        granule.  Then rotate the satellite direction of travel so that north is up'''
        rlon, rlat = self.coordinate_rotate(lon, lat,
                                            centre_lon, centre_lat)
        tlon, tlat = self.coordinate_rotate(travel_lon, travel_lat,
                                            centre_lon, centre_lat)

        # Shift around the satellite direction of travel
        f_angle = np.arctan2(tlon, tlat)

        return (rlon * np.cos(f_angle) - rlat * np.sin(f_angle),
                rlon * np.sin(f_angle) + rlat * np.cos(f_angle))

    def coordinate_rotate(self, lon, lat, centre_lon, centre_lat):
        '''Rotates a set of spherical coordinates so that the centre location
        is on the equator'''
        icentre_lon = np.deg2rad(centre_lon)
        icentre_lat = np.deg2rad(centre_lat)
        ilon = np.deg2rad(lon)
        ilat = np.deg2rad(lat)
        alpha = ilon - icentre_lon
        x1 = np.cos(ilat) * np.cos(alpha)
        x2 = np.cos(ilat) * np.sin(alpha)
        x3 = np.sin(ilat)

        x1p = np.cos(icentre_lat) * x1 + np.sin(icentre_lat) * x3
        x2p = x2
        x3p = -np.sin(icentre_lat) * x1 + np.cos(icentre_lat) * x3

        return np.rad2deg(np.arctan2(x2p,
                                     x1p)), np.rad2deg(np.arcsin(x3p))
