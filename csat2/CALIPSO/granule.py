from datetime import timedelta
from csat2 import misc
from .readfiles import readin_calipso_curtain, get_filename_by_time, filename_to_datetime, DEFAULT_COLLECTION
from .download import download, check
import os


class Granule(object):
    '''CALIPSO granules are defined by time'''

    def __init__(self, year, doy, hour, minute, second, col=DEFAULT_COLLECTION):
        '''orbit - int'''
        super(Granule, self).__init__()
        self.year = year
        self.doy = doy
        self.hour = hour
        self.minute = minute
        self.second = second
        self.locator = None
        self.lonlat = None
        self.col = col

    @classmethod
    def from_filename(cls, filename, col=DEFAULT_COLLECTION):
        dtime = filename_to_datetime(filename)
        _, doy, _ = misc.time.datetime_to_ydh(dtime)
        return cls(dtime.year, doy, dtime.hour, dtime.minute, dtime.second, col=col)

    @classmethod
    def from_datetime(cls, dtime, col=DEFAULT_COLLECTION):
        filename = get_filename_by_time(dtime)
        return cls.from_filename(filename, col)

    def filename(self, product, col=DEFAULT_COLLECTION):
        return get_filename_by_time(self.datetime()+timedelta(minutes=1))
    
    def datetime(self):
        return misc.time.ydh_to_datetime(self.year, self.doy, self.hour+self.minute/60+self.second/3600)
    
    def astext(self):
        return f'CAL.{self.year}{self.doy:0>3}.{self.hour:0>2}{self.minute:0>2}{self.second:0>2}'

    def get_lonlat(self, product, col=None):
        '''Get lon lat data - can specify the product or the collection
        used to the geolocation data from'''
        if not self.lonlat:
            self._read_lonlat(product=product, col=col)
        return self.lonlat
    
    def _read_lonlat(self, product, col=None):
        if not col:
            col = self.col
        data = self.get_variable(
            product,
            varnames=['Longitude', 'Latitude'],
            col=col)
        self.lonlat = [data['Longitude'],
                       data['Latitude']]

    def get_variable(self, product, varnames, col=None):
        if not col:
            col = self.col
        return readin_calipso_curtain(
            product, self.year, self.doy, self.hour, self.minute, sds=varnames, col=col)

    def download_product(self, product, col=None):
        if col is None:
            col = self.col
        if not check(product, self.year, self.doy,
                     self.hour, self.minute,
                     col=col):
            fname = self.filename('LID_L1')
            download(product, self.year, self.doy,
                     self.hour, self.minute,
                     self.second, fname[-4],
                     col=col)

    def __repr__(self):
        return self.astext()

    def increment(self):
        return Granule.from_datetime(self.datetime()+timedelta(hours=1))
