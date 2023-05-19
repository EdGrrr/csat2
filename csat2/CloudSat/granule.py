from datetime import datetime, timedelta
from csat2 import misc
from csat2 import locator
from .readfiles import (
    readin_cloudsat_curtain,
    get_orbit_by_time,
    get_orbit_date,
    get_orbit_datetime,
    DEFAULT_COLLECTION,
)
from .download import download, check
import numpy as np
import os


class Granule(object):
    """Cloudsat granules are defined by the orbit number"""

    def __init__(self, orbit: int, col=DEFAULT_COLLECTION):
        """orbit - int"""
        super(Granule, self).__init__()
        self.orbit = orbit
        self.locator = None
        self.lonlat = None
        self.col = col

    @classmethod
    def fromtext(cls, gran_text, col=DEFAULT_COLLECTION):
        """CS.10056"""
        return cls(int(gran_text[3:8]), col)

    @classmethod
    def fromfilename(cls, filename, col=DEFAULT_COLLECTION):
        basename = os.path.basename(filename)
        parts = basename.split("_")
        orbit = int(parts[1])
        col = [p for p in parts[5:7] if p.startswith("P") or p.startswith("R")]
        col = "_".join(col)
        return cls(orbit, col)

    @classmethod
    def from_datetime(cls, dtime, col=DEFAULT_COLLECTION):
        orbit = get_orbit_by_time(dtime)
        return cls(orbit, col)

    def datetime(self):
        return get_orbit_datetime(self.orbit)

    def astext(self):
        return "CS.{:0>5}".format(self.orbit)

    def get_lonlat(self, product, col=None):
        """Get lon lat data - can specify the product or the collection
        used to the geolocation data from"""
        if not self.lonlat:
            self._read_lonlat(product=product, col=col)
        return self.lonlat

    def _read_lonlat(self, product, col=None):
        if not col:
            col = self.col
        data = self.get_variable(product, varnames=["Longitude", "Latitude"], col=col)
        self.lonlat = [data["Longitude"], data["Latitude"]]

    def get_decimal_times(self, product, col=None):
        if not col:
            col = self.col
        data = self.get_variable(
            product, varnames=["Profile_time", "UTC_start"], col=col
        )
        times = data["UTC_start"] + data["Profile_time"]
        return times / 3600

    def get_datetimes(self, product, col=None):
        if not col:
            col = self.col
        data = self.get_variable(
            product, varnames=["Profile_time", "UTC_start"], col=col
        )
        times = data["UTC_start"] + data["Profile_time"]
        raise NotImplementedError()

    def get_variable(self, product, varnames, col=None):
        if not col:
            col = self.col
        return readin_cloudsat_curtain(product, sds=varnames, col=col, orbit=self.orbit)

    def download_product(self, product, col=None):
        if col is None:
            col = self.col
        year, doy = get_orbit_date(self.orbit)

        if not check(product, year, doy, orbit=self.orbit, col=col):
            download(product, year, doy, orbits=[self.orbit], col=col)

    def locate(self, locs, product, col=None):
        if col is None:
            col = self.col
        cslon, cslat = self.get_lonlat(product, col=col)

        return np.array(
            [
                np.argmin(misc.geo.haversine(loc[0], loc[1], cslon, cslat))
                for loc in locs
            ]
        )

    def __repr__(self):
        return self.astext()

    def increment(self, number=1):
        return Granule(self.orbit + number, self.col)
