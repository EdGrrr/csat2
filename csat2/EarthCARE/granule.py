from datetime import datetime, timedelta
from csat2.EarthCARE.readfiles import readin_earthcare_curtain_filename, get_orbit_datetimes
from csat2.EarthCARE.download import download, download_file_locations, open_maap_stream
from csat2.EarthCARE.utils import DEFAULT_BASELINE, frame_names, lonlat_vars
import os
import numpy as np
from csat2 import misc, locator


class Granule(object):
    """EarthCARE granules are defined by the orbit and frame numbers"""
    def __init__(self, orbit, frame, stream=False):
        self.orbit = orbit                  # int
        self.frame = frame            # str like 'A'
        # These are created empty and filled as needed. There is probably
        # a more pythonic way to do this. Note that the orbit/frame pair
        # defines the the Granule
        self.dtime = None
        # LonLat data
        self.lonlat = None
        self.lonlat_product = None
        # Should we use the MAAP streaming api? This skips any download
        # or disk checking
        self.stream = stream

    def astext(self):
        return f"EC.{self.orbit:05d}{self.frame}"   # "EC" stands for EarthCARE.

    def __repr__(self):
        return self.astext()

    @classmethod
    def fromtext(cls, gran_text):
        """Create a Granule from a string like 'EC.04606A'."""
        text = gran_text.split(".")[1]
        orbit = int(text[:5])
        frame = text[5:]

        if not frame:
            raise ValueError("Granule text must include a frame, e.g., 'EC.04606A'")

        return cls(orbit, frame)

    @classmethod
    def fromfilename(cls, filename):
        """Create a Granule from an EarthCARE filename."""
        basename = os.path.basename(filename)
        parts = basename.split("_")

        orbit_with_id = os.path.splitext(parts[-1])[0]  # e.g., 04606A
        orbit = int(orbit_with_id[:5])
        frame = orbit_with_id[5:]

        return cls(orbit, frame)

    @classmethod
    def from_datetime(cls, dtime, baseline=DEFAULT_BASELINE):
        """
        Create a Granule from a datetime, by finding the closest matching file.
        Queries ATLID 1B as this should exist for all orbits.

        Currently queries ESA server, so slow and requires network access.
        """

        valid_filenames = download_file_locations(
            'ATL_NOM_1B',
            dtime=dtime,
            baseline=baseline)
        return cls.fromfilename(valid_filenames[0]['id'])

    def datetime(self, baseline=DEFAULT_BASELINE):
        """
        Return the datetime for this granule (specific orbit + frame).
        This should not depend on having the file
        downloaded, so will come from the GEOMETA files.
        Currently queries ESA server, so slow and requires network access.
        """

        if self.dtime is None:
            self.dtime = get_orbit_datetimes(
                'ATL_NOM_1B',
                orbit=self.orbit,
                frame=self.frame,
                baseline=baseline)
        return self.dtime

    def download(self, product, baseline=DEFAULT_BASELINE, force_clean=False, force_redownload=False):
        '''We can download a file based on the orbit number, as we will make a query to the
        ESA server anyway. However, we are not current setup to check in advance if we
        need to without the GEOMETA files.'''

        if self.stream:
            return
        else:
            download(product, orbit=self.orbit, frames=self.frame,
                     baselines=baseline, force_clean=force_clean, force_redownload=force_redownload)
        
    def get_variable(self, product, sds, baseline=DEFAULT_BASELINE):
        """
        Retrieve variables from the EarthCARE curtain file.

        Args:
            product (str, optional): Override default product type.
            sds (str or list of str): Name(s) of variables (Scientific Data Sets)
                 to extract from the file (e.g., "latitude", "longitude").
            baseline (str, optional): Override default baseline.

        Returns:
            dict: variable_name -> numpy array
        """

        if self.stream:
            return open_maap_stream(
                product,
                orbit=self.orbit,
                frame=self.frame,
                baseline=baseline)
        else:
            return readin_earthcare_curtain_filename(
                self.get_filename(product=product, baseline=baseline),
                sds=sds
            )

    def get_stream_location(self, product, baseline=DEFAULT_BASELINE):

        dtime = self.dtime or self.datetime()
        valid_filenames = download_file_locations(product, dtime=dtime, baseline=baseline)
        return valid_filenames[0]['maap_h5']

    def get_filename(self, product, baseline=DEFAULT_BASELINE):
        
        dtime = self.dtime or self.datetime()
        files = locator.search(
            "EarthCARE", product,
            year=dtime.year,
            doy=misc.time.datetime_to_ydh(dtime)[1],
            orbit=self.orbit,
            frame=self.frame,
            baseline=baseline,
        )
        if len(files) == 0:
            raise FileNotFoundError(f'No file for {self} - {product}')
        return files[0]

    def get_lonlat(self, product, baseline=DEFAULT_BASELINE):
        """Return longitude and latitude arrays for this granule."""

        if (self.lonlat is None) or (product != self.lonlat_product):
            data = self.get_variable(
                product=product,
                sds=lonlat_vars[product],
                baseline=baseline
            )
            self.lonlat = (data[lonlat_vars[product][0]],
                           data[lonlat_vars[product][1]])
            self.lonlat_product = product
        return self.lonlat

    def get_decimal_times(self, product, baseline=DEFAULT_BASELINE):
        """
        Return time as decimal hours since midnight UTC on granule date (CloudSat-style).
        
        EarthCARE stores 'ScienceData/time' as seconds since 2000-01-01 00:00:00.
        This function mimics CloudSat's UTC_start + Profile_time format.
        """
        
        data = self.get_variable(product=product, sds=["time"], baseline=baseline)
        time_seconds = data["time"]

        # Reference time: 2000-01-01 00:00:00
        base_datetime = datetime(2000, 1, 1)

        # Convert first timestamp to find granule day
        first_time = base_datetime + timedelta(seconds=float(time_seconds.values[0]))
        start_of_day = datetime(first_time.year, first_time.month, first_time.day)

        # Convert all to decimal hours since midnight
        delta = (time_seconds - (start_of_day - base_datetime).total_seconds()) / 3600.0
        return delta

    def get_datetimes(self, product, baseline=DEFAULT_BASELINE):
        """
        Return array of Python datetime objects for this granule.

        EarthCARE stores time as seconds since 2000-01-01 00:00:00.
        """

        data = self.get_variable(product=product, sds=["time"], baseline=baseline)
        time_seconds = data["time"]

        base_datetime = datetime(2000, 1, 1)
        return [base_datetime + timedelta(seconds=round(float(s))) for s in time_seconds]

    def locate(self, product, locs, baseline=DEFAULT_BASELINE):
        """
        Locate the nearest profile index in the granule to each (lon, lat) point.

        Args:
            product, baseline: Optional overrides.
            locs (list of tuples): List of (longitude, latitude) pairs.

        Returns:
            np.ndarray of int indices, one per location.
        """

        lon, lat = self.get_lonlat(product, baseline)
        return np.array([
            np.argmin(misc.geo.haversine(lon0, lat0, lon.values, lat.values))
            for lon0, lat0 in locs
        ])

    def geolocate(self, product, indicies, baseline=DEFAULT_BASELINE):
        """
        Return the lon/lat of an array of indicies
        """
        lon, lat = self.get_lonlat(product, baseline)
        return  lon.values[indicies], lat.values[indicies]

    def increment(self, number=1):
        current_index = frame_names.index(self.frame)
        
        total_index = current_index + number
        new_orbit = self.orbit + total_index // 8
        new_frame = frame_names[total_index % 8]
        
        return Granule(new_orbit, new_frame)

    def next(self, number=1):
        return self.increment(number)
