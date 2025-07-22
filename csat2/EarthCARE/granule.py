from datetime import datetime, timedelta
from csat2.EarthCARE.readfiles import readin_earthcare_curtain_filename
from csat2.EarthCARE.download import download, download_file_locations, check
from csat2.EarthCARE.utils import DEFAULT_VERSION, frame_names, lonlat_vars
import os
import numpy as np
from csat2 import misc, locator


class Granule(object):
    """EarthCARE granules are defined by the orbit and frame numbers"""
    def __init__(self, orbit, frame, version=DEFAULT_VERSION):
        self.orbit = orbit                  # int
        self.frame = frame            # str like 'A'
        # Keep the version here for now, but given the differences in the product
        # version between different files (this is less consistent than MODIS,
        # it might not make sense to have a product version here.
        self.version = version
        # These are created empty and filled as needed. There is probably
        # a more pythonic way to do this. Note that the orbit/frame pair
        # defines the the Granule
        self.dtime = None
        # LonLat data
        self.lonlat = None

    def astext(self):
        return f"EC.{self.orbit:05d}{self.frame}"   # "EC" stands for EarthCARE.

    def __repr__(self):
        return self.astext()

    @classmethod
    def fromtext(cls, gran_text, version=DEFAULT_VERSION):
        """Create a Granule from a string like 'EC.04606A'."""
        text = gran_text.split(".")[1]
        orbit = int(text[:5])
        frame = text[5:]

        if not frame:
            raise ValueError("Granule text must include a frame, e.g., 'EC.04606A'")

        return cls(orbit, frame, version=version)

    @classmethod
    def fromfilename(cls, filename):
        """Create a Granule from an EarthCARE filename."""
        basename = os.path.basename(filename)
        parts = basename.split("_")

        orbit_with_id = os.path.splitext(parts[-1])[0]  # e.g., 04606A
        orbit = int(orbit_with_id[:5])
        frame = orbit_with_id[5:]

        version = parts[1][-2:]

        return cls(orbit, frame, version=version)

    @classmethod
    def from_datetime(cls, dtime, version=DEFAULT_VERSION):
        """
        Create a Granule from a datetime, by finding the closest matching file.
        Queries ATLID 1B as this should exist for all orbits.

        Currently queries ESA server, so slow and requires network access.
        """
        valid_filenames = download_file_locations('ATL_NOM_1B', dtime=dtime)
        return cls.fromfilename(valid_filenames[0])

    def datetime(self):
        """
        Return the datetime for this granule (specific orbit +
        orbit ID). This should not depend on having the file
        downloaded, so will come from the GEOMETA files.

        Currently queries ESA server, so slow and requires network access.
        """
        if self.dtime is None:
            valid_filenames = download_file_locations('ATL_NOM_1B', orbit=self.orbit, frame=self.frame)
            self.dtime = datetime.strptime(valid_filenames[0].split('_')[5],
                                           '%Y%m%dT%H%M%SZ')
        return self.dtime

    def download(self, product, version=DEFAULT_VERSION, force_redownload=False):
        '''We can download a file based on the orbit number, as we will make a query to the
        ESA server anyway. However, we are not current setup to check in advance if we
        need to without the GEOMETA files.'''
        download(product, orbit=self.orbit, frame=self.frame,
                 version=version, force_redownload=force_redownload)

    def get_variable(self, product, sds, version=None):
        """
        Retrieve variables from the EarthCARE curtain file.

        Args:
            varnames (list of str): Full HDF5 paths like 'ScienceData/longitude'.
            product (str, optional): Override default product type.
            version (str, optional): Override default version.

        Returns:
            dict: variable_name -> numpy array
        """
        if version is None:
            version = self.version

        return readin_earthcare_curtain_filename(
            self.get_filename(product=product, version=version),
            sds=sds
        )

    def get_filename(self, product, version=None):
        if version is None:
            version = self.version

        dtime = self.datetime()
        files = locator.search(
            "EarthCARE", product,
            year=dtime.year,
            doy=misc.time.datetime_to_ydh(dtime)[1],
            orbit=self.orbit,
            frame=self.frame,
            version=version,
        )
        if len(files) == 0:
            raise FileNotFoundError(f'No file for {self} - {product}')
        return files[0]

    
    def get_lonlat(self, product, version=None):
        """Return longitude and latitude arrays for this granule."""
        if version is None:
            version = self.version

        if (self.lonlat is None) or (product != self.lonlat_product):
            data = self.get_variable(
                product=product,
                sds=lonlat_vars[product],
                version=version
            )
            self.lonlat = (data[lonlat_vars[product][0]],
                           data[lonlat_vars[product][1]])
            self.lonlat_product = product
        return self.lonlat

    def get_decimal_times(self, product, version=None):
        """
        Return time as decimal hours since midnight UTC on granule date (CloudSat-style).
        
        EarthCARE stores 'ScienceData/time' as seconds since 2000-01-01 00:00:00.
        This function mimics CloudSat's UTC_start + Profile_time format.
        """
        if version is None:
            version = self.version

        data = self.get_variable(["ScienceData/time"], product=product, version=version)
        time_seconds = data["ScienceData/time"]

        # Reference time: 2000-01-01 00:00:00
        base_datetime = datetime(2000, 1, 1)

        # Convert first timestamp to find granule day
        first_time = base_datetime + timedelta(seconds=float(time_seconds[0]))
        start_of_day = datetime(first_time.year, first_time.month, first_time.day)

        # Convert all to decimal hours since midnight
        delta = (time_seconds - (start_of_day - base_datetime).total_seconds()) / 3600.0
        return delta

    def get_datetimes(self, product, version=None):
        """
        Return array of Python datetime objects for this granule.

        EarthCARE stores time as seconds since 2000-01-01 00:00:00.
        """
        if version is None:
            version = self.version

        data = self.get_variable(["ScienceData/time"], product=product, version=version)
        time_seconds = data["ScienceData/time"]

        base_datetime = datetime(2000, 1, 1)
        return [base_datetime + timedelta(seconds=round(float(s))) for s in time_seconds]

    def locate(self, product, locs, version=None):
        """
        Locate the nearest profile index in the granule to each (lon, lat) point.

        Args:
            product, version: Optional overrides.
            locs (list of tuples): List of (longitude, latitude) pairs.

        Returns:
            np.ndarray of int indices, one per location.
        """
        version = version or self.version

        lon, lat = self.get_lonlat(product, version)
        return np.array([
            np.argmin(misc.geo.haversine(lon0, lat0, lon, lat))
            for lon0, lat0 in locs
        ])

    def geolocate(self, product, indicies, version=None):
        """
        Return the lon/lat of an array of indicies
        """
        version = version or self.version
        lon, lat = self.get_lonlat(product, version)
        return  lon[indicies], lat[indicies]

    def increment(self, number=1):
        current_index = frame_names.index(self.frame)
        
        total_index = current_index + number
        new_orbit = self.orbit + total_index // 8
        new_frame = frame_names[total_index % 8]
        
        return Granule(new_orbit, new_frame,
                       version=self.version)

    def next(self, number=1):
        return self.increment(number)
