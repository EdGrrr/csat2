from datetime import datetime, timedelta
from csat2.EarthCARE.readfiles import (
    readin_earthcare_curtain,
    get_orbit_by_time,
    get_orbit_date,
    get_orbit_datetimes
)
from csat2.EarthCARE.download import download, check
from csat2.EarthCARE.utils import get_orbit_date_approx, get_orbit_number_approx, DEFAULT_VERSION
import os
import numpy as np
from csat2 import misc, locator


class Granule(object):
    """EarthCARE granules are defined by the orbit and frame numbers"""
    def __init__(self, orbit, frame, version=DEFAULT_VERSION):
        self.orbit = orbit                  # int
        self.frame = frame            # str like 'A'
        self.version = version
        # These are created empty and filled as needed. There is probably
        # a more pythonic way to do this. Note that the orbit/frame pair
        # defines the the Granule
        self.year = None
        self.doy = None
        self.time = None

    def astext(self):
        return f"EC.{self.orbit:05d}{self.frame}"   #"EC" stands for EarthCARE.

    def __repr__(self):
        return self.astext()

    @classmethod
    def fromtext(cls, gran_text, version=DEFAULT_VERSION):
        """Create a Granule from a string like 'EC.04606A'."""
        text = gran_text.split(".")[1]
        orbit = int(text[:5])
        frame = text[5:]

        if not frame:
            raise ValueError("Granule text must include orbit ID, e.g., 'EC.04606A'")

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

    def datetime(self):
        """
        Return the datetime for this granule (specific orbit +
        orbit ID). This should not depend on having the file
        downloaded, so will come from the GEOMETA files.
        """
        return get_orbit_datetimes(
            orbit_number=self.orbit,
            frame=self.frame)

    def download(self, product, version=DEFAULT_VERSION):
        '''We can download a file based on the orbit number, as we will make a query to the
        ESA server anyway. However, we are not current setup to check in advance if we
        need to without the GEOMETA files.'''
        download(product, orbit=self.orbit, frame=self.frame, version=version)
        #if not check(product, self.year, self.doy, self.orbit, self.frame, version):
        #    download(product, orbit=self.orbit, frame=self.frame, version=version)

    def download_product(self, product=None, version=None):
        """
        Download all EarthCARE granules on the estimated day of this granule.
        Existing files will be skipped automatically.
        """
        product = product or self.product
        version = version or self.version

        candidates = get_orbit_date(self.orbit)
        if not candidates:
            raise ValueError(f"No date estimate found for orbit {self.orbit}")

        for year, doy in candidates:
            date_obj = datetime.strptime(f"{year} {doy}", "%Y %j")
            month = date_obj.month
            day = date_obj.day

            print(f"[INFO] Attempting download for {year}-{month:02d}-{day:02d}")
            download(product, version, year=year, month=month, day=day)
            return  # stop after first candidate

    def get_variable(self, product, varnames, version=None):
        """
        Retrieve variables from the EarthCARE curtain file.

        Args:
            varnames (list of str): Full HDF5 paths like 'ScienceData/longitude'.
            product (str, optional): Override default product type.
            version (str, optional): Override default version.

        Returns:
            dict: variable_name -> numpy array
        """
        if product is None:
            product = self.product
        if version is None:
            version = self.version

        return readin_earthcare_curtain(
            product=product,
            version=version,
            orbit_number=self.orbit,
            frame=self.frame,
            sds=varnames
        )



    def get_lonlat(self, product=None, version=None):
        """Return longitude and latitude arrays for this granule."""
        if product is None:
            product = self.product
        if version is None:
            version = self.version

        data = self.get_variable(
            ["ScienceData/longitude", "ScienceData/latitude"],
            product=product,
            version=version
        )
        return data["ScienceData/longitude"], data["ScienceData/latitude"]



    def get_decimal_times(self, product=None, version=None):
        """
        Return time as decimal hours since midnight UTC on granule date (CloudSat-style).
        
        EarthCARE stores 'ScienceData/time' as seconds since 2000-01-01 00:00:00.
        This function mimics CloudSat's UTC_start + Profile_time format.
        """
        if product is None:
            product = self.product
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


    def get_datetimes(self, product=None, version=None):
        """
        Return array of Python datetime objects for this granule.

        EarthCARE stores time as seconds since 2000-01-01 00:00:00.
        """
        if product is None:
            product = self.product
        if version is None:
            version = self.version

        data = self.get_variable(["ScienceData/time"], product=product, version=version)
        time_seconds = data["ScienceData/time"]

        base_datetime = datetime(2000, 1, 1)
        return [base_datetime + timedelta(seconds=round(float(s))) for s in time_seconds]


    def locate(self, locs, product=None, version=None):
        """
        Locate the nearest profile index in the granule to each (lon, lat) point.

        Args:
            locs (list of tuples): List of (longitude, latitude) pairs.
            product, version: Optional overrides.

        Returns:
            np.ndarray of int indices, one per location.
        """
        product = product or self.product
        version = version or self.version

        lon, lat = self.get_lonlat(product, version)
        return np.array([
            np.argmin(misc.geo.haversine(lon0, lat0, lon, lat))
            for lon0, lat0 in locs
        ])



    # def increment(self, number=1):
    #     """
    #     Return a new Granule with the orbit number incremented.
    #     Note: Orbit ID is preserved.
    #     """
    #     return Granule(self.orbit + number, self.frame,
    #                 product=self.product,
    #                 version=self.version)



#    def increment(self, number=1):
#        """
#        Return a new Granule with incremented orbit and/or frame.
#        EarthCARE frames go from 'A' to 'H' (8 per orbit).
#        """
#        frames = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
#        current_index = frames.index(self.frame)
#        
#        total_index = current_index + number
#        new_orbit = self.orbit + total_index // 8
#        new_frame = frames[total_index % 8]
#
#        return Granule(new_orbit, new_frame,
#                    product=self.product,
#                    version=self.version)






    def increment(self, number: int = 1, max_scan: int = 240):
        """
        Return the *next available* Granule on disk.

        Parameters
        ----------
        number : int, default 1
            1 → next granule that exists;
            2 → skip one existing granule, etc.
            Negative values walk backwards.
        max_scan : int, default 8
            Maximum number of granules to *inspect* while searching.

        Raises
        ------
        FileNotFoundError
            If no on-disk granule is found within ``max_scan`` inspections.
        """
        if number == 0:
            raise ValueError("number must be non-zero")

        frames = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
        direction = 1 if number > 0 else -1
        skip_left = abs(number)
        scanned   = 0

        cand_orbit = self.orbit
        cand_id    = self.frame

        while scanned < max_scan:
            # ── step one granule ───────────────────────────────────────────
            idx        = frames.index(cand_id) + direction
            cand_orbit = cand_orbit + idx // 8
            cand_id    = frames[idx % 8]

            candidate = Granule(
                cand_orbit, cand_id,
                product=self.product,
                version=self.version
            )

            # ── which calendar day(s) contain this orbit? ─────────────────
            try:
                for yr, doy in get_orbit_date(cand_orbit):
                    dt = datetime.strptime(f"{yr} {doy}", "%Y %j")
                    exists = check(
                        candidate.product, candidate.version,
                        year=dt.year, month=dt.month, day=dt.day,
                        orbit=candidate.orbit, frame=candidate.frame
                    )

                    if not exists:
                        # found the day but granule isn’t downloaded
                        print(f"[skip] {candidate.astext()} – expected but not on disk")
                        continue

                    # granule *is* on disk
                    if skip_left == 1:
                        return candidate          # ← success
                    skip_left -= 1                # skip & keep searching

            except FileNotFoundError as err:
                # no geometa / no files for that orbit on that day
                print(f"[skip] orbit {cand_orbit}: {err}")

            scanned += 1                          # inspected one more candidate

        # after max_scan attempts we found nothing
        raise FileNotFoundError(
            f"No local EarthCARE file found within {max_scan} granules "
            f"starting from {self.astext()}"
        )






    @classmethod
    def from_datetime(cls, dtime, version=DEFAULT_VERSION):
        """
        Create a Granule from a datetime, by finding the closest matching file.
        """
        frame_full = get_orbit_by_time(dtime, product=product, version=version)

        orbit_number = int(frame_full[:5])
        frame = frame_full[5:]

        return cls(orbit_number, frame, product=product, version=version)
