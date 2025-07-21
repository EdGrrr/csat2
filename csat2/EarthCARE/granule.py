from datetime import datetime, timedelta
from .readfiles import(
    readin_earthcare_curtain,
    get_orbit_by_time,
    get_orbit_date,
    get_orbit_datetimes
)
from .download import download, check
from .utils import get_orbit_date_approx, get_orbit_number_approx
from .utils import DEFAULT_BASELINE, DEFAULT_PRODUCT_TYPE
import os
import numpy as np
from csat2 import misc, locator



class Granule(object):
    """EarthCARE granules are defined by the orbit number."""
    def __init__(self, orbit, orbit_id, product_type="CPR_CLD_2A", baseline="AB"):
        self.orbit = orbit                  # int
        self.orbit_id = orbit_id            # str like 'A'
        self.product_type = product_type
        self.baseline = baseline
        self.orbit_with_id = f"{orbit:05d}{orbit_id}"  


    def astext(self):
        return f"EC.{self.orbit:05d}{self.orbit_id}"   #"EC" stands for EarthCARE.

    def __repr__(self):
        return self.astext()


    @classmethod
    def fromtext(cls, gran_text, product_type=DEFAULT_PRODUCT_TYPE, baseline=DEFAULT_BASELINE):
        """Create a Granule from a string like 'EC.04606A'."""
        text = gran_text.split(".")[1]
        orbit_number = int(text[:5])
        orbit_id = text[5:]

        if not orbit_id:
            raise ValueError("Granule text must include orbit ID, e.g., 'EC.04606A'")

        return cls(orbit_number, orbit_id, product_type=product_type, baseline=baseline)

    @classmethod
    def fromfilename(cls, filename):
        """Create a Granule from an EarthCARE filename."""
        basename = os.path.basename(filename)
        parts = basename.split("_")

        orbit_with_id = os.path.splitext(parts[-1])[0]  # e.g., 04606A
        orbit_number = int(orbit_with_id[:5])
        orbit_id = orbit_with_id[5:]

        product_type = "_".join(parts[2:5])
        baseline = parts[1][-2:]

        return cls(orbit_number, orbit_id, product_type=product_type, baseline=baseline)




    def datetime(self, product_type=None, baseline=None):
        """
        Return the datetime for this granule (specific orbit + orbit ID).
        """
        product_type = product_type or self.product_type
        baseline = baseline or self.baseline

        return get_orbit_datetimes(
            orbit_number=self.orbit,
            orbit_id=self.orbit_id,
            product_type=product_type,
            baseline=baseline
        )


    def download_product(self, product_type=None, baseline=None):
        """
        Download all EarthCARE granules on the estimated day of this granule.
        Existing files will be skipped automatically.
        """
        product_type = product_type or self.product_type
        baseline = baseline or self.baseline

        candidates = get_orbit_date(self.orbit)
        if not candidates:
            raise ValueError(f"No date estimate found for orbit {self.orbit}")

        for year, doy in candidates:
            date_obj = datetime.strptime(f"{year} {doy}", "%Y %j")
            month = date_obj.month
            day = date_obj.day

            print(f"[INFO] Attempting download for {year}-{month:02d}-{day:02d}")
            download(product_type, baseline, year=year, month=month, day=day)
            return  # stop after first candidate





    def get_variable(self, varnames, product_type=None, baseline=None):
        """
        Retrieve variables from the EarthCARE curtain file.

        Args:
            varnames (list of str): Full HDF5 paths like 'ScienceData/longitude'.
            product_type (str, optional): Override default product type.
            baseline (str, optional): Override default baseline.

        Returns:
            dict: variable_name -> numpy array
        """
        if product_type is None:
            product_type = self.product_type
        if baseline is None:
            baseline = self.baseline

        return readin_earthcare_curtain(
            product_type=product_type,
            baseline=baseline,
            orbit_number=self.orbit,
            orbit_id=self.orbit_id,
            sds=varnames
        )



    def get_lonlat(self, product_type=None, baseline=None):
        """Return longitude and latitude arrays for this granule."""
        if product_type is None:
            product_type = self.product_type
        if baseline is None:
            baseline = self.baseline

        data = self.get_variable(
            ["ScienceData/longitude", "ScienceData/latitude"],
            product_type=product_type,
            baseline=baseline
        )
        return data["ScienceData/longitude"], data["ScienceData/latitude"]



    def get_decimal_times(self, product_type=None, baseline=None):
        """
        Return time as decimal hours since midnight UTC on granule date (CloudSat-style).
        
        EarthCARE stores 'ScienceData/time' as seconds since 2000-01-01 00:00:00.
        This function mimics CloudSat's UTC_start + Profile_time format.
        """
        if product_type is None:
            product_type = self.product_type
        if baseline is None:
            baseline = self.baseline

        data = self.get_variable(["ScienceData/time"], product_type=product_type, baseline=baseline)
        time_seconds = data["ScienceData/time"]

        # Reference time: 2000-01-01 00:00:00
        base_datetime = datetime(2000, 1, 1)

        # Convert first timestamp to find granule day
        first_time = base_datetime + timedelta(seconds=float(time_seconds[0]))
        start_of_day = datetime(first_time.year, first_time.month, first_time.day)

        # Convert all to decimal hours since midnight
        delta = (time_seconds - (start_of_day - base_datetime).total_seconds()) / 3600.0
        return delta


    def get_datetimes(self, product_type=None, baseline=None):
        """
        Return array of Python datetime objects for this granule.

        EarthCARE stores time as seconds since 2000-01-01 00:00:00.
        """
        if product_type is None:
            product_type = self.product_type
        if baseline is None:
            baseline = self.baseline

        data = self.get_variable(["ScienceData/time"], product_type=product_type, baseline=baseline)
        time_seconds = data["ScienceData/time"]

        base_datetime = datetime(2000, 1, 1)
        return [base_datetime + timedelta(seconds=round(float(s))) for s in time_seconds]


    def locate(self, locs, product_type=None, baseline=None):
        """
        Locate the nearest profile index in the granule to each (lon, lat) point.

        Args:
            locs (list of tuples): List of (longitude, latitude) pairs.
            product_type, baseline: Optional overrides.

        Returns:
            np.ndarray of int indices, one per location.
        """
        product_type = product_type or self.product_type
        baseline = baseline or self.baseline

        lon, lat = self.get_lonlat(product_type, baseline)
        return np.array([
            np.argmin(misc.geo.haversine(lon0, lat0, lon, lat))
            for lon0, lat0 in locs
        ])



    # def increment(self, number=1):
    #     """
    #     Return a new Granule with the orbit number incremented.
    #     Note: Orbit ID is preserved.
    #     """
    #     return Granule(self.orbit + number, self.orbit_id,
    #                 product_type=self.product_type,
    #                 baseline=self.baseline)



#    def increment(self, number=1):
#        """
#        Return a new Granule with incremented orbit and/or orbit_id.
#        EarthCARE orbit_ids go from 'A' to 'H' (8 per orbit).
#        """
#        orbit_ids = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
#        current_index = orbit_ids.index(self.orbit_id)
#        
#        total_index = current_index + number
#        new_orbit = self.orbit + total_index // 8
#        new_orbit_id = orbit_ids[total_index % 8]
#
#        return Granule(new_orbit, new_orbit_id,
#                    product_type=self.product_type,
#                    baseline=self.baseline)






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

        orbit_ids = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
        direction = 1 if number > 0 else -1
        skip_left = abs(number)
        scanned   = 0

        cand_orbit = self.orbit
        cand_id    = self.orbit_id

        while scanned < max_scan:
            # ── step one granule ───────────────────────────────────────────
            idx        = orbit_ids.index(cand_id) + direction
            cand_orbit = cand_orbit + idx // 8
            cand_id    = orbit_ids[idx % 8]

            candidate = Granule(
                cand_orbit, cand_id,
                product_type=self.product_type,
                baseline=self.baseline
            )

            # ── which calendar day(s) contain this orbit? ─────────────────
            try:
                for yr, doy in get_orbit_date(cand_orbit):
                    dt = datetime.strptime(f"{yr} {doy}", "%Y %j")
                    exists = check(
                        candidate.product_type, candidate.baseline,
                        year=dt.year, month=dt.month, day=dt.day,
                        orbit=candidate.orbit, orbit_id=candidate.orbit_id
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
    def from_datetime(cls, dtime, product_type=DEFAULT_PRODUCT_TYPE, baseline=DEFAULT_BASELINE):
        """
        Create a Granule from a datetime, by finding the closest matching file.
        """
        orbit_id_full = get_orbit_by_time(dtime, product_type=product_type, baseline=baseline)

        orbit_number = int(orbit_id_full[:5])
        orbit_id = orbit_id_full[5:]

        return cls(orbit_number, orbit_id, product_type=product_type, baseline=baseline)
