"""Provides simple functions for using GOES data"""
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
from csat2 import misc
import csat2.misc.geo
from .util import get_resolution as util_get_resolution
import logging
import xarray as xr
import warnings

log = logging.getLogger(__name__)


def readin_radiances_filename(filename):
    with xr.open_dataset(filename) as ds:
        return ds["Rad"][:]


def available_files(
    product, year, doy, hour, minute, sat, channel, area="RadF", mode=6
):
    if minute == "*":
        minute = "**"
    if hour == "*":
        hour = "**"
    filenames = locator.search(
        "GOES",
        product,
        year=year,
        doy=doy,
        hour=hour,
        minute=minute,
        channel=channel,
        area=area,
        mode=mode,
        sat=sat,
    )
    return filenames


class Granule:
    """A class for handling GOES images. Will increment through images
    and should be able to cope with most common location transformation
    tasks (including parallax correction. Will also download files,
    assuming you place an appropriate API key in your csat2 config
    directory.

    Note, the different scan areas are handled separately, it is
    assumed you wont be converting between them."""

    inc_minutes = {"RadM1": 1, "RadM2": 1, "RadC": 5, "RadF": 10}
    product_names = {
        "L2-CPS": "PSD",
        "L2-COD": "COD",
        "L2-AOD": "AOD",
        "L2-ACTP": "Phase",
    }

    def __init__(self, sat, area, year, doy, hour, minute, locator=None):
        self.sat = sat
        self.area = area
        self.year = year
        self.doy = doy
        self.hour = hour
        # Minute is stored as the start time for that block
        #  Use the inc_minutes value to shift to the start of the block
        self.minute = minute - (minute % self.inc_minutes[area])
        self.locator = locator

    @classmethod
    def fromtext(cls, gran_text):
        m = re.search(
            "G(?P<sat>..).(?P<year>....)(?P<doy>...).(?P<hour>..)(?P<minute>..).(?P<area>.*)",
            gran_text,
        )
        return cls(
            sat="G" + m.group("sat"),
            area=m.group("area"),
            year=int(m.group("year")),
            doy=int(m.group("doy")),
            hour=int(m.group("hour")),
            minute=int(m.group("minute")),
        )

    @classmethod
    def fromfilename(cls, filename):
        filename = os.path.basename(filename)
        m = re.search(
            "OR_ABI-(.*?)-(?P<area>.*?)-M(?P<mode>.)C(?P<channel>..)"
            + "_G(?P<sat>..)_s(?P<year>....)(?P<doy>...)(?P<hour>..)(?P<minute>..)",
            filename,
        )
        return cls(
            sat="G" + m.group("sat"),
            area=m.group("area"),
            year=int(m.group("year")),
            doy=int(m.group("doy")),
            hour=int(m.group("hour")),
            minute=int(m.group("minute")),
        )

    @classmethod
    def fromdatetime(cls, sat, area, dtime):
        return cls(
            sat=sat,
            area=area,
            year=dtime.year,
            doy=dtime.timetuple().tm_yday,
            hour=dtime.hour,
            minute=dtime.minute,
        )

    def astext(self):
        return "{}.{}{:0>3}.{:0>2}{:0>2}.{}".format(
            self.sat, self.year, self.doy, self.hour, self.minute, self.area
        )

    def __repr__(self):
        return self.astext()

    def datetime(self):
        return datetime.datetime(
            year=self.year, month=1, day=1, hour=self.hour, minute=self.minute
        ) + datetime.timedelta(days=self.doy - 1)

    def get_resolution(self, channel, product="L1b", mode="*", fromfile=False):
        if fromfile:
            with Dataset(
                self.get_filename(channel, product=product, mode=mode)
            ) as ncdf:
                return ncdf.variables["Rad"].resolution.split(" ")[1]
        else:
            return util_get_resolution(channel) * 0.000028

    def get_resolution_km(self, channel, product="L1b"):
        return util_get_resolution(channel)

    def _check_locator(self, channel, product="L1b", mode="*"):
        if not self.locator:
            self.locator = GOESLocator(
                self.get_filename(channel=channel, product=product, mode=mode)
            )

    def get_lonlat(self, channel, product="L1b", mode="*"):
        """Proes the lon-lat arrays for a specific channel.
        Note that as the channels have different spatial resolutions,
        ahannel must be specified here."""
        ipc = self.get_index(channel, product, mode)
        llp = self.geolocate(
            np.array(ipc).reshape(2, -1).transpose(), channel, product, mode
        )
        return llp[:, 0].reshape(ipc[0].shape), llp[:, 1].reshape(ipc[0].shape)

    def geolocate(
        self,
        coords,
        channel=None,
        product="L1b",
        mode="*",
        interp=False,
        force_new_locator=False,
    ):
        """Returns the lon/lat of the gridbox indicies passed as coords"""
        if force_new_locator:
            self.locator = GOESLocator(
                self.get_filename(channel=channel, product=product, mode=mode)
            )
        else:
            self._check_locator(channel, product, mode)
        return self.locator.geolocate(coords, interp)

    def locate(
        self,
        coords,
        alt=None,
        channel=None,
        product="L1b",
        mode="*",
        interp=False,
        force_new_locator=False,
    ):
        if force_new_locator:
            self.locator = GOESLocator(
                self.get_filename(channel=channel, product=product, mode=mode)
            )
        else:
            self._check_locator(channel, product, mode)
        return self.locator.locate(coords, alt, interp)

    def points_in_radius(self, loc, dist, channel, product="L1b", mode="*"):
        self._check_locator(channel, product, mode)
        return self.locator.points_in_radius(loc, dist)

    def next(self, number=1, only_downloaded=False, only_exisiting=False):
        """Increment image name"""
        dt = self.datetime()
        dt += datetime.timedelta(minutes=number * self.inc_minutes[self.area])
        year, doy = misc.time.date_to_doy(dt.year, dt.month, dt.day)
        return type(self)(
            self.sat,
            self.area,
            year,
            doy,
            dt.hour,
            dt.minute,
            # Reasonable assumption that coordinates remain the same
            locator=self.locator,
        )

    def get_filename(self, channel=None, product="L1b", mode="*"):
        """This is more complicated that the simple locator, as we
        have to account for the GOES image not being taken at an exact
        time. Even if we can trust the timing, changes in the scan pattern
        can cause a change in the image time.

        Selects an image with a name that puts it within the increment
        timestep."""
        if product == "L2-DMW":
            filenames = locator.search(
                "GOES",
                product,
                year=self.year,
                doy=self.doy,
                hour=self.hour,
                minute=str(int(self.minute) // 10) + "*",
                area=self.area,
                sat=self.sat,
                areas=self.area[-1],
                channel=channel,
                mode=mode,
            )
        elif product.startswith("L2"):
            filenames = locator.search(
                "GOES",
                product,
                year=self.year,
                doy=self.doy,
                hour=self.hour,
                minute=str(int(self.minute) // 10) + "*",
                area=self.area,
                sat=self.sat,
                areas=self.area[-1],
                mode=mode,
            )
        else:
            filenames = locator.search(
                "GOES",
                product,
                year=self.year,
                doy=self.doy,
                hour=self.hour,
                minute=str(int(self.minute) // 10) + "*",
                area=self.area,
                sat=self.sat,
                channel=channel,
                mode=mode,
            )

        minutes = np.array([int(f.split("_")[-3][-5:-3]) for f in filenames])
        ind = np.where(
            np.logical_and(
                ((minutes - self.minute) >= 0),
                ((minutes - self.minute) < self.inc_minutes[self.area]),
            )
        )
        if len(ind[0]) == 0:
            raise IndexError("No matching file")
        if len(ind[0]) > 1:
            raise IndexError("Non-unique filename")
        return filenames[ind[0][0]]

    def get_band_radiance(self, channel, product="L1b", mode="*", refl=False):
        with xr.open_dataset(
            self.get_filename(channel, product=product, mode=mode)
        ) as ds:
            data = ds["Rad"][:]
            if refl:
                data *= ds["kappa0"]
            return data

    def get_band_bt(self, channel, product="L1b", mode="*"):
        """Conversion of radiances to BT using constants in files"""
        with xr.open_dataset(
            self.get_filename(channel, product=product, mode=mode)
        ) as ds:
            data = ds["Rad"][:]
            fk1, fk2 = (ds["planck_fk1"], ds["planck_fk2"])
            bc1, bc2 = (ds["planck_bc1"], ds["planck_bc2"])
            log.info(fk1, fk2, bc1, bc2)
            # TODO: This either throws an xarray warning for an invalid logbook
            # or an error if you try and use np.where to avoid the nans
            # For now, I am going to suppress the warning, which will work as long
            # as the behaviour of nans doesn't change
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                bt = fk2 / (np.log((fk1 / (data)) + 1))
            bt = (bt - bc1) / bc2
            return bt

    def get_product_data(self, product, mode="*", dqf_filter=None):
        """Use for readin in derived data. For DMW, use get_dmw_data.
        dqf_filter is a lambda that determines returns true if the dqf satisfies the reading requirements"""
        with xr.open_dataset(self.get_filename(product=product, mode=mode)) as ds:
            data = ds[self.product_names[product]][:].data
            if dqf_filter is not None:
                data = np.where(dqf_filter(ds["DQF"]), data, np.nan)
            return data

    def get_dmw_data(self, product="L2-DMW", mode="*", channel=None):
        data = {}
        with xr.open_dataset(
            self.get_filename(product=product, mode=mode, channel=channel)
        ) as ds:
            for name in ["lat", "lon", "wind_speed", "wind_direction", "pressure"]:
                data[name] = ds[name][:]
            return data

    def get_shape(self, channel, product="L1b", mode="*"):
        with xr.open_dataset(
            self.get_filename(channel, product=product, mode=mode)
        ) as ds:
            return ds["Rad"].shape

    def get_index(self, channel, product="L1b", mode="*"):
        granshape = self.get_shape(channel, product, mode)
        return np.array(
            np.meshgrid(np.arange(0, granshape[1]), np.arange(0, granshape[0]))
        )

    def get_viewangles(self, channel, product="L1b", mode="*"):
        with xr.open_dataset(
            self.get_filename(channel, product=product, mode=mode)
        ) as ds:
            return np.array(np.meshgrid(ds["x"][:], ds["y"][:]))

    def get_satza(self, channel):
        self._check_locator(channel)
        return np.fromfunction(
            lambda i, j: np.sqrt(
                self.locator.x[i.astype("int")] ** 2
                + self.locator.y[j.astype("int")] ** 2
            ),
            (len(self.locator.x), len(self.locator.y)),
        )

    def get_llcoord(self, channel, product="L1b", mode="*"):
        with xr.open_dataset(
            self.get_filename(channel, product=product, mode=mode)
        ) as ds:
            return ds["x"][0], ds["y"][-1]

    def check(self, channel=None, product="L1b", mode="*"):
        """Is there a filename that satisfies these criteria?"""
        try:
            self.get_filename(channel=channel, product=product, mode=mode)
            return True
        except:
            return False

    def download(
        self,
        channel=None,
        product="L1b",
        mode="*",
        force_redownload=False,
        source="google",
    ):
        if (not force_redownload) and self.check(channel, product, mode):
            log.info("{} Ch:{} exists".format(self.astext(), channel))
            return
        local_folder = locator.get_folder(
            "GOES",
            product,
            year=self.year,
            doy=self.doy,
            hour=self.hour,
            minute=self.minute,
            channel=channel,
            area=self.area,
            mode=mode,
            sat=self.sat,
        )

        # Create the local folder (if required)
        try:
            os.makedirs(local_folder)
        except FileExistsError:
            pass

        dl_minute = str(int(self.minute) // 10) + "*"

        if source == "google":
            bucket = "gcp-public-data-goes-{sat}".format(sat=self.sat[1:])
            if product == "L1b":
                pattern = (
                    "ABI-{product}-{area:.4}/{year}/{doy:0>3}/"
                    + "{hour:0>2}/OR_ABI-{product}-{area}-M{mode}C{channel:0>2}_"
                    + "G{sat}_s{year}{doy:0>3}{hour:0>2}{minute:0>2}*"
                ).format(
                    product=product,
                    year=self.year,
                    doy=self.doy,
                    hour=self.hour,
                    minute=dl_minute,
                    channel=channel,
                    area=self.area,
                    mode=mode,
                    sat=self.sat[1:],
                )
            elif product.startswith("L2"):
                pattern = (
                    "ABI-{product}{area}/{year}/{doy:0>3}/"
                    + "{hour:0>2}/OR_ABI-{product}{area}-M{mode}_"
                    + "G{sat}_s{year}{doy:0>3}{hour:0>2}{minute:0>2}*"
                ).format(
                    product=product,
                    year=self.year,
                    doy=self.doy,
                    hour=self.hour,
                    minute=dl_minute,
                    area=self.area[-1],
                    mode=mode,
                    sat=self.sat[1:],
                )
            prefix = os.path.dirname(pattern)
            log.debug(prefix)

            # You will need to create a credential file for a blank project
            # There are no permissions here, but otherwise this fails as you
            # have no authentication file. It may be possible to do this with
            # a blank file?
            storage_client = storage.Client.from_service_account_json(
                os.path.expandvars("${HOME}/.csat2/goes-service-key.json")
            )
            log.debug(storage_client)
            files = storage_client.get_bucket(bucket).list_blobs(prefix=prefix)

            for fl in files:
                log.debug(fl.name)
                if fnmatch(fl.name, pattern):
                    re_minutes = int(fl.name.split("_")[-3][-5:-3])
                    if ((re_minutes - self.minute) >= 0) and (
                        (re_minutes - self.minute) < self.inc_minutes[self.area]
                    ):
                        log.debug(os.path.basename(fl.name))
                        fl.download_to_filename(
                            os.path.join(local_folder, os.path.basename(fl.name))
                        )

        elif source == "amazon":
            raise NotImplementedError("AWS source not yet implmented")


class GOESLocator:
    """For converting between geodetic and satellite
    coordinates, following
    https://www.goes-r.gov/users/docs/PUG-L1b-vol3.pdf"""

    def __init__(self, filename):
        m = re.search(
            "OR_ABI-(.*?)-(?P<area>.*?)-M(?P<mode>.)C(?P<channel>..)"
            + "_G(?P<sat>..)_s(?P<year>....)(?P<doy>...)(?P<hour>..)(?P<minute>..)",
            filename,
        )
        self.res = util_get_resolution(int(m.group("channel")))
        with Dataset(filename) as ncdf:
            ncdf.set_auto_mask(False)
            gvar = ncdf.variables["goes_imager_projection"]
            self.req = gvar.semi_major_axis
            self.invf = gvar.inverse_flattening
            self.rpol = gvar.semi_minor_axis
            self.e = 0.0818191910435
            self.pph = gvar.perspective_point_height
            self.H = self.pph + self.req
            self.lam0 = np.deg2rad(gvar.longitude_of_projection_origin)
            # note, these values are the bi centres, not edges
            self.x = ncdf.variables["x"][:].astype("float64")
            self.y = ncdf.variables["y"][:].astype("float64")
            # These are the upper/left edges
            self.xl = self.x - (self.x[1] - self.x[0]) / 2
            self.yl = self.y - (self.y[1] - self.y[0]) / 2
            self.xint = None
            self.corrdx = None

    def _create_interp(self):
        # These assume the bin centre for x and y
        self.xint = scipy.interpolate.interp1d(
            range(len(self.x)), self.x, fill_value="extrapolate"
        )
        self.yint = scipy.interpolate.interp1d(
            range(len(self.y)), self.y, fill_value="extrapolate"
        )
        # The 'r' interpolators allow you to return from a position
        # in ABI grid to the index in the image. Will return errval
        # if you are ouside the grid area.
        self.rxint = scipy.interpolate.interp1d(
            self.x, range(len(self.x)), bounds_error=False, fill_value=np.nan
        )
        self.ryint = scipy.interpolate.interp1d(
            self.y, range(len(self.y)), bounds_error=False, fill_value=np.nan
        )

    def _sat_to_geodetic(self, x, y):
        """x, y - the longitude and latitude in the satellite
        frame (in radians)"""
        a = np.sin(x) ** 2 + np.cos(x) ** 2 * (
            np.cos(y) ** 2 + (self.req**2 / self.rpol**2) * np.sin(y) ** 2
        )
        b = -2 * self.H * np.cos(x) * np.cos(y)
        c = self.H**2 - self.req**2
        disc = b**2 - 4 * a * c
        if isinstance(x, float):
            if disc < 0:
                return np.nan, np.nan
            rs = (-b - np.sqrt(disc)) / (2 * a)
        else:
            # Is an array, not single point
            rs = (
                -b - np.sqrt(disc, out=np.empty(disc.shape) * np.nan, where=(disc >= 0))
            ) / (2 * a)
            # rs = np.where(disc >= 0, (-b - np.sqrt(disc))/(2*a), np.nan)
        sx = rs * np.cos(x) * np.cos(y)
        sy = -rs * np.sin(x)
        sz = rs * np.cos(x) * np.sin(y)

        lon = self.lam0 - np.arctan(sy / (self.H - sx))
        lat = np.arctan(
            (self.req**2 * sz)
            / (self.rpol**2 * np.sqrt((self.H - sx) ** 2 + sy**2))
        )
        return lon, lat

    def _geodetic_to_sat(self, lon, lat):
        phic = np.arctan((self.rpol / self.req) ** 2 * np.tan(lat))
        rc = self.rpol / np.sqrt(1 - (self.e * np.cos(phic)) ** 2)
        sx = self.H - rc * np.cos(phic) * np.cos(lon - self.lam0)
        sy = -rc * np.cos(phic) * np.sin(lon - self.lam0)
        sz = rc * np.sin(phic)

        x = np.arcsin(-sy / np.sqrt(sx**2 + sy**2 + sz**2))
        y = np.arctan(sz / sx)
        if isinstance(x, float):
            if np.isnan(sx + sy + sz):  # Catch nans
                return np.nan, np.nan
            elif (self.H * (self.H - sx)) < (
                sy**2 + (self.req * sz / self.rpol) ** 2
            ):
                return np.nan, np.nan
            else:
                return x, y
        # Mask invisible points - we are catching nans
        # here already, so the invalid value in less is ignored
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            mask = np.where(
                np.isnan(sx + sy + sz)
                | (
                    (self.H * (self.H - sx))
                    < (sy**2 + (self.req * sz / self.rpol) ** 2)
                )
            )
        x[mask] = np.nan
        y[mask] = np.nan
        return x, y

    def _rebin(self, locs, bins):
        binned = np.digitize(locs, bins) - 1
        # These are designed to be large enough that they will throw an
        # exception if you try and use them to locate a pixel in GOES data
        errval = -999999
        # Locs is sometimes a nan. Here we make use of a comparison against a
        # nan always being false. Unfortunately there is not currently a nice
        # way to do this without turning off warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            binned[np.isfinite(locs) & (locs < min(bins))] = errval
            binned[np.isfinite(locs) & (locs > max(bins))] = errval
            binned[np.isnan(locs)] = errval
        return binned

    def _create_alt_correction(self):
        # Angles as in logbook
        y, x = np.array(np.meshgrid(self.y, self.x, indexing="ij"))
        ipc = np.array(
            np.meshgrid(np.arange(0, len(self.x)), np.arange(0, len(self.y)))
        )
        llp = self.geolocate(np.array(ipc).reshape(2, -1).transpose(), interp=False)
        lon = llp[:, 0].reshape(ipc[0].shape)
        theta = np.deg2rad(lon) - self.lam0
        phi = np.deg2rad(llp[:, 1].reshape(ipc[0].shape))

        # Correction for the satellite parallax
        # (km x-y per km height)
        self.corrdx = (np.sin(x) * np.cos(theta) + np.cos(x) * np.sin(theta)) * np.cos(
            phi
        )
        self.corrdy = (
            np.cos(x) * np.sin(y) * np.cos(theta) * np.cos(phi)
            - np.sin(x) * np.sin(y) * np.sin(theta) * np.cos(phi)
            - np.cos(y) * np.sin(phi)
        )

    def geolocate(self, coords, interp=False):
        """Returns the lon/lat of the gridbox indicies passed as coords"""
        if interp:
            if not self.xint:
                self._create_interp()
            x = self.xint(np.clip(coords[:, 0].astype("int"), 0, len(self.x) - 1))
            y = self.yint(np.clip(coords[:, 1].astype("int"), 0, len(self.y) - 1))
        else:
            x = self.x[np.clip(coords[:, 0].astype("int"), 0, len(self.x) - 1)]
            y = self.y[np.clip(coords[:, 1].astype("int"), 0, len(self.y) - 1)]
        lon, lat = self._sat_to_geodetic(x, y)
        return np.concatenate(
            (np.rad2deg(lon)[:, None], np.rad2deg(lat)[:, None]), axis=1
        )

    def locate(self, coords, alt=None, interp=False):
        """Get the satellite coordinates for N lon-lat coords, shape - [N, 2]
        Performs a simple altitude parallax adjustment, passing altitude as km

        If the input location is outside the grid region, a large negative value
        is returned (so you can't use it as a grid coordinate).

        The interp flag will interpolate the corrdinate for a more accurate
        location, returning np.nan if outside the grid. Nan is preferred as it
        propagates through any following calculations."""
        x, y = self._geodetic_to_sat(np.deg2rad(coords[:, 0]), np.deg2rad(coords[:, 1]))
        if interp:
            if not self.xint:
                self._create_interp()
            output = np.concatenate(
                (self.rxint(x)[:, None], self.ryint(y)[:, None]), axis=1
            )
        else:
            output = np.concatenate(
                (self._rebin(x, self.xl)[:, None], self._rebin(y, self.yl)[:, None]),
                axis=1,
            )
        if alt is not None:
            # Correct for altitude parallax
            # Following method laid out in notebook
            # Vector dot product
            if self.corrdx is None:
                self._create_alt_correction()
            # store locations outside of grid
            ctlpos = output.astype("float")
            ctlpos[(ctlpos < 0).sum(axis=-1) > 0] = np.nan
            ctlmask = np.isnan(ctlpos)

            ctl_ind = output.astype("int")
            ctl_ind[ctlmask] = 0

            ctlpos_corrdx = (alt / self.res) * self.corrdx[
                ctl_ind[:, 1].astype("int"), ctl_ind[:, 0].astype("int")
            ]
            ctlpos_corrdx[ctlmask[:, 0]] = np.nan
            ctlpos_corrdy = (alt / self.res) * self.corrdy[
                ctl_ind[:, 1].astype("int"), ctl_ind[:, 0].astype("int")
            ]
            ctlpos_corrdy[ctlmask[:, 1]] = np.nan
            ctlpos_corr = np.concatenate(
                (ctlpos_corrdx[:, None], ctlpos_corrdy[:, None]), axis=-1
            )
            output = ctlpos + ctlpos_corr
        return output

    def points_in_radius(self, point, radius):
        pos = self.locate(np.array([point]))[0]
        offset = int(radius / (self.res))

        lon_bounds = np.clip((pos[0] - offset, pos[0] + offset), 0, len(self.x) - 1)
        lat_bounds = np.clip((pos[1] - offset, pos[1] + offset), 0, len(self.y) - 1)

        ipc = np.array(
            np.meshgrid(
                np.arange(lon_bounds[0], lon_bounds[1]),
                np.arange(lat_bounds[0], lat_bounds[1]),
            )
        )
        llp = self.geolocate(np.array(ipc).reshape(2, -1).transpose(), interp=False)
        lon_box = llp[:, 0].reshape(ipc[0].shape)
        lat_box = llp[:, 1].reshape(ipc[0].shape)

        # lon_box = lon[lat_bounds[0]:lat_bounds[1], lon_bounds[0]:lon_bounds[1]]
        # lat_box = lat[lat_bounds[0]:lat_bounds[1], lon_bounds[0]:lon_bounds[1]]

        dist = csat2.misc.geo.haversine(point[0], point[1], lon_box, lat_box)
        dlocs = np.where(dist < radius)
        return np.array(list(zip(dlocs[0] + lat_bounds[0], dlocs[1] + lon_bounds[0])))
