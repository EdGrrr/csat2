from csat2.GOES import Granule
import numpy as np
from csat2.GOES.util import get_resolution
from csat2 import locator
from netCDF4 import Dataset


def create_advection_files(sat="G16", area="RadC", gran=None):
    # These advection changes will need to be cached, they take about 30 sec to calculate
    if gran == None:
        gran = Granule.fromtext("G16.2018005.1430.RadC")

    bands = [2, 1, 4]
    bands = list(zip(bands, map(get_resolution, bands)))
    for band, res in [bands[-1]]:
        gran.download(band)
        # Get positions of granule
        lon, lat = gran.get_lonlat(band)
        # Initially calculate for a 10km shift
        dlon = 10 / (111.111 * np.cos(np.deg2rad(lat)))
        dlat = 10 / (111.111)
        # Will need two coordinate changes, for both uwind and vwind separately.
        # Assuming that can probably add the results though
        newpos = gran.locate(
            np.array(list(zip((lon + dlon).ravel(), lat.ravel()))), interp=True
        )
        oldpos = gran.get_index(13, mode=3)

        # dxlon is the change in dx caused by a 10km change in the longitude direction
        newpos_dxlon = newpos[:, 0].reshape(lon.shape) - oldpos[0]
        newpos_dylon = newpos[:, 1].reshape(lon.shape) - oldpos[1]

        newpos = gran.locate(
            np.array(list(zip((lon).ravel(), (lat + dlat).ravel()))), interp=True
        )
        newpos_dxlat = newpos[:, 0].reshape(lon.shape) - oldpos[0]
        newpos_dylat = newpos[:, 1].reshape(lon.shape) - oldpos[1]

        filename = locator.get_folder(
            "GOES", "advection", sat=sat, area=area, res=res
        ) + "/ABI_adv.{area}.{res}km.nc".format(area=area, res=res)
        with Dataset(filename, "w", format="NETCDF4") as ncdf:
            ncdf.createDimension("x", lon.shape[1])
            ncdf.createDimension("y", lon.shape[0])

            ncdf.llcx, ncdf.llcy = gran.get_llcoord(band, mode="*")

            for varname, outdata in [
                ["newpos_dxlat", newpos_dxlat],
                ["newpos_dxlon", newpos_dxlon],
                ["newpos_dylat", newpos_dylat],
                ["newpos_dylon", newpos_dylon],
            ]:
                var = ncdf.createVariable(varname, "f", ("y", "x"), zlib=True)
                var[:] = outdata


class GOESAdvector:
    def __init__(self, gran=None, sat="G16", area="RadC", interp=True):
        """Builds a advector object for a GOES granule. Currently
        only interpolates from the 2km data"""
        if interp == False:
            raise (ValueError, "Not currently implemented")
        if gran:
            self.gran = gran
        else:
            self.gran = Granule.fromtext("G16.2018005.1430.RadC")
        filename = locator.search(
            "GOES", "advection", sat=self.gran.sat, area=self.gran.area, res=2
        )[0]
        with Dataset(filename) as ncdf:
            # Check to see if gran is in same position as cached files
            if gran:
                llcx, llcy = ncdf.llcx, ncdf.llcy
                ll_gran = gran.get_llcoord(channel=13, mode="*")
                if ll_gran[0] != llcx or ll_gran[1] != llcy:
                    raise (ValueError, "Granule doesn't match cache")
            self.data = {}
            for name in [
                "newpos_dxlat",
                "newpos_dxlon",
                "newpos_dylat",
                "newpos_dylon",
            ]:
                self.data[name] = np.ma.filled(ncdf.variables[name][:], np.nan)

    def advect_field(self, data, udist, vdist, band):
        """Advects the data field udist km in the longitude direction
        and vdist km in the latitude direction"""
        xinc = (udist / 10) * self.data["newpos_dxlon"] + (vdist / 10) * self.data[
            "newpos_dxlat"
        ]
        yinc = (udist / 10) * self.data["newpos_dylon"] + (vdist / 10) * self.data[
            "newpos_dylat"
        ]

        x_ind = np.fromfunction(lambda x, y: x, data.shape, dtype=np.int)
        y_ind = np.fromfunction(lambda x, y: y, data.shape, dtype=np.int)

        nx_ind = x_ind + xinc
        ny_ind = y_ind + yinc
        clipmask = (
            (nx_ind < 0)
            * (ny_ind < 0)
            * (nx_ind > data.shape[1])
            * (ny_ind > data.shape[0])
        )
        output = data[
            np.clip(nx_ind, 0, data.shape[1]), np.clip(ny_ind, 0, data.shape[0])
        ]
        output[clipmask] = np.nan
        return output

    def advect_pixels(self, pixels, udist, vdist):
        xinc = (udist / 10) * self.data["newpos_dxlon"][pixels[0], pixels[1]] + (
            vdist / 10
        ) * self.data["newpos_dxlat"][pixels[0], pixels[1]]
        yinc = (udist / 10) * self.data["newpos_dylon"][pixels[0], pixels[1]] + (
            vdist / 10
        ) * self.data["newpos_dylat"][pixels[0], pixels[1]]
        return pixels[0] + yinc, pixels[1] + xinc
