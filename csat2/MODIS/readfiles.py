from __future__ import print_function, division
from csat2 import locator
from csat2.misc.time import doy_to_date, ydh_to_datetime
import numpy as np
from netCDF4 import Dataset
from csat2.misc import hdf
import logging
import xarray as xr

log = logging.getLogger(__name__)

####################
# Collection setup #
####################

DEFAULT_COLLECTION = "61"
log.info("MODIS default is collection {}".format(DEFAULT_COLLECTION))


def setCollection(col):
    global DEFAULT_COLLECTION
    DEFAULT_COLLECTION = col
    log.info("MODIS default is collection {}".format(DEFAULT_COLLECTION))


###################
# Readin function #
###################


def readin(product, *args, **kwargs):
    """Reads in a single day of MODIS data from the specified product.
    Currently this will work with M*D08_D3, M*DATML2, M*D04_L2, M*D06_L2
    the GEOMETA data and and Ed's subset of MODIS data ('subset')

    Usage:
    readin( string - MODIS product name
            int - year
            int - day of year (Jan 1st = 1)
            list - list of Datasets required, Lat, Lon and Scan time
                                are automatic)"""
    if product == "MOD08_D3":
        return readin_MODIS_L3(*args, sat="terra", **kwargs)
    elif product == "MYD08_D3":
        return readin_MODIS_L3(*args, sat="aqua", **kwargs)
    elif product in [
        "MODATML2",
        "MYDATML2",
        "MOD021KM",
        "MYD021KM",
        "MOD02HKM",
        "MYD02HKM",
        "MOD02QKM",
        "MYD02QKM",
        "MOD03",
        "MYD03",
        "MOD04_L2",
        "MYD04_L2",
        "MOD06_L2",
        "MYD06_L2",
        "MODATML2_NRT",
        "MYDATML2_NRT",
        "MOD021KM_NRT",
        "MYD021KM_NRT",
        "MOD02HKM_NRT",
        "MYD02HKM_NRT",
        "MOD02QKM_NRT",
        "MYD02QKM_NRT",
        "MOD03_NRT",
        "MYD03_NRT",
        "MOD04_L2_NRT",
        "MYD04_L2_NRT",
        "MOD06_L2_NRT",
        "MYD06_L2_NRT",
    ]:
        return readin_MODIS_L2(product, *args, **kwargs)
    elif product == "subset":
        return readin_MODIS_subset(*args, **kwargs)
    elif product == "cdnc_best":
        return readin_MODIS_cdnc_best(*args, **kwargs)
    elif product == "GEOMETA":
        return readin_MODIS_GEOMETA(*args, **kwargs)
    else:
        print("Product not implemented yet")


def readin_MODIS_L3(year, doy, sds, sat="aqua", col=DEFAULT_COLLECTION):
    """Reads in MODIS L3 data from the collection 5 daily files on cloud"""
    product = {"aqua": "MYD08_D3", "terra": "MOD08_D3"}[sat]
    filename = locator.search("MODIS", product, year=year, doy=doy, collection=col)[0]
    log.info(filename)
    return readin_MODIS_L3_filename(filename, sds)


def readin_MODIS_L3_filename_xarray(filename, names):
    """Note this is currently three times as slow as the old csat implmentation.
    Consider re-writing to base more on that"""
    # Keep the attributes, even through scaling!
    with xr.set_options(keep_attrs=True):
        # Make sure to turn off mask and scale. MODIS does this backwards!
        indata = xr.open_dataset(filename, decode_cf=False, mask_and_scale=False)
        outdata = indata[names]
        for name in names:
            outdata[name] = outdata[name].where(
                outdata[name] != outdata[name]._FillValue
            )
            try:
                outdata[name] = (outdata[name] - outdata[name].add_offset) * outdata[
                    name
                ].scale_factor
            except AttributeError:
                pass
        dim_names = {}
        for dim in outdata.dims:
            ndim = dim.replace(":mod08", "")
            outdata[ndim] = indata[ndim]
            dim_names[dim] = ndim
        outdata = (
            outdata.rename_dims(dim_names)
            .set_coords(dim_names.values())
            .rename({"XDim": "lon", "YDim": "lat"})
        )
    return outdata


def readin_MODIS_L3_filename(filename, names):
    """Builds up the xarray dataset, approximately twice as fast
    Misses out some attributes, but could be added in as required"""
    with hdf.Dataset(filename) as ncdf:
        ds = xr.Dataset()  # New output dataset
        tdims = []
        for name in names:
            var = ncdf.variables[name]
            var.set_auto_scale(False)
            dims = ["time"] + [a.replace(":mod08", "") for a in var.dimensions]
            units = var.units
            try:
                indata = (var[:] - var.add_offset) * var.scale_factor
            except AttributeError:  # No scale factors
                indata = var[:]
            indata = np.ma.filled(indata, np.nan)
            ds[name] = xr.DataArray(indata[None, :, :], dims=dims)
            ds[name].attrs["units"] = units
            tdims.extend(dims[1:])
        tdims = set(tdims)
        for tdim in tdims:
            ds[tdim] = xr.DataArray(ncdf.variables[tdim][:], dims=(tdim,))
        return ds


def readin_MODIS_subset(year, doy, sds=None, sat="aqua", col=DEFAULT_COLLECTION):
    sat = {"aqua": "", "terra": ".terra"}[sat]
    datafilename = locator.search(
        "MODIS", "subset", year=year, doy=doy, collection=col, sat=("new" + sat)
    )[0]
    indata = {}
    with Dataset(datafilename, "r") as ncdf:
        ds = xr.Dataset()  # New output dataset
        tdims = []
        for name in sds:
            var = ncdf.variables[name]
            var.set_auto_scale(False)
            dims = var.dimensions
            indata = np.ma.filled(var[:], np.nan)
            ds[name] = xr.DataArray(indata, dims=dims)
            tdims.extend(dims)
        tdims = set(tdims)
        for tdim in tdims:
            if tdim == "time":
                ds["time"] = xr.DataArray(
                    np.array([ydh_to_datetime(year, doy, 0)]), dims=("time",)
                )
            else:
                ds[tdim] = xr.DataArray(ncdf.variables[tdim][:], dims=(tdim,))
        return ds


def readin_MODIS_cdnc_best(year, doy, sds=None, version="1", sat="aqua", resolution=1):
    filename = locator.search(
        "MODIS",
        "cdnc_best",
        year=year,
        doy=doy,
        version=version,
        sat={"aqua": "A", "terra": "T"}[sat],
        res={1: "", 0.25: "q"}[resolution],
    )[0]
    with Dataset(filename, "r") as ncdf:
        ds = xr.Dataset()  # New output dataset
        tdims = []

        for name in sds:
            var = ncdf.variables[name]
            dims = var.dimensions
            try:
                # fill = getattr(var, '_FillValue')
                indata = var[:].filled(np.nan)
            except ValueError:
                indata = var[:]
            except AttributeError:
                # Trying to put a nan in an interger array are we?
                indata = var[:]
            ds[name] = xr.DataArray(indata, dims=dims)
            tdims.extend(dims)
        tdims = set(tdims)
        for tdim in tdims:
            if tdim == "time":
                ds["time"] = xr.DataArray(
                    np.array([ydh_to_datetime(year, doy, 0)]), dims=("time",)
                )
            elif tdim == "lat":
                ds[tdim] = xr.DataArray(ncdf.variables[tdim][:][::-1], dims=(tdim,))
            else:
                ds[tdim] = xr.DataArray(ncdf.variables[tdim][:], dims=(tdim,))
    return ds


def readin_MODIS_L2(
    product, year, doy, time, sds=None, col=DEFAULT_COLLECTION, varind=None
):
    """Reads in single day of MODIS L2 data"""
    if sds == None:
        sds = []
    sds += ["Latitude", "Longitude"]

    filename = locator.search(
        "MODIS", product, year=year, doy=doy, collection=col, time=time
    )[0]
    log.debug(filename)

    if len(filename) == 0:
        raise IOError("No files for " + str(year) + " " + str(doy))

    if varind is not None:
        return readin_MODIS_L2_filename_fast(filename, sds, varind)
    else:
        return readin_MODIS_L2_filename(filename, sds)


def readin_metadata(product, year, doy, time, col, sds):
    if type(time) == type(1):
        time = "{:0>4}".format(time)
    filename = locator.search(
        "MODIS", product, year=year, doy=doy, collection=col, time=time
    )[0]
    if isinstance(sds, type("")):
        log.debug("Listifying")
        sds = [sds]
    metadata = {}
    with hdf.Dataset(filename) as ncdf:
        for varname in sds:
            metadata[varname] = {}
            for metname in ncdf.variables[varname].ncattrs():
                metadata[varname][metname] = ncdf.variables[varname].getncattr(metname)
    return metadata


def remove_dim_suffix(dim):
    if ":" in dim:
        return ":".join(dim.split(":")[:-1])
    else:
        return dim


def readin_MODIS_L2_filename(filename, names):
    log.debug(filename)
    with hdf.Dataset(filename) as ncdf:
        ds = xr.Dataset()
        tdims = []
        for name in names:
            var = ncdf.variables[name]
            var.set_auto_scale(False)
            var.set_auto_mask(False)
            dims = [
                a.replace(":mod08", "").replace(":MODIS_SWATH_Type_L1B", "")
                for a in var.dimensions
            ]
            try:
                indata = (var[:] - var.add_offset) * var.scale_factor
            except AttributeError:  # No scale factors
                indata = var[:]
            try:
                indata = np.where(indata == var._Fillvalue, np.nan, indata)
            except AttributeError:  # No fill value
                pass
            ds[name] = xr.DataArray(indata, dims=dims)
            try:
                ds[name].attrs["units"] = var.units
            except AttributeError:
                pass
            tdims.extend(dims)
        tdims = set(tdims)
        for tdim in tdims:
            try:
                ds[tdim] = xr.DataArray(ncdf.variables[tdim][:], dims=(tdim,))
            except KeyError:
                pass
    return ds


def readin_MODIS_L2_filename_fast(filename, names, varind=None):
    # If the input is a string, put it into a list
    log.debug(filename)
    with hdf.Dataset(filename) as ncdf:
        ds = xr.Dataset()
        tdims = []
        for name in names:
            var = ncdf.variables[name]
            var.set_auto_scale(False)
            var.set_auto_mask(False)
            if (varind is not None) and len(ncdf.variables[name].shape) > 2:
                vdata = var[varind]
                dims = [
                    a.replace(":mod08", "").replace(":MODIS_SWATH_Type_L1B", "")
                    for a in var.dimensions[1:]
                ]
            else:
                vdata = var[:]
                dims = [
                    a.replace(":mod08", "").replace(":MODIS_SWATH_Type_L1B", "")
                    for a in var.dimensions
                ]
            try:
                indata = (vdata - var.add_offset) * var.scale_factor
            except AttributeError:  # No scale factors
                indata = vdata
            try:
                indata = np.where(indata == var._Fillvalue, np.nan, indata)
            except AttributeError:  # No fill value
                pass
            ds[name] = xr.DataArray(indata, dims=dims)
            try:
                ds[name].attrs["units"] = var.units
            except AttributeError:
                pass
            tdims.extend(dims)
        tdims = set(tdims)
        for tdim in tdims:
            try:
                ds[tdim] = xr.DataArray(ncdf.variables[tdim][:], dims=(tdim,))
            except KeyError:
                pass
        return ds


def readin_MODIS_GEOMETA(year, doy, sat="aqua"):
    if sat not in ["aqua", "terra"]:
        raise IOError("Not a suitable satellite (choose 'aqua' or 'terra')")
    _, mon, day = doy_to_date(year, doy)
    filenames = locator.search(
        "MODIS", "GEOMETA", year=year, mon=mon, day=day, sat=sat.upper()
    )
    if len(filenames) == 0:
        raise IOError("No files for " + str(year) + " " + str(doy))
    return np.genfromtxt(
        filenames[0],
        skip_header=3,
        delimiter=",",
        dtype=[
            ("GranID", "U50"),
            ("StartTime", "U20"),
            ("Archive", "i8"),
            ("Orbit", "i8"),
            ("DayNight", "U1"),
            ("EastBC", "f8"),
            ("NorthBC", "f8"),
            ("SouthBC", "f8"),
            ("WestBC", "f8"),
            ("GRLon1", "f8"),
            ("GRLon2", "f8"),
            ("GRLon3", "f8"),
            ("GRLon4", "f8"),
            ("GRLat1", "f8"),
            ("GRLat2", "f8"),
            ("GRLat3", "f8"),
            ("GRLat4", "f8"),
        ],
    )


def list_available_times(product, year, doy, col=DEFAULT_COLLECTION):
    """Returns the available times for a specific product and day"""
    filenames = sorted(
        locator.search("MODIS", product, year=year, doy=doy, collection=col, time="*")
    )
    return list(map(lambda x: x.split("/")[-1].split(".")[2], filenames))


def intersect_times(
    year,
    doy,
    sat,
    bbox,
    day_night=None,
):
    """Returns the times where a MODIS granule intersects a bounding box (bbox), defined as
    bbox = (lat_max,lat_min,lon_max,lon_min)
    lat in range[-90,90]
    lon in range[-180,180]
    sat in ['aqua','terra']
    day_night - string combination of 'D','N','B' - e.g. 'DB'"""
    metadata = readin_MODIS_GEOMETA(year, doy, sat)
    if day_night == None:
        times = [
            str(entry[0]).split(".")[2]
            for entry in metadata
            if _intersects(
                (entry["NorthBC"], entry["SouthBC"], entry["EastBC"], entry["WestBC"]),
                bbox,
            )
        ]
    else:
        times = [
            str(entry[0]).split(".")[2]
            for entry in metadata
            if _intersects(
                (entry["NorthBC"], entry["SouthBC"], entry["EastBC"], entry["WestBC"]),
                bbox,
            )
            & (entry["DayNight"] in day_night)
        ]
    return times


def _intersects(entry, bbox):
    """Returns true if the entry (latmax,latmin,lonmax,lonmin)
    intersects the bbox"""
    return (
        (entry[0] > bbox[1])
        & (entry[1] < bbox[0])
        & (entry[2] > bbox[3])
        & (entry[3] < bbox[2])
    )


def geoloc_interpolate(geoloc):
    """Interpolates geolocation data, intended for use
    only with MODIS_L2_regrid"""
    geolocexpand = np.zeros(5 * len(geoloc))
    dtln_flag = (abs(max(geoloc)) + abs(min(geoloc))) / 2
    for x in range(0, len(geoloc) - 1):
        for y in range(0, 5):
            geolocexpand[5 * x + y] = geoloc[x] + (geoloc[x + 1] - geoloc[x]) / 5 * y
            if (geoloc[x + 1] - geoloc[x]) > dtln_flag:
                geolocexpand[5 * x + y] = geoloc[x] + (
                    geoloc[x + 1] - geoloc[x] - 360
                ) / (5 * y)
    x = x + 1
    for y in range(0, 5):
        geolocexpand[5 * x + y] = geoloc[x] + (geoloc[x] - geoloc[x - 1]) / 5 * y
    return geolocexpand


def field_interpolate(data, factor=5, dateline_check=False):
    """Interpolates the given 2D field by the factor, edge pixels are defined
    by the ones in the centre, odd factords only!

    dateline_check will now catch cases where the granule crosses the dateline
    and handle this.  Only intended for MODIS granules though! Still not
    suitable for granules that cross the pole.
    """
    dateline_flag = False

    if dateline_check:
        # Only three signs need to be checked at the dateline must
        # enter and leave
        dmax, dmin = data.max(), data.min()
        if (dmax > 175) and (dmin < -175):
            """-180 is minimum longitude"""
            data[data < 0] += 360
            dateline_flag_add = True
            # print('Dateline')
            dateline_flag = True
        elif (dmax > 355) and (dmin < 5):
            """0 is minimum longitude"""
            data[data > 180] -= 360
            dateline_flag_add = False
            # print('Dateline')
            dateline_flag = True

    output = np.zeros((factor * data.shape[0], factor * data.shape[1])) * np.nan
    output[int(factor // 2) :: factor, int(factor // 2) :: factor] = data
    for i in range(1, factor + 1):
        output[(int(factor // 2) + i) : (-1 * factor // 2 + 1) : factor, :] = (
            i
            * (
                (
                    output[int(factor // 2) + factor :: factor, :]
                    - output[int(factor // 2) : (-1 * factor) : factor, :]
                )
                / float(factor)
            )
            + output[int(factor // 2) : (-1 * factor) : factor, :]
        )

    prediff = output[(int(factor // 2) + 1), :] - output[(int(factor // 2)), :]
    for i in range(1, int(factor // 2) + 1):
        output[(int(factor // 2) - i), :] = output[int(factor // 2), :] - prediff * i

    postdiff = output[-1 * (factor // 2 + 1), :] - output[-1 * (factor // 2 + 2), :]
    for i in range(1, int(factor / 2) + 1):
        output[-1 * (factor // 2 + 1) + i, :] = (
            output[-1 * (factor // 2 + 1), :] + postdiff * i
        )

    for i in range(1, factor + 1):
        output[:, (int(factor // 2) + i) : (-1 * factor // 2 + 1) : factor] = (
            i
            * (
                (
                    output[:, int(factor // 2) + factor :: factor]
                    - output[:, int(factor // 2) : (-1 * factor) : factor]
                )
                / float(factor)
            )
            + output[:, int(factor // 2) : (-1 * factor) : factor]
        )

    prediff = output[:, (int(factor // 2) + 1)] - output[:, (int(factor // 2))]
    for i in range(1, int(factor // 2) + 1):
        output[:, (int(factor // 2) - i)] = output[:, int(factor // 2)] - prediff * i

    postdiff = output[:, -1 * (factor // 2 + 1)] - output[:, -1 * (factor // 2 + 2)]
    for i in range(1, int(factor // 2) + 1):
        output[:, -1 * (factor // 2 + 1) + i] = (
            output[:, -1 * (factor // 2 + 1)] + postdiff * i
        )

    if dateline_flag:
        if dateline_flag_add is True:
            output[output > 180] -= 360
        elif dateline_flag_add is False:
            output[output < 0] += 360
    return output
