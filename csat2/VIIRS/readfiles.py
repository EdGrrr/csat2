import numpy as np
from netCDF4 import Dataset
import logging
from csat2 import locator
import csat2.misc.time
import csat2.misc.geo

log = logging.getLogger(__name__)

DEFAULT_COLLECTION = "5200"

print("Default VIIRS collection: {}".format(DEFAULT_COLLECTION))

R_E = csat2.misc.geo.R_E

########################################################################
# Functions for reading in Satellite data, VIIRS products
########################################################################
#
#  ERI Gryspeerdt, Imperial College London, 2020
########################################################################


def readin(product, *args, **kwargs):
    """Reads in a single granule of VIIRS data from the specified product.

    Usage:
    readin( string - VIIRS product name
            int - year
            int - day of year (Jan 1st = 1)
            list - list of Datasets required)"""
    if product == "GEOMETA":
        return readin_VIIRS_GEOMETA(*args, **kwargs)
    else:
        return readin_VIIRS_L2(product, *args, **kwargs)


def dictkeysappend(dict1, dict2, axis=0):
    """Appends data in dict2 to already existing keys in a dictionary dict1.
    If that key does not exist, it is created"""
    for key in dict2.keys():
        try:
            if axis == "longest":
                axis = np.argmax(dict1[key].shape)
                dict1[key] = np.concatenate((dict1[key], dict2[key]), axis=axis)
                axis = "longest"
            else:
                dict1[key] = np.concatenate((dict1[key], dict2[key]), axis=axis)
        except KeyError:
            dict1[key] = dict2[key]
    return dict1


def readin_VIIRS_L2(
    product, year, doy, sds=None, time=None, col=None, varind=None, scale=True
):
    """Reads in single day of VIIRS L2 data"""
    if col is None:
        col = DEFAULT_COLLECTION
    if sds is None:
        sds = []
    if isinstance(sds, str):
        sds = [sds]
    if isinstance(time, type("")):
        filename = locator.search(
            "VIIRS", product, year=year, doy=doy, collection=col, time=time
        )[0]
    else:
        filename = locator.search(
            "VIIRS", product, year=year, doy=doy, collection=col, time=time[0]
        )[0]
    log.info(filename)
    return readin_VIIRS_L2_filename(filename, sds, scale)


def readin_metadata(product, year, doy, time, col, sds):
    if isinstance(time, int):
        time = "{:0>4}".format(
            time
        )  # time should be a string to cope with leading zeros
    filename = locator.search(
        "VIIRS", product, year=year, doy=doy, collection=col, time=time
    )[0]
    with Dataset(filename) as ncdf:
        group = netcdf_search_groups(ncdf, sds)
        return {
            a: ncdf.groups[group].variables[sds].getncattr(a)
            for a in ncdf.groups[group].variables[sds].ncattrs()
        }


def netcdf_search_groups(ncdf, varname):
    """Determines which group a requested variable is in
    Args:
        ncdf - netcdf Dataset object
        varname - str - variable name
    Returns
        str - group name"""
    groups = ncdf.groups.keys()
    log.debug(groups)
    for group in groups:
        if varname in ncdf.groups[group].variables.keys():
            return group


def readin_VIIRS_L2_filename(filename, names, scale=True):
    """Reads the requested variable names from a specified file into a dictionary

    Args:
        filename - str - filename
        names - list[str] - list of variable names
        scale - use netcdf4 autoscale. Set to False to read integers
    Returns:
        dict - output data arrays, key is variable name"""
    with Dataset(filename) as ncdf:
        ncdf.set_auto_scale(False)
        indata = {}
        for name in names:
            group = netcdf_search_groups(ncdf, name)
            log.debug(group, name)
            indata[name] = ncdf.groups[group].variables[name][:]
            indata[name] = indata[name].astype("float")
            if scale:
                indata[name][indata[name] > 65530] = np.nan
                try:
                    indata[name] *= ncdf.groups[group].variables[name].scale_factor
                except AttributeError:
                    pass
                try:
                    indata[name] += ncdf.groups[group].variables[name].add_offset
                except AttributeError:
                    pass
    return indata


def readin_VIIRS_GEOMETA(year, doy, sat="N"):
    if sat not in ["N", "J"]:
        raise IOError("Not a suitable satellite (N, J)")
    satname = {"N": "NPP", "J": "JPSS1"}[sat]
    _, mon, day = csat2.misc.time.doy_to_date(year, doy)
    filenames = locator.search(
        "VIIRS", "GEOMETA", year=year, mon=mon, day=day, sat=satname
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


def list_available_times(product, year, doy, col=None):
    """Returns the available times for a specific product and day"""
    if col is None:
        col = DEFAULT_COLLECTION
    filenames = sorted(
        locator.search("VIIRS", product, year=year, doy=doy, collection=col, time="*")
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
    metadata = readin_VIIRS_GEOMETA(year, doy, sat)
    if day_night is None:
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
