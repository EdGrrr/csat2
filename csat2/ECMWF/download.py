import cdsapi
import csat2
from csat2.ECMWF.variables import _convert_vname_to_cds, _convert_cds_to_vname, _get_levelstr
from csat2.ECMWF.ECMWF import variable_names, readin_ERA
from csat2.ECMWF.unpack_nc import _unpack_nc
import os
import pkg_resources
import glob
import logging
from netCDF4 import Dataset, num2date
import xarray as xr
import numpy as np

log = logging.getLogger(__name__)

def check(
    year,
    month,
    variables,
    level,
    resolution,
    days=None,
    doys=None,
    time="timed",
    product="ERA5",
):
    """Check that files exist for the year, month, level and variables in question
    Does not do any validation of the contents at the moment.

    Returns a tuple:
        boolean - all files are present
        list - missing [variable, doy] pairs"""

    if days:
        doys = [csat2.misc.time.date_to_doy(year, month, day)[1] for day in days]
    elif doys:
        doys = doys
    else:
        if month < 12:
            doys = np.arange(
                csat2.misc.time.date_to_doy(year, month, 1)[1],
                csat2.misc.time.date_to_doy(year, month + 1, 1)[1],
            )
        else:
            doys = np.arange(
                csat2.misc.time.date_to_doy(year, 12, 1)[1],
                csat2.misc.time.date_to_doy(year, 12, 31)[1] + 1,
            )

    missing = []
    exist = True
    for variable in variables:
        for doy in doys:
            nfiles = len(
                csat2.locator.search(
                    "ECMWF",
                    product,
                    year=year,
                    doy=doy,
                    variable=variable,
                    resolution=resolution,
                    time=time,
                    level=level,
                )
            )
            if nfiles == 0:
                missing.append([variable, doy])
                exist = False
    return exist, missing


def download(
    year,
    month,
    variables,
    level,
    resolution,
    days=None,
    times=8,
    lst=True,
    lst_times=["0730", "1030", "1330", "1630"],
    force_redownload=False,
):
    """Downloads ECMWF ERA5 data for a single month. Multiple variables can be selected,
    but only a single level at each time.

    days - download specific days, leave a None to get the whole month.
        You typically only need this for NRT data.

    level - the pressure level of the variable in hPa (or 'surf' for
        surface data)

    resolution - Stored resolution of the downloaded data. Currently only
        '1grid' and '0.25grid' (native) are supported.

    times - Either a list of the times to be downloaded, or an integer that
        specifies the number of times to be used per day

    lst - If True, the data is transferred to a LST grid for time times
        specified in lst_times. In this case, times, should be specified
        as an integer.

    force_redownload - forces the redownload of the data, even if it already
        exists. Use this if you are changing the internal properties of the
        files (such as the number of times), as it depends on the check function,
        which doesn't do any validation of the file contents"""

    if resolution not in ["1grid", "0.25grid"]:
        raise ValueError("Resolution: {} not yet implmented".format(resolution))

    if (
        not force_redownload
        and check(year, month, variables, level, resolution, days=days)[0]
    ):
        log.info("All files exist")
        return

    if isinstance(variables, str):
        variables = [variables]
    cds_variables = [_convert_vname_to_cds(vname) for vname in variables]

    if not days:
        days = range(1, 32)
    days = ["{:0>2}".format(day) for day in days]

    if isinstance(times, str):
        times = [times]
    elif isinstance(times, int):
        interval = 24 // times
        times = ["{:0>2}:00".format(t) for t in range(0, 24, interval)]

    levelstr = _get_levelstr(level)

    # Select a download location - file will be removed
    # anyway, so not a huge issue
    output_folder = csat2.locator.get_folder(
        "ECMWF",
        "ERA5",
        year=year,
        variable=variables[0],  # Variables is a list
        level=levelstr,
        resolution=resolution,
        time="timed",
    )

    try:
        os.makedirs(output_folder)
    except FileExistsError:
        pass

    download_file = output_folder + "/download.nc"

    # Currently setup to block waiting for the download
    # In future it may be better just to submit the request,
    # but that can wait
    c = cdsapi.Client()
    log.info("Submitting CDS request")
    if level == "surf":
        c.retrieve(
            "reanalysis-era5-single-levels",
            {
                "variable": cds_variables,
                "product_type": "reanalysis",
                "year": str(year),
                "month": ["{:0>2}".format(month)],
                "day": days,
                "time": times,
                "format": "netcdf",
            },
            download_file,
        )
    else:
        c.retrieve(
            "reanalysis-era5-pressure-levels",
            {
                "variable": cds_variables,
                "pressure_level": level,
                "product_type": "reanalysis",
                "year": str(year),
                "month": ["{:0>2}".format(month)],
                "day": days,
                "time": times,
                "format": "netcdf",
            },
            download_file,
        )

    # Regrid the downloaded file (if required)
    if resolution == "1grid":
        gridfile = pkg_resources.resource_filename("csat2", "data/d1_grid")
        os.system("cdo remapbil,{0} {1} {1}.regrid".format(gridfile, download_file))
        os.system("mv {0}.regrid {0}".format(download_file))
    elif resolution == "0.25grid":
        # This should be the native download resolution, but best to check.
        with Dataset(download_file) as ncdf:
            vname = variable_names(ncdf)
            vshape = ncdf.variables[vname].shape
            if (vshape[-1] != 1440) or (vshape[-2] != 721):
                raise ValueError("Invalid data size for 0.25grid - {}".format(vshape))

    _unpack_nc(download_file, year, levelstr, resolution, len(times), lst=lst, lst_times=lst_times)


