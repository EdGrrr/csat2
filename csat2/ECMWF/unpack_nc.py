import os
import glob
import logging
import csat2
from csat2.ECMWF.download import variable_names, _convert_cds_to_vname, _get_levelstr
from csat2.ECMWF.ECMWF import readin_ERA
import numpy as np
import xarray as xr

from netCDF4 import Dataset, num2date


log = logging.getLogger(__name__)


def import_netcdf(netcdf_file, year, level, resolution, n_times, **kwargs):
    """
    Unpacks a CDS netcdf file into csat2 file storage.

    WARNING: Only use this if strictly necessary!
    Check that the supplied parameters are accurate to avoid mislabelled files, damaging the data store.

    To activate usage, supply the keyword argument confirm=True
    """
    # TODO: Detect and split on level, year. Detect n_times - may reuqire xarray. Issue: level is not given in the netcdf file?

    if not kwargs.get("confirm", False):
        raise NotImplementedError(
            "Do not use without due caution! Read the documentation."
        )

    if resolution != "0.25grid":
        raise ValueError("Only 0.25grid resolution is supported")

    levelstr = _get_levelstr(level)
    _unpack_nc(netcdf_file, year, levelstr, resolution, n_times, lst=False)


def _unpack_nc(
    netcdf_file,
    year,
    levelstr,
    resolution,
    n_times,
    lst=True,
    lst_times=["0730", "1030", "1330", "1630"],
):
    """
    Unpack a netcdf file (downloaded or raw) into the csat2 file storage
    """
    output_folder = os.path.dirname(netcdf_file)

    # We need to split by variables!
    slice_prefix = output_folder + "/temp_"
    os.system("cdo -b 32 splitname {} {}".format(netcdf_file, slice_prefix))

    os.remove(netcdf_file)

    # These are the files we have
    variable_files = glob.glob(slice_prefix + "*")
    log.info("Created files: " + " ".join(variable_files))

    for vfile in variable_files:
        # Work out the relevant longname
        with Dataset(vfile) as ncdf:
            vname = variable_names(ncdf)
            long_name = ncdf.variables[vname].long_name.lower().replace(" ", "_")
            prod = _convert_cds_to_vname(long_name)
            if "valid_time" in ncdf.variables:
                rename_time = True
        
        if rename_time:
            with Dataset(vfile, "a") as ncdf:
                ncdf = ncdf.renameVariable("valid_time", "time")

        # The output folder may not exist yet
        var_folder = csat2.locator.get_folder(
            "ECMWF",
            "ERA5",
            year=year,
            variable=prod,
            level=levelstr,
            resolution=resolution,
            time="timed",
        )

        try:
            os.makedirs(var_folder)
        except FileExistsError:
            pass

        file_pattern = "{var_folder}/{prod}_{levelstr}_cdosplit_".format(
            var_folder=var_folder, prod=prod, levelstr=levelstr
        )
        os.system(
            "cdo splitsel,{nsteps:.0f} {vfile} {file_pattern}".format(
                nsteps=n_times, vfile=vfile, file_pattern=file_pattern
            )
        )

        # Remove the variable file
        os.remove(vfile)

        # We need to sort out the names foreach of these files
        # They are currnetly names by the day of the month
        # We want to switch that to the day of the year
        files = glob.glob(file_pattern + "*")
        files.sort()

        doys = []
        for file in files:
            with Dataset(file) as ncdf:
                dtime = num2date(
                    ncdf.variables["time"][0], ncdf.variables["time"].units
                )
                doy = dtime.timetuple().tm_yday
            index = str(doy).rjust(3, "0") + ".nc"
            new = (file_pattern + index).replace("cdosplit_", "")
            os.rename(file, new)
            doys.append(doy)

        # Do the LST convert (if required)
        if lst:
            create_lst_files(year, doys, prod, levelstr, resolution, lst_times)


def create_lst_files(
    year,
    doys,
    variable,
    levelstr,
    resolution,
    lst_times=["0730", "1030", "1330", "1630"],
):
    """Creates local solar time files for a given variable and set of doys

    Typically you want to run this after you have downloaded a whole month,
    otherwise you will end upwith missing data at the end of the day. However,
    it can be used for individual days (if you want).

    year - integer, year for data

    doys - list of doys to create LST data for. Should be a list even if you
        are only using one day of data

    variable - local name of the variable to used

    levelstr - level string, either integer for pressure level, string for
        pressure level (can end in hPa) or 'surf' for a single level variable

    resolution - only tested with '1grid'

    lst_times - either a list of floats giving the hours since midnight to
        calculate LST data for or a list of strings giving the time in hours
        minutes - e.g. '0730' as half past 7

    Output - creates the LST files in the location specified in the csat2
        machine file"""
    # Create the output folder
    var_folder = csat2.locator.get_folder(
        "ECMWF",
        "ERA5",
        year=year,
        variable=variable,
        level=levelstr,
        resolution=resolution,
        time="LST",
    )

    try:
        os.makedirs(var_folder)
    except FileExistsError:
        pass

    if isinstance(lst_times[0], str):
        lst_times_float = [int(t[:2]) + int(t[2:4]) / 60 for t in lst_times]
    else:
        lst_times_float = lst_times

    for doy in doys:
        try:
            tdata = readin_ERA(
                year,
                doy,
                variable,
                level=levelstr,
                product="ERA5",
                time="timed",
                resolution=resolution,
            )
            try:
                tdata_next = readin_ERA(
                    *csat2.misc.time.doy_step(year, doy, 1),
                    variable,
                    level=levelstr,
                    product="ERA5",
                    time="timed",
                    resolution=resolution
                )
                tdata = xr.concat((tdata, tdata_next[0]), dim="time")
            except:
                blankdata = xr.DataArray(
                    np.zeros((1, tdata.shape[1], tdata.shape[2])) * np.nan,
                    coords={
                        "time": np.array(
                            [
                                csat2.misc.time.ydh_to_datetime(
                                    *csat2.misc.time.doy_step(year, doy, 1), 0
                                )
                            ]
                        ),
                        "lon": tdata.lon,
                        "lat": tdata.lat,
                    },
                    dims=tdata.dims,
                )
                tdata = xr.concat((tdata, blankdata), dim="time")

            # Convert from ns to hours
            times = (tdata.time - tdata.time[0]).values.astype("float") / 1e9 / 3600
            longitudes = tdata.lon.values

            # DO the actual localtime calculation
            gdata = tdata.values.swapaxes(0, 2)  # .swapaxes(0, 1)
            output_tdata = [
                csat2.misc.time.toLocalSolarTime(ltime, times, longitudes, gdata)[
                    None, :, :
                ]
                for ltime in lst_times_float
            ]
            output_tdata = np.concatenate(output_tdata, axis=0)

            fill_value = -1000.0
            filename = "{}/{}_{}_{:0>3}.nc".format(var_folder, variable, levelstr, doy)
            with Dataset(filename, "w") as ncdf:
                ncdf.createDimension("lat", len(tdata.lat))
                ncdf.createDimension("lon", len(tdata.lon))
                ncdf.createDimension("time", len(lst_times))
                latitude = ncdf.createVariable("lat", "f", ("lat",))
                latitude[:] = tdata.lat
                longitude = ncdf.createVariable("lon", "f", ("lon",))
                longitude[:] = tdata.lon

                setattr(
                    ncdf,
                    "title",
                    "ERA5 "
                    + variable
                    + " data on a regular 1degree grid, retimed to LST",
                )
                setattr(
                    ncdf,
                    "institution",
                    "Space and Atmospheric Physics Group, Imperial College London",
                )
                setattr(ncdf, "user", "gryspeerdt")
                setattr(ncdf, "source", "ERA5 data from ECMWF/Copernicus")
                setattr(ncdf, "_FillValue", fill_value)

                Vtime = ncdf.createVariable("time", "f", ("time",))
                Vtime[:] = np.array(lst_times).astype("float32")

                newvar = ncdf.createVariable(variable, "f", ("time", "lat", "lon"))
                newvar[:] = output_tdata.astype("float32")
                setattr(newvar, "longname", variable.replace("_", " "))
            log.info("Created LST: {}".format(os.path.basename(filename)))
        except FileNotFoundError:
            continue
        except Exception as e:
            log.error("Readin failed on day ", doy, year)
            log.error(e)
            # print("Readin failed on day ", doy, year)
