import cdsapi
import csat2
from csat2.ECMWF.ECMWF import variable_names, readin_ERA
import os
import pkg_resources
import glob
import logging
from netCDF4 import Dataset, num2date
import xarray as xr
import numpy as np


# This translates between local names and the cdsapi names
# Partly to ensure independence if they decide to change names
# in the future
# The list is ['local_name', 'cds_name']
CDS_NAME_TABLE = [['U-wind-component', 'u_component_of_wind'],
                  ['V-wind-component', 'v_component_of_wind'],
                  ['SST', 'sea_surface_temperature']]
LOCAL_NAMES, CDS_NAMES = list(zip(*CDS_NAME_TABLE))


def _convert_vname_to_cds(vname):
    try:
        return CDS_NAMES[LOCAL_NAMES.index(vname)]
    except:
        return vname.lower()


def _convert_cds_to_vname(cds_name):
    try:
        return LOCAL_NAMES[CDS_NAMES.index(cds_name)]
    except:
        return cds_name.capitalize()


def check(year, month, variables, level, resolution, days=None, doys=None, time='timed', product='ERA5'):
    '''Check that files exist for the year, month, level and variables in question
    Does not do any validation of the contents at the moment.

    Returns a tuple:
        boolean - all files are present
        list - missing [variable, doy] pairs'''

    if days:
        doys = [csat2.misc.time.date_to_doy(year, month, day)[1] for day in days]
    elif doys:
        doys = doys
    else:
        if month < 12:
            doys = np.arange(csat2.misc.time.date_to_doy(year, month, 1)[1],
                             csat2.misc.time.date_to_doy(year, month+1, 1)[1])
        else:
            doys = np.arange(csat2.misc.time.date_to_doy(year, 12, 1)[1],
                             csat2.misc.time.date_to_doy(year, 12, 31)[1]+1)

    missing = []
    exist = True
    for variable in variables:
        for doy in doys:
            nfiles = len(csat2.locator.search(
                'ECMWF', product, year=year, doy=doy,
                variable=variable, resolution=resolution,
                time=time, level=level))
            if nfiles == 0:
                missing.append([variable, doy])
                exist = False
    return exist, missing


def download(year, month, variables, level, resolution,
             days=None,
             times=8,
             lst=True, lst_times=['0730', '1030', '1330', '1630'],
             force_redownload=False):
    '''Downloads ECMWF data for a single month. Multiple variables can be selected,
    but only a single level at each time.

    days - download specific days, leave a None to get the whole month.
        You typically only need this for NRT data.

    level - the pressure level of the variable in hPa (or 'surf' for 
        surface data)

    resolution - Stored resolution of the downloaded data. Currently only 
        '1grid' is supported.

    times - Either a list of the times to be downloaded, or an integer that 
        specifies the number of times to be used per day

    lst - If True, the data is transferred to a LST grid for time times 
        specified in lst_times. In this case, times, should be specified
        as an integer.

    force_redownload - forces the redownload of the data, even if it already 
        exists. Use this if you are changing the internal properties of the
        files (such as the number of times), as it depends on the check function,
        which doesn't do any validation of the file contents'''

    if resolution is not '1grid':
        raise ValueError('Resolution: {} not yet implmented'.format(resolution))

    if (not force_redownload and check(
            year, month, variables, level, resolution, days=days)[0]):
        logging.info('All files exist')
        return
    
    if isinstance(variables, str):
        variables = [variables]
    cds_variables = [_convert_vname_to_cds(vname)
                     for vname in variables]

    if not days:
        days = range(1, 32)
    days = ['{:0>2}'.format(day) for day in days]

    if isinstance(times, str):
        times = [times]
    elif isinstance(times, int):
        interval = 24//times
        times = ['{:0>2}:00'.format(t)
                 for t in range(0, 24, interval)]

    if level == 'surf':
        levelstr = 'surf'
    else: # level is some kind of pressure level
        if isinstance(level, str):
            level = level.replace('hPa', '')
        if isinstance(level, int):
            level = '{.0f}'.format(level)
        levelstr = '{}hPa'.format(level)

    # Select a download location - file will be removed
    # anyway, so not a huge issue
    output_folder = csat2.locator.get_folder(
        'ECMWF', 'ERA5',
        year=year,
        variable=variables[0], # Variables is a list
        level=levelstr,
        resolution=resolution,
        time='timed')

    try:
        os.makedirs(output_folder)
    except FileExistsError:
        pass

    download_file = output_folder + '/download.nc'

    # Currently setup to block waiting for the download
    # In future it may be better just to submit the request,
    # but that can wait
    c = cdsapi.Client()
    logging.info('Submitting CDS request')
    if level == 'surf':
        c.retrieve(
            'reanalysis-era5-single-levels',
            {
                'variable': cds_variables,
                'product_type': 'reanalysis',
                'year': str(year),
                'month': ['{:0>2}'.format(month)],
                'day': days,
                'time': times,
                'format': 'netcdf',
            },
            download_file)
    else:
        c.retrieve(
            'reanalysis-era5-pressure-levels',
            {
                'variable': cds_variables,
                'pressure_level': level,
                'product_type': 'reanalysis',
                'year': str(year),
                'month': ['{:0>2}'.format(month)],
                'day': days,
                'time': times,
                'format': 'netcdf',
            },
            download_file)

    # Regrid the downloaded file (if required)
    if resolution == '1grid':
        gridfile = pkg_resources.resource_filename(
            'csat2', 'data/d1_grid')
        os.system('cdo remapbil,{0} {1} {1}.regrid'.format(
            gridfile,
            download_file))
        os.system('mv {0}.regrid {0}'.format(download_file))

    # We need to split by variables!
    slice_prefix = output_folder + '/temp_'
    os.system('cdo -b 32 splitname {} {}'.format(
        download_file,
        slice_prefix))

    os.remove(download_file)

    # These are the files we have
    variable_files = glob.glob(slice_prefix+'*')
    logging.info('Created files: '+' '.join(variable_files))

    for vfile in variable_files:
        # Work out the relevant longname
        with Dataset(vfile) as ncdf:
            vname = variable_names(ncdf)
            long_name = ncdf.variables[vname].long_name.lower().replace(' ', '_')
            prod = _convert_cds_to_vname(long_name)

        # The output folder may not exist yet
        var_folder = csat2.locator.get_folder(
            'ECMWF', 'ERA5',
            year=year,
            variable=prod,
            level=levelstr,
            resolution=resolution,
            time='timed')

        try:
            os.makedirs(var_folder)
        except FileExistsError:
            pass

        file_pattern = '{var_folder}/{prod}_{levelstr}_cdosplit_'.format(
            var_folder=var_folder,
            prod=prod,
            levelstr=levelstr)
        os.system(
            'cdo splitsel,{nsteps:.0f} {vfile} {file_pattern}'.format(
                nsteps=24//interval, vfile=vfile,
                file_pattern=file_pattern))

        # Remove the variable file
        os.remove(vfile)

        # We need to sort out the names foreach of these files
        # They are currnetly names by the day of the month
        # We want to switch that to the day of the year
        files = glob.glob(file_pattern+'*')
        files.sort()

        doys = []
        for file in files:
            with Dataset(file) as ncdf:
                dtime = num2date(
                    ncdf.variables['time'][0],
                    ncdf.variables['time'].units)
                doy = dtime.timetuple().tm_yday
            index = str(doy).rjust(3, '0')+'.nc'
            new = (file_pattern+index).replace('cdosplit_', '')
            os.rename(file, new)
            doys.append(doy)

        # Do the LST convert (if required)
        if lst:
            create_lst_files(year, doys, prod, levelstr,
                             resolution, lst_times)


def create_lst_files(year, doys, variable, levelstr, resolution,
                     lst_times=['0730', '1030', '1330', '1630']):
    '''Creates local solar time files for a given variable and set of doys
    
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
        machine file
'''
    # Create the output folder
    var_folder = csat2.locator.get_folder(
        'ECMWF', 'ERA5',
        year=year,
        variable=variable,
        level=levelstr,
        resolution=resolution,
        time='LST')
    
    try:
        os.makedirs(var_folder)
    except FileExistsError:
        pass

    if isinstance(lst_times[0], str):
        lst_times_float = [int(t[:2])+int(t[2:4])/60 for t in lst_times]
    else:
        lst_times_float = lst_times
        
    for doy in doys:
        try:
            tdata = readin_ERA(year, doy,
                               variable, level=levelstr,
                               product='ERA5', time='timed',
                               resolution=resolution)
            try:
                tdata_next = readin_ERA(*csat2.misc.time.doy_step(year, doy, 1),
                                        variable, level=levelstr,
                                        product='ERA5', time='timed',
                                        resolution=resolution)
                tdata = xr.concat((tdata, tdata_next[0]), dim='time')
            except:
                blankdata = xr.DataArray(
                    np.zeros((1, tdata.shape[1], tdata.shape[2]))*np.nan,
                    coords={'time': np.array([
                        csat2.misc.time.ydh_to_datetime(
                            *csat2.misc.time.doy_step(year, doy, 1), 0)]),
                            'lon': tdata.lon,
                            'lat': tdata.lat},
                    dims=tdata.dims)
                tdata = xr.concat((tdata, blankdata), dim='time')

            # Convert from ns to hours
            times = (tdata.time-tdata.time[0]).values.astype('float')/1e9/3600
            longitudes = tdata.lon.values

            # DO the actual localtime calculation
            gdata = tdata.values.swapaxes(0, 2)#.swapaxes(0, 1)
            output_tdata = [csat2.misc.time.toLocalSolarTime(
                ltime, times,
                longitudes, gdata)[None, :, :]
                            for ltime in lst_times_float]
            output_tdata = np.concatenate(output_tdata, axis=0)
            
            fill_value = -1000.0
            filename = '{}/{}_{}_{:0>3}.nc'.format(var_folder, variable, levelstr, doy)
            with Dataset(filename,'w') as ncdf:
                ncdf.createDimension('lat', len(tdata.lat))
                ncdf.createDimension('lon', len(tdata.lon))
                ncdf.createDimension('time', len(lst_times))
                latitude = ncdf.createVariable('lat', 'f', ('lat',))
                latitude[:] = tdata.lat
                longitude = ncdf.createVariable('lon', 'f', ('lon',))
                longitude[:] = tdata.lon

                setattr(ncdf, 'title', 'ERA5 '+variable+' data on a regular 1degree grid, retimed to LST')
                setattr(ncdf, 'institution', 'Space and Atmospheric Physics Group, Imperial College London')
                setattr(ncdf, 'user', 'gryspeerdt')
                setattr(ncdf, 'source', 'ERA5 data from ECMWF/Copernicus')
                setattr(ncdf, '_FillValue', fill_value)

                Vtime = ncdf.createVariable('time', 'f', ('time',))
                Vtime[:] = np.array(lst_times).astype('float32')

                newvar = ncdf.createVariable(variable, 'f', ('time', 'lat', 'lon'))
                newvar[:] = output_tdata.astype('float32')
                setattr(newvar, 'longname', variable.replace('_', ' '))
            logging.info('Created LST: {}'.format(os.path.basename(filename)))
        except FileNotFoundError:
            continue
        except:
            print('Readin failed on day ', doy, year)
