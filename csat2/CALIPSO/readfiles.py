from csat2 import locator
import csat2.misc.time
import csat2.misc.geo
import csat2.misc.hdf
from pyhdf.error import HDF4Error
import numpy as np
import logging
from datetime import datetime, timedelta

DEFAULT_COLLECTION = 'V4-10'

print('Default CALIPSO collection: {}'.format(DEFAULT_COLLECTION))

########################################################################
# Functions for reading in Satellite data, CloudSat data
########################################################################
#
#  ERI Gryspeerdt, Imperial College London, 2022
#   - Based on code originally written for csat library
########################################################################


def files_available(product, year, doy, col=DEFAULT_COLLECTION):
    filenames = locator.search('CALIPSO',
                               product, year=year, doy=doy,
                               hour='**', minute='**',
                               col=col)
    return sorted(filenames)


def variables_available(product, year, doy, col=DEFAULT_COLLECTION):
    filename = locator.search('CALIPSO',
                              product, year=year, doy=doy,
                              hour='**', minute='**',
                              col=col)
    if len(filename) == 0:
        raise IOError('No files for that date')
    filename = np.sort(filename)[0]
    data = csat2.misc.hdf.read_hdf4(filename, vdata=True)
    return data


def filename_to_datetime(fname):
    return datetime.strptime(fname.split('.')[-2][:-2], '%Y-%m-%dT%H-%M-%S')


def get_filename_by_time(dtime):
    geometa_file = locator.search('CALIPSO', 'GEOMETA',
                                  year=dtime.year)[0]
    timestr = dtime.strftime('%Y-%m-%dT%H-%M-%S')
    with open(geometa_file, 'r') as f:
        #Gets the id of the first orbit after the time string
        fname = f.readline()
        newline = f.readline()
        while newline.split('.')[-2][:-2] < timestr:
            fname = newline
            newline = f.readline()
    return fname.strip()


def readin_calipso_curtain(product, year, doy, hour, minute, second, sds=None,
                           col=DEFAULT_COLLECTION):
    '''Reads in data from CALIPSO curtain files.'''
    _, mon, day = csat2.misc.time.doy_to_date(year, doy)

    filename = locator.search(
        'CALIPSO', product, year=year,
        doy=doy, mon=mon, day=day, hour=hour, minute=minute, col=col)[0]

    if isinstance(sds, str):
        sds = [sds]

    return csat2.misc.hdf.read_hdf4(filename, sds, vdata=True)
