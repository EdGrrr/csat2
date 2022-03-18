from csat2 import locator
import csat2.misc.time
import csat2.misc.geo
import csat2.misc.hdf
from pyhdf.error import HDF4Error
import numpy as np
import logging
from datetime import datetime, timedelta

DEFAULT_COLLECTION = 'P2_R05'

print('Default CloudSat collection: {}'.format(DEFAULT_COLLECTION))

########################################################################
# Functions for reading in Satellite data, CloudSat data
########################################################################
#
#  ERI Gryspeerdt, Imperial College London, 2022
#   - Based on code originally written for csat library
########################################################################


def files_available(product, year, doy, rev):
    filenames = locator.search('CLOUDSAT',
                               product, year=year, doy=doy,
                               rev=rev, orbit='*****')
    return sorted(filenames)


def variables_available(product, year, doy, rev):
    filename = locator.search('CLOUDSAT',
                              product, year=year, doy=doy,
                              rev=rev, orbit='*****')
    if len(filename) == 0:
        raise IOError('No files for that date')
    filename = np.sort(filename)[0]
    data = csat2.misc.hdf.read_hdf4(filename, vdata=True)
    return data


def available_orbits(product, year, doy, rev):
    filename = locator.search('CLOUDSAT', product,
                              year=year, doy=doy, rev=rev, orbit='*****')
    if len(filename) == 0:
        raise IOError('No files for that date')
    filename = np.sort(filename)
    return list(map(lambda x: os.path.basename(x).split('_')[1], filename))


def get_orbit_date_approx(orbit):
    index_date = datetime(2007, 1, 1)
    index_orbit = 3607
    orbits_per_day = 14.56424
    test_orb = (index_date +
                timedelta(days=(
                    (orbit-index_orbit) /
                    orbits_per_day)))
    return test_orb.year, test_orb.timetuple().tm_yday


def get_orbit_filename(orbit):
    year, _ = get_orbit_date_approx(orbit)
    try:
        geometa_file = locator.search('CLOUDSAT', 'GEOMETA',
                                      year=year)[0]
        with open(geometa_file, 'r') as f:
            fname = next(obj for obj in f if int(obj[14:19]) == int(orbit))
    except StopIteration:
        geometa_file = locator.search('CLOUDSAT', 'GEOMETA',
                                      year=year+1)[0]
        with open(geometa_file, 'r') as f:
            fname = next(obj for obj in f if int(obj[14:19]) == int(orbit))
    return fname.strip()


def get_orbit_datetime(orbit):
    return filename_to_datetime(get_orbit_filename(orbit))


def filename_to_datetime(fname):
    return datetime.strptime(fname.split('_')[0], '%Y%j%H%M%S')


def get_orbit_date(orbit):
    year, _ = get_orbit_date_approx(orbit)
    try:
        geometa_file = locator.search('CLOUDSAT', 'GEOMETA',
                                      year=year)[0]
        with open(geometa_file, 'r') as f:
            orbtime = next(obj for obj in f if int(obj[14:19]) == int(orbit)).split('_')[0]
    except StopIteration:
        geometa_file = locator.search('CLOUDSAT', 'GEOMETA',
                                      year=year+1)[0]
        with open(geometa_file, 'r') as f:
            orbtime = next(obj for obj in f if int(obj[14:19]) == int(orbit)).split('_')[0]        
    return int(orbtime[:4]), int(orbtime[4:7])


def get_orbit_by_time(dtime):
    geometa_file = locator.search('CLOUDSAT', 'GEOMETA',
                                  year=dtime.year)[0]
    timestr = dtime.strftime('%Y%j%H%M%S')
    with open(geometa_file, 'r') as f:
        #Gets the id of the first orbit after the time string
        orb = int(next(obj for obj in f if obj[:13] > timestr).split('_')[1])
    # Subtract 1 to get the orbit containing the timestamp
    return orb - 1


def readin_cloudsat_curtain(product, year=None, doy=None, orbit=None, sds=None,
                            col=DEFAULT_COLLECTION):
    '''Reads in data from CloudSat curtain files.'''
    if year is None or doy is None:
        year, doy = get_orbit_date(orbit)
    _, mon, day = csat2.misc.time.doy_to_date(year, doy)

    filename = locator.search(
        'CLOUDSAT', product, year=year,
        doy=doy, mon=mon, day=day, col=col, orbit=orbit)[0]

    if isinstance(sds, str):
        sds = [sds]

    return csat2.misc.hdf.read_hdf4(filename, sds, vdata=True)
