from csat2 import locator
import csat2.misc.time
from .readfiles import DEFAULT_COLLECTION
import sys
import os
import os.path
import json
import logging
import json
import ftplib
import pysftp

with open(os.path.expandvars('${HOME}/.csat2/cloudsat_cira.json')) as f:
    cloudsat_auth = json.load(f)

cloudsat_server = 'www.cloudsat.cira.colostate.edu'

# product = '2B-GEOPROF-LIDAR'
# year, doy = 2008, 100
# orbits = [10366]
# DEFAULT_COLLECTION = 'P2_R05'


def download_file_locations(product, year, doy, orbits=None,
                            col=DEFAULT_COLLECTION):
    _, mon, _ = csat2.misc.time.doy_to_date(year, doy)
    
    folder = f'Data/{product}.{col}/{year}/{doy:0>3}/'

    srv = pysftp.Connection(host=cloudsat_server,
                            username=cloudsat_auth['username'],
                            private_key=cloudsat_auth['keyfile'])
    srv.chdir(folder)

    files = srv.listdir()
    srv.close()

    if orbits is not None:
        return [f for f in files if int(f.split('_')[1]) in orbits]
    else:
        return files


def download(product, year, doy, orbits=None, col=DEFAULT_COLLECTION):
    files = download_file_locations(product, year, doy, orbits, col)

    local_folder = locator.get_folder('CLOUDSAT', product,
                                      year=year, doy=doy,
                                      col=col)

    try:
        os.makedirs(local_folder)
    except FileExistsError:
        pass

    if len(files) == 0:
        raise ValueError('No files for {} on {}, {}'.format(product, year, doy))

    # Assumes there is only one valid folder for this collection of files
    folder = f'Data/{product}.{col}/{year}/{doy:0>3}/'

    srv = pysftp.Connection(host=cloudsat_server,
                            username=cloudsat_auth['username'],
                            private_key=cloudsat_auth['keyfile'])
    srv.chdir(folder)

    for filename in files:
        newfile = local_folder + '/' + os.path.basename(filename)
        if not os.path.exists(newfile):
            srv.get(filename, newfile)
            # with open(newfile, 'w+b') as fh:
            #     ftp.retrbinary(f'RETR {filename}', fh.write)
            print(f'{filename} complete')
        else:
            logging.info('Skipping {}'.format(os.path.basename(filename)))
    srv.close()


def check(product, year, doy, orbit, col=DEFAULT_COLLECTION):
    '''Does a product already exist for a specific time/date/collection'''
    filename = locator.search(
        'CLOUDSAT',
        product,
        year=year,
        doy=doy,
        col=col,
        orbit=orbit)
    if len(filename) == 1:
        return True
    else:
        return False
