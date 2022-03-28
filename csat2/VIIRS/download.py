from csat2 import locator
import csat2.misc.time
from .readfiles import DEFAULT_COLLECTION
import sys
import os
import os.path
import json
import logging
from csat2.download.earthdata import get_token, geturl
log = logging.getLogger(__name__)

    
def download_file_locations(product, year, doy, times=None,
                            col=DEFAULT_COLLECTION):
    base_url = 'https://ladsweb.modaps.eosdis.nasa.gov/'
    laads_folder = (base_url +
                    f'archive/allData/{col}/{product}/{year}/{doy:0>3}/')

    token = get_token()
    files = [laads_folder + a['name'] for a in
             json.loads(geturl(laads_folder+'.json', token).decode('utf-8'))]
    return files


def download_file_locations_nrt(product, year, doy, times=None,
                                col=DEFAULT_COLLECTION):
    base_url = 'https://nrt4.modaps.eosdis.nasa.gov'
    laads_folder = (base_url +
                    f'/api/v2/content/details/allData/{col}/{product}/{year}/{doy:0>3}/')

    token = get_token()
    files = [base_url+a['downloadsLink'] for a in
             json.loads(geturl(laads_folder+'?fields=all&formats=json', token).decode('utf-8'))['content']]
    return files


def download(product, year, doy, times=None, col='5110'):
    if product.endswith('NRT'):
        files = download_file_locations_nrt(product, year, doy, times, col)
    else:
        files = download_file_locations(product, year, doy, times, col)

    local_folder = locator.get_folder('VIIRS', product,
                                      year=year, doy=doy,
                                      collection=col)

    try:
        os.makedirs(local_folder)
    except FileExistsError:
        pass

    if len(files) == 0:
        raise ValueError('No files for {} on {}, {}'.format(product, year, doy))

    if times is not None:
        files = [a for a in files if str(os.path.basename(a).split('.')[2]) in times]

    for filename in files:
        newfile = local_folder + '/' + os.path.basename(filename)
        if not os.path.exists(newfile):
            with open(newfile, 'w+b') as fh:
                geturl(filename, get_token(), fh)
        else:
            log.info('Skipping {}'.format(os.path.basename(filename)))


def download_geometa(year, doy, sat, col=DEFAULT_COLLECTION, force_redownload=False):
    _, mon, day = csat2.misc.time.doy_to_date(year, doy)
    laads_file = (
        'https://ladsweb.modaps.eosdis.nasa.gov/' +
        'archive/geoMetaJPSS/' +
        '{col}/{sat}/{year}/{si}03MOD_{year}-{mon:0>2}-{day:0>2}.txt'.format(
            col=col,
            sat=sat,
            year=year,
            mon=mon,
            day=day,
            si={'NPP': 'VNP',
                'JPSS1': 'VJ1'}[sat]))

    local_file = locator.format_filename('VIIRS', 'GEOMETA',
                                         year=year, doy=doy,
                                         sat=sat)

    try:
        os.makedirs(os.path.dirname(local_file))
    except FileExistsError:
        pass

    token = get_token()
    if (not os.path.exists(local_file)) or force_redownload:
        if force_redownload:
            os.remove(local_file)
        with open(local_file, 'w+b') as fh:
            geturl(laads_file, token, fh)
    else:
        log.info('Skipping {}'.format(os.path.basename(local_file)))


def check(product, year, doy, time, col=DEFAULT_COLLECTION):
    '''Does a product already exist for a specific time/date/collection'''
    filename = locator.search(
        'VIIRS',
        product,
        year=year,
        doy=doy,
        collection=col,
        time=time)
    if len(filename) == 1:
        return True
    else:
        return False
