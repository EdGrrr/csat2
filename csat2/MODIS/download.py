from __future__ import print_function, division
import csat2
from csat2 import locator
from .readfiles import DEFAULT_COLLECTION
import sys
import os
import os.path
import json
import logging


USERAGENT = 'csat2/download_v{}'.format(csat2.__version__).replace('\r', '')
TOKENFILE = os.environ['HOME']+'/.csat2/laadsdaacrc'


def _geturl(url, token=None, out=None):
    headers = {'user-agent': USERAGENT}
    if token is not None:
        headers['Authorization'] = 'Bearer ' + token
        try:
            import ssl
            CTX = ssl.SSLContext(ssl.PROTOCOL_TLSv1_2)
            from urllib.request import urlopen, Request, HTTPError, URLError
            try:
                response = urlopen(Request(url, headers=headers), context=CTX)
                if out is None:
                    return response.read().decode('utf-8')
                else:
                    if response.getcode() != 200:
                        raise HTTPError('No data available')
                    length = response.getheader('content-length')
                    if length:
                        length = int(length)
                        blocksize = max(4096, length//100)
                    else:
                        raise HTTPError('No content-length from LAADS')
                    size = 0
                    while True:
                        buf1 = response.read(blocksize)
                        if not buf1:
                            break
                        out.write(buf1)
                        size += len(buf1)
                        if length:
                            print('Downloading {} {:6.2f}MB {:6.0f}%\r'.format(
                                os.path.basename(url), length/(1024*1024), 100*size/length), end='')
                    print()
            except HTTPError as e:
                raise IOError('Data not available: {}'.format(e))
            except URLError:
                raise IOError('Data not available')

        except AttributeError:
            # OS X Python 2 and 3 don't support tlsv1.1+ therefore... curl
            import subprocess
            try:
                args = ['curl', '--fail', '-sS', '-L', '--get', url]
                for (k, v) in headers.items():
                    args.extend(['-H', ': '.join([k, v])])
                if out is None:
                    # python3's subprocess.check_output returns stdout as a byte string
                    result = subprocess.check_output(args)
                    return result.decode('utf-8') if isinstance(result, bytes) else result
                else:
                    subprocess.call(args, stdout=out)
            except subprocess.CalledProcessError as e:
                logging.warn('curl GET error message: %' +
                             (e.message if hasattr(e, 'message') else e.output),
                             file=sys.stderr)
        return None


def get_token():
    try:
        with open(TOKENFILE) as f:
            return f.read().strip()
    except FileNotFoundError:
        raise IOError(
            'Place your LAADS API key in the file {}'.format(TOKENFILE))


def download(product, year, doy, times=None, col=DEFAULT_COLLECTION, force_redownload=False):
    if times:
        '''If times are supplied, check to see if they all exist. If not,
        assume that we are trying to download lots of things and ask 
        the server'''
        times = [t for t in times if not check(product, year, doy, t, col=col)]
        if len(times) == 0: # All files exist
            return
        
    laads_folder = ('https://ladsweb.modaps.eosdis.nasa.gov/' +
                    'archive/allData/{col}/{product}/{year}/{doy:0>3}'.format(
                        col=col,
                        product=product,
                        year=year,
                        doy=doy))

    local_folder = locator.get_folder('MODIS', product,
                                      year=year, doy=doy,
                                      collection=col)

    try:
        os.makedirs(local_folder)
    except FileExistsError:
        pass

    logging.debug(laads_folder)
    token = get_token()
    files = [a['name'] for a in
             json.loads(_geturl(laads_folder+'.json', token))]
    if times is not None:
        files = [a for a in files if str(a.split('.')[2]) in times]
    
    if len(files) == 0:
        raise ValueError(
            'No files for {} on {}, {}'.format(product, year, doy))

    for filename in files:
        laads_loc = laads_folder + '/' + filename
        newfile = local_folder + '/' + filename
        if (not os.path.exists(newfile)) or force_redownload:
            if force_redownload:
                os.remove(newfile)
            with open(newfile, 'w+b') as fh:
                _geturl(laads_loc, token, fh)
        else:
            logging.info('Skipping {}'.format(filename))


def download_geometa(year, doy, sat, col=DEFAULT_COLLECTION, force_redownload=False):
    _, mon, day = csat2.misc.time.doy_to_date(year, doy)
    laads_file = (
        'https://ladsweb.modaps.eosdis.nasa.gov/' +
        'archive/geoMeta/' +
        '{col}/{sat}/{year}/M{si}D03_{year}-{mon:0>2}-{day:0>2}.txt'.format(
            col=col,
            sat=sat,
            year=year,
            mon=mon,
            day=day,
            si={'AQUA': 'Y',
                'TERRA': 'O'}[sat]))

    local_file = locator.format_filename('MODIS', 'GEOMETA',
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
            _geturl(laads_file, token, fh)
    else:
        logging.info('Skipping {}'.format(os.path.basename(local_file)))
    
            

def check(product, year, doy, time=None, col=DEFAULT_COLLECTION):
    '''Does a product already exist for a specific time/date/collection'''
    filename = locator.search(
        'MODIS',
        product,
        year=year,
        doy=doy,
        collection=col,
        time=time)
    if len(filename) == 1:
        return True
    else:
        return False
