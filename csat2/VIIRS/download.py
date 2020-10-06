from __future__ import print_function, division
from csat2 import locator
import csat2.misc.time
import sys
import os
import os.path
import json
import logging

USERAGENT = 'csat2/download_v0.1' + \
    sys.version.replace('\n', '').replace('\r', '')
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
                fh = urlopen(Request(url, headers=headers), context=CTX)
                if out is None:
                    return fh.read().decode('utf-8')
                else:
                    length = fh.getheader('content-length')
                    if length:
                        length = int(length)
                        blocksize = max(4096, length//100)
                    else:
                        raise HTTPError('No content-length from LAADS')
                    size = 0
                    while True:
                        buf1 = fh.read(blocksize)
                        if not buf1:
                            break
                        out.write(buf1)
                        size += len(buf1)
                        if length:
                            print('Downloading {} {:6.2f}MB {:6.0f}%\r'.format(os.path.basename(url), length/(1024*1024), 100*size/length), end='')
                    print()
            except HTTPError:
                return IOError('Data not available')
            except URLError:
                return IOError('Data not available')

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
                print('curl GET error message: %' +
                      (e.message if hasattr(e, 'message') else e.output),
                      file=sys.stderr)
        return None


def get_token():
    try:
        with open(TOKENFILE) as f:
            return f.read().strip()
    except:
        raise IOError('Place your LAADS API key in {}'.format(TOKENFILE))


def download(product, year, doy, times=None, col='5110'):
    laads_folder = ('https://ladsweb.modaps.eosdis.nasa.gov/' +
                    'archive/allData/{}/{}/{}/{:0>3}'.format(
                        col,
                        product,
                        year,
                        doy))

    local_folder = locator.get_folder('VIIRS', product,
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

    if len(files) == 0:
        raise ValueError('No files for {} on {}, {}'.format(product, year, doy))

    if times is not None:
        files = [a for a in files if str(a.split('.')[2]) in times]

    for filename in files:
        laads_loc = laads_folder + '/' + filename
        newfile = local_folder + '/' + filename
        if not os.path.exists(newfile):
            with open(newfile, 'w+b') as fh:
                _geturl(laads_loc, token, fh)
        else:
            logging.info('Skipping {}'.format(filename))
