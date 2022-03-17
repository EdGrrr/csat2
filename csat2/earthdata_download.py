from __future__ import print_function, division
import csat2
import sys
import os
import os.path
import logging


USERAGENT = 'csat2/download_v{}'.format(csat2.__version__).replace('\r', '')
TOKENFILE = os.environ['HOME']+'/.csat2/laadsdaacrc'


def geturl(url, token=None, out=None):
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
