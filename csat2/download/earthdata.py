import csat2
import sys
import os
import os.path
import logging
import requests
from tqdm import tqdm

USERAGENT = 'csat2/download_v{}'.format(csat2.__version__).replace('\r', '')
TOKENFILE = os.environ['HOME']+'/.csat2/laadsdaacrc'


def geturl(url, token=None, out=None):
    session = requests.Session()
    response = session.get(url, stream=True)
    headers = {'user-agent': USERAGENT}
    if token is not None:
        headers['Authorization'] = 'Bearer ' + token
    session = requests.Session()
    session.headers = headers

    if response.status_code != 200:
        raise ValueError(f'Download failed status:{response.status_code}')
    total = int(response.headers.get('content-length', 0))
    if out:
        with tqdm(
                desc=os.path.basename(url),
                total=total,
                unit='iB',
                unit_scale=True,
                unit_divisor=1024,
        ) as bar:
            blocksize = max(4096, total//100)
            for data in response.iter_content(chunk_size=blocksize):
                size = out.write(data)
                bar.update(size)
    else:
        return response.content


def get_token():
    try:
        with open(TOKENFILE) as f:
            return f.read().strip()
    except FileNotFoundError:
        raise IOError(
            'Place your LAADS API key in the file {}'.format(TOKENFILE))
