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
import typing
from html.parser import HTMLParser
from requests import Session
from tqdm import tqdm

session = Session()


class LinkParser(HTMLParser):
    '''Based on example code from NASA Langley ASDC'''
    def __init__(self, *args, **kwargs):
        self.hrefs = []

    def handle_starttag(self, tag, attrs):
        if tag == "a":
            for (attr, value) in attrs:
                if attr == "href":
                    self.hrefs.append(value)

    def get_hrefs(self, text):
        self.feed(text)
        return self.hrefs


def asdc_download_file(url, output_fname=None):
    response = session.get(url, stream=True)
    if response.status_code != 200:
        raise ValueError(f'Download failed status:{response.status_code}')
    total = int(response.headers.get('content-length', 0))
    if output_fname:
        with open(output_fname, 'wb') as f, tqdm(
                desc=os.path.basename(output_fname),
                total=total,
                unit='iB',
                unit_scale=True,
                unit_divisor=1024,
        ) as bar:
            for data in response.iter_content(chunk_size=1024):
                size = f.write(data)
                bar.update(size)
    else:
        return response.content


def download_file_locations(product, year, mon,
                            col=DEFAULT_COLLECTION, full_location=False):

    folder = f'https://asdc.larc.nasa.gov/data/CALIPSO/{product}-Standard-{col}/{year}/{mon:0>2}/'

    is_child_href = lambda href: (
        not href.startswith('http') and
        not (href.startswith('/') or
             href.startswith('#') or
             href.startswith('mailto')))

    linkpage = asdc_download_file(folder).decode('utf-8')
    parser = MyHTMLParser()
    hrefs = [link for link in parser.get_hrefs(linkpage) if is_child_href(link) and link.endswith('hdf')]

    #Remove duplicates
    hrefs = list(set(hrefs))
    hrefs.sort()

    return hrefs


def download(product, year, doy, hour, minute, second, daynight, col=DEFAULT_COLLECTION):
    local_folder = locator.get_folder('CALIPSO', product,
                                      year=year, doy=doy,
                                      col=col)

    _, mon, day = csat2.misc.time.doy_to_date(year, doy)

    try:
        os.makedirs(local_folder)
    except FileExistsError:
        pass

    asdc_url = (f'https://asdc.larc.nasa.gov/data/CALIPSO/{product}-Standard-{col}/{year}/{mon:0>2}/' +
                f'CAL_{product}-Standard-{col}.{year}-{mon:0>2}-{day:0>2}T{hour:0>2}-{minute:0>2}-{second:0>2}Z{daynight}.hdf')

    localname = os.path.join(local_folder, os.path.basename(asdc_url))

    asdc_download_file(asdc_url, output_fname=localname)


def check(product, year, doy, hour, minute, col=DEFAULT_COLLECTION):
    '''Does a product already exist for a specific time/date/collection.
Note that the seconds is not required to uniquely specify a CALIPSO granule'''
    filename = locator.search(
        'CALIPSO',
        product,
        year=year,
        doy=doy,
        hour=hour,
        minute=minute,
        col=col)
    if len(filename) == 1:
        return True
    else:
        return False
