from calendar import c
from functools import cached_property
from gettext import find
import logging

from tqdm import tqdm
from csat2 import GOES
from csat2.GOES.granule import GOESLocator
import re
import json
import os
import csat2
import requests
import datetime
import warnings
import numpy as np
import zipfile
import shutil
import satpy  # new requirement -->  a pretty mature codebase at this stage...
import eumdac
from satpy.scene import (
    Scene,
)  # new requirement -->  this is a simple package that provides a simple interface to the EUMETSAT Data Centre API. Could probably reverse engineer these parts of the code to work with the EUMETSAT API directly
from satpy import find_files_and_readers
import xarray as xr

CREDENTIALFILE = os.environ["HOME"] + "/.csat2/eumetsat_key.json"

with open(CREDENTIALFILE, "r") as jfile:
    credentials = json.load(jfile)


def get_chunk_files(product):
    return sorted([e for e in product.entries if "FD--CHK-BODY" in e], key=lambda e: e[-7:-3])

def download_fci(directory, selected_product: eumdac.product.Product, load_chunks_iis=None):
    directory = os.path.join(directory, str(selected_product))
    try:
        os.makedirs(directory)
    except FileExistsError:
        pass
    
    chunks = get_chunk_files(selected_product)
    if load_chunks_iis is not None:
        chunks = [chunks[i] for i in load_chunks_iis]
    for chunk in tqdm(chunks):
        try:
            with selected_product.open(chunk) as fsrc, open(
                os.path.join(directory, fsrc.name), mode="wb"
            ) as fdst:
                shutil.copyfileobj(fsrc, fdst)
            # print(f"Download of product {selected_product} finished.")
        except eumdac.product.ProductError as error:
            print(
                f"Error related to the product '{selected_product}' while trying to download it: '{error.msg}'"
            )
        except requests.exceptions.ConnectionError as error:
            print(f"Error related to the connection: '{error.msg}'")
        except requests.exceptions.RequestException as error:
            print(f"Unexpected error: {error}")


def get_eumetsat_ids(dt, high_res=False, limit=1, max_interval_min=10):
    high_res_collection_id = "EO:EUM:DAT:0665"  # 0 degree high res
    normal_res_collection_id = "EO:EUM:DAT:0662"  # 0 degree normal res
    collection_id = high_res_collection_id if high_res else normal_res_collection_id

    token = eumdac.AccessToken(
        (credentials["consumer_key"], credentials["consumer_secret"])
    )
    datastore = eumdac.DataStore(token)
    collection = datastore.get_collection(collection_id)

    result = collection.search(
        dtstart=dt,
        dtend=dt + datetime.timedelta(minutes=max_interval_min),
    )
    if result.total_results == 0:
        raise ValueError("No products found")
    if result.total_results > limit:
        warnings.warn(f"{result.total_results} found")
    return result.first()


def get_eumdac_product(gran):
    products = get_eumetsat_ids(gran.datetime(), high_res=gran.resolution == "HRFI")
    return products



def read_fci_scene(fname):
    files = find_files_and_readers(
        base_dir=fname, reader="fci_l1c_nc"
    )
    return Scene(filenames=files)

    
def read_fci_file(fname, datasets):
    scn = read_fci_scene(fname)
    scn.load(datasets, upper_right_corner="NE")
    ds = scn.to_xarray_dataset().compute() # doing dask and lazy loading would be nice - but needs us to stop the tempfiles.

    return ds

def get_available_products(fname):
    scn = read_fci_scene(fname)
    avail_data = scn.available_dataset_names()
    avail_comp = scn.available_composite_names()
    avail_data.extend(avail_comp)
    return avail_data


# def fci_datafile_name(
#     high_res, process_time, start_time, end_time, repeat_cycle_in_day
# ):
#     imager_n = 1  # i.e. MTGi (Meteosat-12 is MTG1)
#     subtype = "HRFI" if high_res else "FDHSI"
#     return f"W_XX-EUMETSAT-Darmstadt,IMG+SAT,MTI{imager_n}+FCI-1C-RRAD-{subtype}-FD--x-x---x_C_EUMT_{process_time.strftime('%Y%m%d%H%M%S')}_IDPFI_OPE_{start_time.strftime('%Y%m%d%H%m%S')}_{end_time.strftime('%Y%m%d%H%m%S')}_N__O_{repeat_cycle_in_day}_0000.nc"


class Granule(GOES.Granule):
    # Valid satellite names are MTG1, MTG2 etc.
    # Valid resolutions are HRFI (high spatial res) and (FDHSI) full disk high spectral res

    inc_minutes = {"FD": 10}

    def __init__(self, *args, resolution="FDHSI", **kwargs):
        super().__init__(*args, **kwargs)
        # FCI files (like SEVIRI) cannot be band separated, so all bands are read together
        # and stored. There is not an easy way around this other than unzipping
        # the files
        self.resolution = resolution
        self.data = None

    def download(self, force_redownload=False):
        local_folder = csat2.locator.get_folder(
            "FCI",
            "L1c",
            year=self.year,
            doy=self.doy,
            hour=self.hour,
            area=self.area,
            sat=self.sat,
            resolution=self.resolution,
        )

        try:
            os.makedirs(local_folder)
        except FileExistsError:
            pass

        product = get_eumdac_product(self)

        download_fci(local_folder, product)

    @classmethod
    def fromtext(cls, gran_text):
        m = re.search(
            "M(?P<sat>...)\\.(?P<year>....)(?P<doy>...)\\.(?P<hour>..)(?P<minute>..)\\.(?P<area>.*)\\.(?P<res>.*)",
            gran_text,
        )
        return cls(
            sat="M" + m.group("sat"),
            area=m.group("area"),
            year=int(m.group("year")),
            doy=int(m.group("doy")),
            hour=int(m.group("hour")),
            minute=int(m.group("minute")),
            resolution=m.group("res"),
        )

    def next(self, number=1, only_downloaded=False, only_exisiting=False):
        """Increment image name"""
        dt = self.datetime()
        dt += datetime.timedelta(minutes=number * self.inc_minutes[self.area])
        year, doy = csat2.misc.time.date_to_doy(dt.year, dt.month, dt.day)
        return self.__class__(
            sat=self.sat,
            area=self.area,
            year=year,
            doy=doy,
            hour=dt.hour,
            minute=dt.minute,
            resolution=self.resolution,
            # Reasonable assumption that coordinates remain the same
            locator=self.locator,
        )

    def get_filename(self, product="L1b", *args, **kwargs):
        """Selects an image with a mane that puts it within the increment
        timestep."""
        filenames = csat2.locator.search(
            "FCI",
            "L1c",
            year=self.year,
            doy=self.doy,
            hour=self.hour,
            area=self.area,
            sat=self.sat,
            minute="**",
            resolution=self.resolution,
        )
        product = get_eumdac_product(self)

        files = [f for f in filenames if str(product) in f]
        if len(files) == 0:
            raise IndexError("No matching file")
        if len(files) > 1:
            raise IndexError("Non-unique filename")
        return files[0]

    def _read_file(self, datasets):
        fname = self.get_filename()
        self.data = read_fci_file(fname, datasets)

    def get_band_radiance(self, channel):
        datasets = [channel] if isinstance(channel, str) else channel
        self._read_file(datasets=datasets)
        return self.data[channel]
    
    @cached_property
    def available_products(self):
        fname = self.get_filename()
        return get_available_products(fname)
        
    def get_band_bt():
        raise NotImplementedError("Not implemented yet")
        

    # def get_shape(self, channel=None, product="L1b", *args, **kwargs):
    #     return (3712, 3712)

    # # def get_viewangles(self, channel=None, product="L1b"):
    # #     self._check_locator(product)
    # #     return np.array(np.meshgrid(self.locator.x, self.locator.y))

    # # def get_llcoord(self, channel=None, product="L1b", mode="*"):
    # #     self._check_locator(product)
    # #     return self.locator.x[0], self.locator.y[-1]

    # # def _check_locator(self, channel=None, product="L1b", *args, **kwargs):
    # #     if not self.locator:
    # #         self.locator = SeviriLocator(self.get_filename(product=product))

    def get_satpy_scene(self):
        fname = self.get_filename()
        scn =  read_fci_scene(fname)
        return scn

    def check(self, chunks=None):
        try:
            fname = self.get_filename()
        except IndexError:
            return False
        
        possible_chunks = get_chunk_files(get_eumdac_product(self))
        if chunks is None:
            chunks = list(range(len(possible_chunks)))
        
        files_to_check = [os.path.join(fname, possible_chunks[i]) for i in chunks]
        return all([os.path.exists(f) for f in files_to_check])
