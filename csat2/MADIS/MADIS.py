from csat2 import locator, misc
from netCDF4 import Dataset, num2date, chartostring
import os
import requests
from tqdm import tqdm
import logging

log = logging.getLogger(__name__)


def readin(product, year, doy, hour, sds, only_ships=False):
    def _file_output(ncdf, sds, only_ships):
        ncdf.set_auto_mask(False)
        if only_ships:
            # value 1 is floating buoy or ship
            # value 0 is moored buoy or CMAN
            ship_flag = ncdf.variables["dataPlatformType"][:] == 1
        for name in sds:
            output[name] = ncdf.variables[name][:]
            if only_ships:
                output[name] = output[name][ship_flag]
            if name == "timeObs":
                output[name] = num2date(output[name], ncdf.variables[name].units)
            if ncdf.variables[name].dtype.str.startswith("|S"):
                output[name] = chartostring(output[name])
        return output

    filename = locator.search("MADIS", product, year=year, doy=doy, hour=hour)[0]

    output = {}
    if product == "metar" or filename.endswith(".gz"):
        # Gzipped output stored, unzip to temporary file
        import gzip

        with gzip.open(filename) as gz:
            with Dataset("dummy", mode="r", memory=gz.read()) as ncdf:
                return _file_output(ncdf, sds, only_ships)
    else:
        # No longer any need to un-gzip netcdf files as can use memory-mapped
        with Dataset(filename) as ncdf:
            return _file_output(ncdf, sds, only_ships)


def variables(product, year, doy, hour):
    filename = locator.search("MADIS", product, year=year, doy=doy, hour=hour)[0]
    if product == "metar" or filename.endswith(".gz"):
        import gzip

        with gzip.open(filename) as gz:
            with Dataset("dummy", mode="r", memory=gz.read()) as ncdf:
                for name in ncdf.variables.keys():
                    try:
                        lname = ncdf.variables[name].long_name
                    except:
                        lname = ""
                    try:
                        units = ncdf.variables[name].units
                    except:
                        units = ""
                    print(name, lname, units)
    else:
        with Dataset(filename) as ncdf:
            for name in ncdf.variables.keys():
                try:
                    lname = ncdf.variables[name].long_name
                except:
                    lname = ""
                try:
                    units = ncdf.variables[name].units
                except:
                    units = ""
                print(name, lname, units)


def check(product, year, doy, hour):
    filename = locator.search("MADIS", product, year=year, doy=doy, hour=hour)
    if len(filename) > 0:
        return True
    else:
        return False


def download(product, year, doy, hour, force_redownload=False):
    if check(product, year, doy, hour) and not force_redownload:
        return

    _, mon, day = misc.time.doy_to_date(year, doy)
    if product.endswith("-nrt"):
        nproduct = product.replace("-nrt", "")
        url = "ftp://madis-data.ncep.noaa.gov/point/{product}/netcdf/{year}{mon:0>2}{day:0>2}_{hour:0>2}00.gz".format(
            product=nproduct, year=year, mon=mon, day=day, hour=hour
        )
    else:
        url = "ftp://madis-data.ncep.noaa.gov/archive/{year}/{mon:0>2}/{day:0>2}/point/{product}/netcdf/{year}{mon:0>2}{day:0>2}_{hour:0>2}00.gz".format(
            product=product, year=year, mon=mon, day=day, hour=hour
        )

    local_folder = locator.get_folder("MADIS", product, year=year, doy=doy)
    log.debug("Downloading to: " + local_folder)
    try:
        os.makedirs(local_folder)
    except FileExistsError:
        pass

    newfile = local_folder + "/" + url.split("/")[-1]

    if force_redownload:
        try:
            os.remove(newfile)
        except FileNotFoundError:
            pass

    if not os.path.exists(newfile):
        import urllib.request

        urllib.request.urlretrieve(url, newfile)
    else:
        raise ValueError("File exists")
