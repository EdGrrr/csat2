import csat2
from csat2 import locator
import csat2.misc.time
import sys
import os
import os.path
import logging
import requests
from .utils import check_valid_args



def get_file_from_url(url, local_path):
    '''helper function to download a file from a URL to a local path'''
    try:
        response = requests.get(url, stream=True)
        response.raise_for_status()  # Raise an error for bad responses
        with open(local_path, 'wb') as file:
            for chunk in response.iter_content(chunk_size=8192):
                file.write(chunk)
        logging.info(f"Downloaded {local_path} from {url}")
    except requests.exceptions.RequestException as e:
        logging.error(f"Failed to download {url}: {e}")
        sys.exit(1)
    except Exception as e:
        logging.error(f"Failed to download {url}: {e}")
        raise



def download_files(year:int ,doy:int,times:int,product: str = 'isccp-basic',version: str = 'hgg',force_redownload = False):
    '''Download file locations for ISCCP data.
    Args: 
        product (str): either 'isccp' or 'isccp-basic', basic is just a subset of the full isccp data and generally contains all the data you need, basic is the default
        version (str): one of hgg,hgh,hgm,hgx (hgg is the default, hgx is unavailable for ISCCP Basic)
        year (int): the year of the data to download
        doy (int): the day of year of the data to download, 1-365 (or 366 for leap years)
        times (int): the times are UTC in three hourly intervals, e.g. 0,3,6,9,12,15,18,21

    Returns:
        downloads the ISCCP data files to the local filesystem in a structured way.
    '''
    check_valid_args(year, doy, times, product, version)  # Validate the input arguments
    logging.info(f"Downloading ISCCP data for {year}-{doy:03d} at {times:02d}00 UTC, product: {product}, version: {version}")
    

    base_url = f"https://www.ncei.noaa.gov/data/international-satellite-cloud-climate-project-isccp-h-series-data/access"

    product_uppercase = {'isccp-basic': 'ISCCP-Basic', 'isccp': 'ISCCP'} ## converts the product name as it appears in two different ways, uppercase version in the remote filename and the local file naming convention


    version_uppercase = version.upper()  # Convert version to uppercase as it appears uppercase in the remote filename
    local_dir = locator.get_folder(product_uppercase[product],version_uppercase,year=year,doy = doy)

    year, month,dom = csat2.misc.time.doy_to_date(year, doy) # returns a tuple of integer values for year, month and day of month
    doy_str = f"{doy:03d}"  # Format day of year as a three-digit string
    time_str = f"{times:02d}" +'00'   # Format time as a four digit string (HHMM), e.g. 0000, 0300, 0600, etc.
    year_str = str(year)  # Convert year to string
    month_str = f"{month:02d}"  # Format month as a two-digit string
    dom_str = f"{dom:02d}"  # Format day of month as a two-digit string

    ## the data is stored in a directory structure like this:
    # https://www.ncei.noaa.gov/data/international-satellite-cloud-climate-project-isccp-h-series-data/access/isscp-basic/hgg/YEARMM/files


    os.makedirs(local_dir, exist_ok=True)  # Create the local directory if it doesn't exist


    local_filename =  locator.format_filename(product_uppercase[product],version_uppercase,year=year,doy = doy,time = time_str) ## This is the expected local filename format for ISCCP data,
    local_path = os.path.join(local_dir, local_filename)

    
    remote_dir = f"{base_url}/{product}/{version}/{year_str}{month_str}/"
    remote_filename = f'{product_uppercase[product]}.{version_uppercase}.v01r00.GLOBAL.{year}.{month_str}.{dom_str}.{time_str}.GPC.10KM.CS00.EA1.00.nc'  # Construct the remote filename


    if (not os.path.exists(local_path)) or force_redownload:
        remote_url = f"{remote_dir}{remote_filename}"
        logging.info(f"Downloading from {remote_url} to {local_path}")
        get_file_from_url(remote_url, local_path)
        return local_path  # Return the local path of the downloaded file

    else:
        # If the file already exists, skip the download
        logging.info(f"File {local_filename} already exists in {local_dir}. Skipping download.")
        return local_path
    


    




    

    
