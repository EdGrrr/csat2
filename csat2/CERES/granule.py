import os
from csat2 import locator
from csat2 import misc
from csat2.CERES.download import download_files



class Granule:
    '''An object for accessing and managing a specific CERES granule
    Note that this is only currently set up to deal with the CERES SYN - hourly product, and may have to be substancially editted if
    you want the data at different temporal frequencies.'''
    def __init__(self,year:int, doy :int, time:int):
        self.year = year
        _,month,dom = misc.time.doy_to_date(year, doy)  # Get month and day of month (dom) from year and day of year
        self.month = month
        self.dom = dom
        self.doy = doy
        self.time = time ## this is the UTC time in hours (0-23)
        

        
    def get_fileloc(self):
        ''' Return the local directory, note that this does not require that the file actually exists there.'''
        
        fileloc = locator.format_filename('CERES','SYN-hourly',year = self.year, month = self.month, dom  = self.dom)

        return fileloc
  

    def check(self):
        ''' Check if the file exists locally.'''
        fileloc = self.get_fileloc()
        if os.path.exists(fileloc):
            return True
        else:
            return False
        
    def download(self, force_redownload: bool = False):
            '''Download the file locally'''
            fileloc = self.get_fileloc()
            if not self.check() or force_redownload:
                download_files(
                    year=self.year,
                    month=self.month,  # Now passes month instead of doy
                    dom=self.dom,
                    local_path=fileloc
                )
                return f'Downloaded file to {fileloc}'
            else:
                return f'File already exists at {fileloc}, not redownloading.'
        
    # def get_variable(self,varname:str):
    #     ''' Get a specific variable from the granule. Will download the file if it does not already exist locally.
    #    Note that the files are saved on the server in .hdf format'''

        
    #     # Ensure the file is downloaded
    #     self.download()
        
    #     fileloc = self.get_fileloc()
        
    #     # Open the netCDF file and extract the variable
    #     with nc.Dataset(fileloc, 'r') as dataset:
    #         variable_data = dataset.variables[varname][:]
        
    #     return variable_data