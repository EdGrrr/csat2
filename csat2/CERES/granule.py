class Granule
    '''An object for accessing and managing a specific CERES granule
    Note that this is only currently set up to deal with the CERES SYN - hourly product, and may have to be substancially editted if
    you want the data at different temporal frequencies.'''
    def __init__(self,year:int, doy :int, time:int):
        self.year = year
        self.doy = doy
        self.time = time ## this is the UTC time in hours (0-23)

        
    def get_fileloc(self):
        ''' Return the local directory'''
  