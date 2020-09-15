from __future__ import print_function
from glob import glob
import os
import logging
import socket
import yaml
import pkgutil
import fnmatch
from csat2.misc.time import doy_to_date


class FileLocator:
    '''csat2 locator class. A simple database type object that returns available
    files on the system mathing a given pattern.

    Can deal with multiple potential locations for files through the specification
    of multiple patterns for a given (origin, product) combination.

    The file format is (indentation is not important):

    M macro_name=path_to_folder
    Include a macro (only at the start of the file)
    Macros are declared before they are used.  They are replaced as the paths
    are read in (only once).  Note that 'product' is automatically a macro

    [source] (e.g. MODIS, CERES, DARDAR)
    - [product1|product2] (products following a specific pattern, e.g. MASK, CLOUD,
       split multiple with pipes - the actual product is automatically replaced where {product} is used)

    e.g.
    [CALIOP]
    -[05kmALay|05kmCLay|01kmCLay|333mCLay|VFM]
      {csat_folder}/CALIOP/{product}/{ver}/{year}/CAL*{product}*{year}-{mon:0>2}-{day:0>2}T*.hdf


    Path syntax is the same as the python string formatting syntax (what a coincidence...)
    Use {variable} to make a substitution - can be a macro or a variable from the readin function
    {doy:0>3} left pads the variable doy with 0s to a width of three characters
    * can be used as a wildcard (e.g. for MODIS files that include the processing date
        The first file is often the one that is returned (depending on how locator.search is called)
    Note that * can also be passed as an argument to locator.search, e.g. if you want to return all the
     MODIS granules for a specific day
'''

    def __init__(self, searchfile=None, data=''):
        self.search_paths = {}
        self.macros = {}
        if searchfile:
            with open(searchfile) as f:
                logging.debug('Loading {}'.format(searchfile))
                searchfiledata = f.readlines()
            self.load_search_paths(searchfiledata)
        elif data:
            self.load_search_paths(data)

    def search(self, origin, product, exit_first=True, *args, **kwargs):
        '''Searches the available patterns for the given origin and product.
        Note that a pattern can contain '*' to return all available paths

        exit_first - return the results from the first pattern with results'''
        patterns = self.get_patterns(origin, product, **kwargs)
        filenames = []
        for p in patterns:
            pat = p.format(**kwargs)
            logging.debug('Pattern: '+pat)
            filenames.extend(glob(pat))
            if exit_first and (len(filenames) > 0):
                return filenames
        return filenames

    def get_patterns(self, origin, product, **kwargs):
        sd = kwargs
        if ('year' in sd.keys()) and ('doy' in sd.keys()):
            # Get the month and day if only the year and doy are provided
            _, sd['mon'], sd['day'] = doy_to_date(sd['year'], sd['doy'])
        try:
            patterns = self.search_paths[origin][product]
        except:
            raise
        return patterns

    def paths(self):
        return self.search_paths

    def load_search_paths(self, searchfiledata):
        search_paths = {}
        macro = {}
        origin = None
        products = None

        for line in searchfiledata:
            line = line.strip()
            if line == '' or line[0] == '#':
                # Comment
                continue
            elif line[0] == 'M':
                temp = line.split('=')
                macro[temp[0][2:]] = temp[1].strip()
            elif line[0] == '[':
                # Origin (MODIS, TRMM etc.)
                origin = line[1:-1]
                search_paths[origin] = {}
            elif line[0] == '-':
                # Products
                products = line[2:-1].split('|')
                for product in products:
                    search_paths[origin][product] = []
            else:
                for product in products:
                    macro['product'] = product
                    (search_paths[origin][product]
                     .append(partial_format(line.strip(), macro)))
        self.search_paths = search_paths
        del(macro['product'])
        self.macro = macro

    def get_folder(self, origin, product, return_primary=True, *args, **kwargs):
        '''Get the folder locations from the search paths dictionary.
        return_primary - return the first matching path (assumed main directory)'''
        patterns = self.get_patterns(origin, product, **kwargs)
        patterns = [os.path.dirname(p).format(**kwargs)
                    for p in patterns]
        if return_primary:
            return patterns[0]
        else:
            return patterns


def partial_format(istr, format_dict):
    '''Operates in the same way as the starndard python format for string
    interpolation, but doesn't fail if a key is missing'''
    s = 0
    outstr = ''
    while True:
        try:
            i = istr.index('{', s)
            j = istr.index('}', i)
            outstr += istr[s:i]
            key = istr[i+1:j]
            if key in format_dict.keys():
                outstr += format_dict[key]
            else:
                outstr += ('{'+key+'}')
            s = j+1
        except ValueError:
            outstr += istr[s:]
            return outstr

def configloader(user_location=None, default_fallback=True):
    '''Searches for and loads a search path configuration
    user_location - specified location for the config file
    default_fallback - use the package supplied defaults'''
    _hostname = socket.getfqdn()
    _package_directory = os.path.dirname(os.path.abspath(__file__))

    if user_location:
        config = yaml.safe_load(open(user_location))

    if (not user_location) or default_fallback:
        logging.info('Appending default config')
        data = pkgutil.get_data(__name__, "config/config.cfg")
        newconfig = yaml.safe_load(data)
        if not config:
            config = newconfig
        else:
            for name in newconfig['machines'].keys():
                config['machines'][name] = newconfig['machines'][name]

    if (not user_location) and (not default_fallback):
        raise ValueError('Must supply user_location or set default feedback as true')

    logging.debug('machine_config{}'.format(config))
    
    # Match the current machine with the config file
    for machine in config['machines'].keys():
        if fnmatch.fnmatch(_hostname, machine):
            logging.info('Current machine is: {}'.format(machine))
            spfile = config['machines'][machine]
            mfile = os.path.join(os.path.dirname(user_location), spfile)
            # First check the override directory and use in preference
            if os.path.exists(mfile):
                try:
                    locator = FileLocator(mfile)
                except:
                    raise IOError(mfile + ' is not a valid config file')
            elif default_fallback:
                mdata = pkgutil.get_data(__name__, "config/{}".format(spfile)).decode('ascii')
                mdata = mdata.split('\n')
                print(mdata)
                try:
                    locator = FileLocator(data=mdata)
                except:
                    # This means that the correct source file has been found
                    # but there was no suitable instrument entry
                    #
                    # Assumes that the distributed files are valid!
                    # Make sure to test!
                    raise IOError('No suitable default file')
            else:
                raise IOError('No matching machine file in {}'.format(user_location))
            break # Break on a match
    else:
        logging.warning('No matching machine found, empty locator created')
        locator = FileLocator()

    return locator
