from pyhdf import (
    SD,
    HDF,
    VS,
)  # VS is not used in this code, but not importing it generates issues elsewhere
import numpy as np
import copy
import logging
import subprocess
import netCDF4

log = logging.getLogger(__name__)


def check_netcdf_hdf4():
    try:
        result = subprocess.run(["nc-config", "--has-hdf4"], capture_output=True, text=True, check=True)
        hdf4_support = result.stdout.strip()
        if hdf4_support == 'yes':
            return True
        else:
            return False
    except FileNotFoundError:
        return False
    except subprocess.CalledProcessError as e:
        return False


# Check if the installed version of netCDF4 can open hdf4 files
if  check_netcdf_hdf4():
    Dataset = netCDF4.Dataset
else:
    class Variable:
        def __init__(self, name, datafile):
            self.name = name
            self.sds = datafile.select(name)
            self.attrs = self.sds.attributes()
            self.auto_scale = True
            self.auto_mask = True
            self.dimensions = datafile.datasets()[name][0]
            self.shape = datafile.datasets()[name][1]
            try:
                self._Fillvalue = self.attrs["_FillValue"]
            except KeyError:
                pass

        def set_auto_scale(self, val):
            self.auto_scale = val

        def set_auto_mask(self, val):
            self.auto_mask = val

        def ncattrs(self):
            return self.attrs.keys()

        def getncattr(self, name):
            return self.attrs[name]

        def __getitem__(self, index):
            data = self.sds[index]  # Allow indexing like a NumPy array
            if self.auto_mask:
                data = np.ma.array(data, mask=(data == self._FillValue))
            if self.auto_scale:
                try:
                    data = data * self.scale_factor + self.add_offset
                except AttributeError:
                    pass
            return data

        def __setitem__(self, index, value):
            raise AttributeError("Cannot edit data")

        def __getattr__(self, name):
            try:
                return self.attrs[name]
            except:
                raise AttributeError("Not a valid atribute")

    class Variables:
        def __init__(self, datafile):
            self.datafile = datafile
            self.names = list(self.datafile.datasets().keys())

        def __getitem__(self, name):
            if name in self.names:
                return Variable(name, self.datafile)
            else:
                raise KeyError('Not a valid variable')

        def keys(self):
            return self.datafile.datasets().keys()

    class Dataset:
        def __init__(self, filename, mode="r"):
            if mode != "r":
                raise NotImplementedError('Read only for hdf files')
            self.datafile = SD.SD(filename)
            self.variables = Variables(self.datafile)

        def __enter__(self):
            return self
            
        # def __del__(self):
        #     self.datafile.end()

        def __exit__(self, exc_type, exc_value, traceback):
            if self.datafile is not None:
                self.datafile.end()


def read_hdf4(filename, names=None, vdata=None, fill_missing=True, datadict=None):
    """Reads data from a HDF4 file into a dictionary."""

    # Is set to true if we are reading in all variables
    nonameflag = False

    if not vdata:
        datafile = SD.SD(filename)
        if not names:
            names = datafile.datasets()
        if not datadict:
            datadict = {}

        for name in names:
            sds = datafile.select(name)
            data = sds.get()
            attributes = sds.attributes()
            if fill_missing:
                try:
                    fillvalue = attributes["_FillValue"]  # Missing data.
                    data = np.where(data == fillvalue, np.nan, data)
                except KeyError:
                    pass
            datadict[name] = data

    else:
        # HDF vdata format
        datafile = HDF.HDF(filename)
        vs = datafile.vstart()
        if not names:
            names = vs.vdatainfo()
            names = list(zip(*names))
            names = names[0]
            nonameflag = True

        datadict = {}
        failednames = []
        for name in names:
            try:
                vd = vs.attach(name)
                vdnames = [v[0] for v in vd.fieldinfo()]
                vddata = vd.read(vd._nrecs)
                if len(vdnames) == 1:  # Cloudsat data
                    data = np.array(vddata).squeeze()
                elif len(vdnames) > 1:  # CALIPSO metadata
                    data = {vdnames[i]: vddata[0][i] for i in range(len(vdnames))}
                if fill_missing:
                    try:  # Deal with missing data
                        missingval = vd.attrinfo()["missing"][2]
                        data = np.where(data == missingval, np.nan, data)
                    except KeyError:
                        pass
                datadict[name] = data
                vd.detach()
            except:
                failednames.append(name)
        vs.end()
        datafile.close()
        if nonameflag:
            datadict = read_hdf4(filename, names=None, datadict=datadict)
        elif len(failednames) > 0:
            datadict = read_hdf4(filename, failednames, datadict=datadict)

    return datadict


def read_hdf4_metadata(filename, names, vdata=True):
    "Retrieves long name and units from hdf4 file, for a single variable."
    output_metadata = {}
    vnames = copy.copy(names)
    if vdata:
        datafile = HDF.HDF(filename)
        vs = datafile.vstart()
        for name in vnames:
            try:
                vd = vs.attach(name)
                attributes = vd.attrinfo()
                var_attr = {}
                for attribute in attributes:
                    var_attr[attribute] = attributes.get(attribute)
                output_metadata[name] = var_attr
                names.remove(name)
                vd.detatch()
            except:
                pass
        vs.end()
        datafile.close()

    datafile = SD.SD(filename)
    for name in names:
        sds = datafile.select(name)
        attributes = sds.attributes()
        var_attr = {}
        for attribute in attributes:
            var_attr[attribute] = attributes.get(attribute)
        output_metadata[name] = var_attr
    return output_metadata


def read_hdf5(filename, sds=None):
    """
    Read specified datasets from an EarthCARE HDF5 file into a flat dict.

    Parameters
    ----------
    filename : str
        Path to the .h5 file
    sds : list of str, optional
        Full dataset paths to include (e.g. 'ScienceData/latitude').
        If None, loads all datasets.

    Returns
    -------
    dict : mapping dataset name -> numpy array
    """
    try:
        data = {}
        with h5py.File(filename, 'r') as h5:
            def visitor(name, obj):
                if isinstance(obj, h5py.Dataset):
                    if sds is None or name in sds:
                        data[name] = obj[()]
            h5.visititems(visitor)
        if not data:
            raise IOError(f"No datasets found in {filename}")
        return data
    except Exception as e:
        raise IOError(f"Failed to read HDF5 file {filename}: {e}")
