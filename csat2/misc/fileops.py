import struct
import numpy as np
from netCDF4 import Dataset
import datetime
import pandas as pd

####################
# NetCDF functions #
####################


def nc_dump(filename, data, dnames=None, vname=None, format=None):
    '''Dumps and array into a netcdf file to store for later'''
    nc4flgs = ['NETCDF4', 'nc4', 'NC4']
    if format in nc4flgs:
        ncdf = Dataset(filename, 'w', format='NETCDF4')
    else:
        ncdf = Dataset(filename, 'w')
    ndims = len(data.shape)

    if not dnames:  # Use numbers as dimension names if none are given
        dnames = tuple(str(i) for i in range(ndims))

    for i, sdim in zip(dnames, data.shape):
        ncdf.createDimension(str(i), sdim)

    # Use 'var' if no variable name is given
    if not vname:
        vname = 'var'
    if not type(dnames) == type((1, 2)):
        dnames = tuple(dnames)

    if format in nc4flgs:
        Var = ncdf.createVariable(vname, 'f', dnames, zlib=True)
    else:
        Var = ncdf.createVariable(vname, 'f', dnames)
    Var[:] = data.astype('float')
    ncdf.close()
    return


def nc4_dump(*args, **kwargs):
    '''Dumps and array into a netcdf file to store for later'''
    return nc_dump(*args, format='NETCDF4', **kwargs)


def nc_load(filename, vname=None):
    ncdf = Dataset(filename, 'r', format='NETCDF4')
    if vname:
        data = ncdf.variables[vname][:]
    else:
        data = {}
        for name in ncdf.variables.keys():
            data[name] = ncdf.variables[name][:]
    ncdf.close()
    return data


def nc4_dump_mv(filename, data, dnames=None):
    '''Dumps and array into a netcdf file to store for later.  
    Takes a dictionary of arrays of the same partial shape where the key is the variable name.
    i.e. (100,20,35) is the same as (100,20).'''
    ncdf = Dataset(filename, 'w', format='NETCDF4')
    names = list(data.keys())
    ndims = max([len(data[name].shape) for name in names])

    if not dnames:  # Use numbers as dimension names if none are given
        dnames = tuple(str(i) for i in range(ndims))

    fullshape = data[[name for name in names
                      if len(data[name].shape) == ndims][0]].shape
    for i, sdim in zip(dnames, fullshape):
        ncdf.createDimension(str(i), int(sdim))

    # Use 'var' if no variable name is given
    #if not vname: vname = 'var'
    if not type(dnames) == type((1, 2)):
        dnames = tuple(dnames)

    for vname in names:
        Var = ncdf.createVariable(
            vname, 'f', dnames[:len(data[vname].shape)], zlib=True)
        Var[:] = data[vname].astype('float')
    ncdf.close()
    return

#####################
# Fortran functions #
#####################


def fortran_read(file_obj, ftype):
    '''Read a FORTRAN formatted datafile'''
    length = struct.unpack('>i', file_obj.read(4))[0]
    outdata = file_obj.read(length)
    upstring = '>'+str(length/struct.calcsize(ftype))+ftype
    outdata = struct.unpack(upstring, outdata)
    endlength = struct.unpack('>i', file_obj.read(4))[0]
    if (length == endlength):
        return length, outdata
    else:
        raise IOError('Malformed data record')


def fortran_write(file_obj, data, ftype):
    '''Write a FORTRAN formatted datafile'''
    if (type(data) == type(1)) or (type(data) == type(1.)):
        data = np.array([data]).astype(ftype)
    if (ftype == 'c') and (type(data) == type('string')):
        tdata = []
        tdata[:0] = data
        data = np.array(tdata)
    else:
        data = np.array(data).astype(ftype)
    length = len(data.ravel())
    number = length*struct.calcsize(ftype)
    file_obj.write(struct.pack('>i', number))
    # Write the data
    upstring = '>'+str(int(length))+ftype
    file_obj.write(struct.pack(upstring, *list(data.ravel())))
    # Put the ending length marker
    file_obj.write(struct.pack('>i', number))
    return


def chunks(l, n):
    '''Split a list 'l' into chunks of length 'n' '''
    for i in range(0, len(l), n):
        yield l[i:i+n]

######################
# Other file formats #
######################

def readin_icartt(filename):
    '''ICARTT data format (Langley flight data)
    based on the format description here 
    https://www-air.larc.nasa.gov/missions/etc/IcarttDataFormat.htm'''
    with open(filename, 'rb') as f:
        header = [f.readline().strip().decode('ascii')]
        header_count, _ = header[0].split(',')
        for i in range(int(header_count)-2):
            header.append(f.readline().strip().decode('ascii'))

        colnames = f.readline()
        colnames = colnames.strip().decode('ascii').split(',')

        data = []
        while True:
            datarow = f.readline().strip().decode('ascii').split(',')
            if len(datarow)>1:
                data.append(datarow)
            else:
                break

        data = np.array(data, dtype=np.float)

        starttime = datetime.datetime.strptime(filename.split('_')[-2], '%Y%m%d')
        
        data = pd.DataFrame(data, columns=colnames)
        for name in data.keys():
            if 'UTC' in name:
                data[name] = pd.to_datetime(data[name], unit='s', origin=starttime)
        return data
