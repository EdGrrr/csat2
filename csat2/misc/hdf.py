from pyhdf import SD, HDF, VS # VS is not used in this code, but not importing it generates issues elsewhere
import numpy as np
import copy


def read_hdf4(filename, names=None, vdata=None, fill_missing=True, datadict=None):
    '''Reads data from a HDF4 file into a dictionary.'''

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
                    fillvalue = attributes['_FillValue']  # Missing data.
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
                data = vd.read(nRec=vd.inquire()[0])
                for x in range(0, len(data)):
                    data[x] = data[x][0]
                data = np.array(data)
                if fill_missing:
                    try:  # Deal with missing data
                        missingval = vd.attrinfo()['missing'][2]
                        wfillmask = 0
                        wfillmask = np.where(data == missingval, np.nan, 1)
                        data = data * wfillmask
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
