'''
Functions for reading in lists of data into dictionaries - eg. satellite data
Either:
    Use dictkeysappend - slow
    Use l_dictkeysappend and l_toarray after finishing the readin loop - faster
'''
import numpy as np
import logging
import xarray as xr
log = logging.getLogger(__name__)


def dictkeysappend(dict_main, dict_addition, axis=-1):
    '''Adds data from keys in dict_addition to dict_main, if that keys doesn't exist, it is created.
    dict1 is modified inplace'''
    if axis != 'longest':
        naxis = axis
    for key in dict_addition.keys():
        try:
            if axis == 'longest':
                naxis = np.argmax(dict_main[key].shape)
            dict_main[key] = np.concatenate((dict_main[key], dict_addition[key]), axis=naxis)
        except KeyError:
            dict_main[key] = dict_addition[key]


def l_dictkeysappend(dict_main, dict_addition):
    '''Add the data in dict_addition to the end of the lists in dict_main'''
    for key in dict_addition.keys():
        if key in dict_main.keys():
            dict_main[key].extend([dict_addition[key]])
        else:
            dict_main[key] = [dict_addition[key]]


def l_toarray(dict1, axis=None):
    '''Converts the lists in dict1 to arrays
    Modifies dict1 inplace'''
    names = list(dict1.keys())
    if isinstance(dict1[names[0]][0], (xr.DataArray)):
        log.debug('xarray merging')
        if axis is None:
            # If axis is not specified, use the zeroth axis
            axis = 0
        for key in names:
            # Should allow specification of a dimension names
            # e.g. 'time'
            if isinstance(dict1[key], (list)):
                if axis in dict1[key][0].dims:
                    dict1[key] = xr.concat(
                        dict1[key],
                        dim=axis)
                else:
                    dict1[key] = xr.concat(
                        dict1[key],
                        dim=dict1[names[0]][0].dims[axis])
    else:
        log.debug('numpy merging')
        if not axis:
            axis = 0
        for key in names:
            dict1[key] = np.concatenate(dict1[key], axis=axis)
