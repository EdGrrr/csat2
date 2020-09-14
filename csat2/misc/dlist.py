'''
Functions for reading in lists of data into dictionaries - eg. satellite data
Either:
    Use dictkeysappend - slow
    Use l_dictkeysappend and l_toarray after finishing the readin loop - faster
'''
import numpy as np


def dictkeysappend(dict1, dict2, axis=-1):
    '''Adds data from keys in dict 2 to dict 1, if that keys doesn't exist, it is created.
    dict1 is modified inplace'''
    if axis != 'longest':
        naxis = axis
    for key in dict2.keys():
        try:
            if axis == 'longest':
                naxis = np.argmax(dict1[key].shape)
            dict1[key] = np.concatenate((dict1[key], dict2[key]), axis=naxis)
        except KeyError:
            dict1[key] = dict2[key]


def l_dictkeysappend(dict1, dict2):
    '''Add the data in dict2 to the end of the lists in dict1'''
    for key in dict2.keys():
        if key in dict1.keys():
            dict1[key].extend([dict2[key]])
        else:
            dict1[key] = [dict2[key]]


def l_toarray(dict1, axis=0):
    '''Converts the lists in dict1 to arrays
    Modifies dict1 inplace'''
    for key in dict1.keys():
        dict1[key] = np.concatenate(dict1[key], axis=axis)
