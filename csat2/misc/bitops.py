'''Bit operation functions'''

import numpy as np


def bit_select_pattern(data, pattern):
    '''Takes an integer array and then removes only the data specified 
    in the pattern eg' 01100000', shifting it so the lowest required bit is now 0/1'''
    outdata = np.bitwise_and(data.astype('int'), int(pattern, 2))
    shift_cnt = len(pattern) - pattern.rfind('1') - 1
    return np.right_shift(outdata, shift_cnt)


def bit_select_locations(data, locations):
    '''Select bit locations, rather than providing the pattern.
    Locations must be a monotonically increasing sequence with no gaps
    e.g. locations=(3,4,5) or locations=3'''
    if isinstance(locations, int):
        return np.mod(np.right_shift(data, locations), 2)
    elif len(locations) == 1:
        return np.mod(np.right_shift(data, locations[0]), 2)
    elif np.all(np.diff(locations)==1):
        shift_cnt = min(locations)
        return np.mod(np.right_shift(data, shift_cnt), 2**len(locations))
    else:
        raise ValueError('Locations must be monotonically increasing with no gaps')
