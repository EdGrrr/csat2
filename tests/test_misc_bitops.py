from csat2.misc import bitops
import numpy as np
import pytest


def test_bitopsSelectPattern():
    assert bitops.bit_select_pattern(np.array([9]), '1000') == 1
    assert bitops.bit_select_pattern(np.array([9]), '0100') == 0
    assert bitops.bit_select_pattern(np.array([9]), '0010') == 0
    assert bitops.bit_select_pattern(np.array([9]), '0001') == 1
    assert bitops.bit_select_pattern(np.array([13]), '1100') == 3


def test_bitopsSelectLocations():
    assert bitops.bit_select_locations(np.array([9]), (0, 1)) == 1
    assert bitops.bit_select_locations(np.array([9]), (1, 2)) == 0
    assert bitops.bit_select_locations(np.array([9]), (2, 3)) == 2
    assert bitops.bit_select_locations(np.array([13]), (2, 3)) == 3


def test_bitopsSelectLocationsFailure():
    with pytest.raises(ValueError):
        bitops.bit_select_locations(np.array([13]), (1, 3))

    with pytest.raises(ValueError):
        bitops.bit_select_locations(np.array([13]), (3, 2))

        
