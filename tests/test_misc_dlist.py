from csat2.misc import dlist
import unittest
import numpy as np
import copy
import xarray as xr


class TestDirectoryKeys(unittest.TestCase):
    def test_dictkeysappend_axis_m1(self):
        dict1 = {}
        dict2 = {"test": np.array([[1, 2]])}
        dlist.dictkeysappend(dict1, dict2)
        self.assertEqual(dict1, dict2)
        dlist.dictkeysappend(dict1, dict2)
        self.assertEqual(len(dict1.keys()), 1)
        self.assertSequenceEqual(dict1["test"].shape, (1, 4))
        self.assertEqual(list(dict2.keys()), ["test"])

    def test_dictkeysappend_axis_zero(self):
        dict1 = {}
        dict2 = {"test": np.array([[1, 2]])}
        dlist.dictkeysappend(dict1, dict2, axis=0)
        self.assertEqual(dict1, dict2)
        dlist.dictkeysappend(dict1, dict2, axis=0)
        self.assertEqual(len(dict1.keys()), 1)
        self.assertSequenceEqual(dict1["test"].shape, (2, 2))
        self.assertEqual(list(dict2.keys()), ["test"])

    def test_dictkeysappend_axis_longest(self):
        dict1 = {}
        dict2 = {"test": np.array([[1, 2]])}
        dlist.dictkeysappend(dict1, dict2, axis=0)
        self.assertEqual(dict1, dict2)
        dlist.dictkeysappend(dict1, dict2, axis="longest")
        self.assertEqual(len(dict1.keys()), 1)
        self.assertSequenceEqual(dict1["test"].shape, (1, 4))
        self.assertEqual(list(dict2.keys()), ["test"])


class TestDirectoryList(unittest.TestCase):
    def test_ldictkeysappend(self):
        dict1 = {}
        dict2 = {"test": np.array([[1, 2]])}
        dlist.l_dictkeysappend(dict1, dict2)
        self.assertEqual(dict1.keys(), dict2.keys())
        self.assertEqual(len(dict1["test"]), 1)

        dlist.l_dictkeysappend(dict1, dict2)
        self.assertEqual(dict1.keys(), dict2.keys())
        self.assertEqual(len(dict1["test"]), 2)

        dtemp = copy.deepcopy(dict1)
        dlist.l_toarray(dtemp, axis=-1)
        self.assertEqual(dtemp.keys(), dict2.keys())
        self.assertEqual(dtemp["test"].shape, (1, 4))

        dtemp = copy.deepcopy(dict1)
        dlist.l_toarray(dtemp)
        self.assertEqual(dtemp.keys(), dict2.keys())
        self.assertEqual(dtemp["test"].shape, (2, 2))

        dlist.l_toarray(dict1, axis=0)
        self.assertEqual(dict1.keys(), dict2.keys())
        self.assertEqual(dict1["test"].shape, (2, 2))


class TestDirectoryListXarray(unittest.TestCase):
    def test_ldictkeysappend_xarray(self):
        dict1 = {}
        dict2 = {
            "test": xr.DataArray(
                np.array([[1, 2]]),
                coords={"x": np.array([0]), "y": np.array([0, 1])},
                dims=("x", "y"),
            )
        }
        dlist.l_dictkeysappend(dict1, dict2)
        self.assertEqual(dict1.keys(), dict2.keys())
        self.assertEqual(len(dict1["test"]), 1)

        dlist.l_dictkeysappend(dict1, dict2)
        self.assertEqual(dict1.keys(), dict2.keys())
        self.assertEqual(len(dict1["test"]), 2)

        # Unspecified axis
        dtemp = copy.deepcopy(dict1)
        dlist.l_toarray(dtemp)
        self.assertEqual(dtemp.keys(), dict2.keys())
        self.assertEqual(dtemp["test"].shape, (2, 2))

        # Specified axis
        dtemp = copy.deepcopy(dict1)
        dlist.l_toarray(dtemp, axis=0)
        self.assertEqual(dtemp.keys(), dict2.keys())
        self.assertEqual(dtemp["test"].shape, (2, 2))

        dtemp = copy.deepcopy(dict1)
        dlist.l_toarray(dtemp, axis=-1)
        self.assertEqual(dtemp.keys(), dict2.keys())
        self.assertEqual(dtemp["test"].shape, (1, 4))

        # Named axis
        dlist.l_toarray(dict1, axis="x")
        self.assertEqual(dict1.keys(), dict2.keys())
        self.assertEqual(dict1["test"].shape, (2, 2))
