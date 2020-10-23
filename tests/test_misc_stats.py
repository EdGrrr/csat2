import unittest
from csat2.misc import stats
import scipy.stats
import numpy as np


class TestNans(unittest.TestCase):
    def test_nanfunctions(self):
        testx = np.array([1, 0, 2, np.nan, 0.5])

        self.assertEqual(stats.nanmean(testx), np.mean([0, 1, 2, 0.5]))
        self.assertEqual(stats.nanstd(testx), np.std([0, 1, 2, 0.5]))
        self.assertEqual(stats.nanmax(testx), 2)
        self.assertEqual(stats.nanmin(testx), 0)
        self.assertEqual(stats.zero_nans(testx)[-2], 0)
        self.assertEqual(stats.zero_nans(testx)[0], 1)
        self.assertEqual(stats.lin_av(testx[:3])[0], 0.5)
        self.assertEqual(stats.lin_av(testx[:3])[1], 1)

        testy = np.array([np.nan, 0, 2, 1, 1])
        self.assertEqual(stats.nanlinregress(testx, testy),
                         scipy.stats.linregress(testx[[1, 2, 4]], testy[[1, 2, 4]]))

        assert stats.nanmask(testx).mask is not None


class TestLatWeights(unittest.TestCase):
    def setUp(self):
        self.data = np.fromfunction(lambda x, y: x, (10, 1))
        self.data2 = np.array([[0], [0], [1]])

    def test_latweightedav(self):
        assert np.isclose(stats.lat_weighted_av(self.data, [-90, 90]),  4.5)
        assert np.isclose(stats.lat_weighted_av(self.data % 5, [-90, 90]),  2)
        assert np.isclose(stats.lat_weighted_av(self.data % 5, [-60, 60]),  2)
        assert np.isclose(stats.lat_weighted_av(self.data % 5, [-30, 30]),  2)

        assert stats.lat_weighted_av(
            self.data, (-90, 0)) > stats.lat_weighted_av(self.data, (0, 90))
        assert np.isclose(stats.lat_weighted_av(self.data2, (-90, 0)), 0.5)
        assert np.isclose(stats.lat_weighted_av(
            self.data2, (0, 90)), 1-np.sqrt(3)/2)
