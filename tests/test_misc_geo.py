from csat2.misc import geo
import unittest
import numpy as np


class TestGeoFunctions(unittest.TestCase):
    def setUp(self):
        self.lat = np.arange(-70, 70, 10)
        self.lon = np.zeros(self.lat.shape)

    def test_haversine_array(self):
        # Make sure that the distances along a meridian are equal
        output = geo.haversine(
            self.lat[:-1], self.lon[:-1], self.lat[1:], self.lon[1:])
        self.assertTrue(np.mean(np.diff(output)) < 10e-10)

    def test_haversine_absolute(self):
        self.assertEqual(geo.haversine(
            -86.67, 36.12, -118.40, 33.94), 2887.2599506071106)

    def test_bearing(self):
        assert np.isclose(geo.bearing(0, 0, 180, 1), 0)
        assert np.isclose(geo.bearing(180, 1, 0, 0), 0)
        assert np.isclose(geo.bearing(0, 0, 90, 0), 90)

    def test_coordinate_rotate(self):
        # Not clear what the ideal behaviour is here, but this makes sure
        # it doesn't crash...
        assert np.all(np.isclose(geo.coordinate_rotate(0, 0, 0, 0), (0, -90)))
        assert np.all(np.isclose(geo.coordinate_rotate(0, 0, 0, 90), (0, -90)))
        assert np.all(np.isclose(geo.coordinate_rotate(0, 0, 90, 0), (-90, 0)))
