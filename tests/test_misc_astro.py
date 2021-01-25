from csat2.misc import astro, time
import unittest
import numpy as np
from datetime import datetime


class TestSzaCalc(unittest.TestCase):
    def setUp(self):
        self.data = [
            # Test values from NOAA calcualtor https://www.esrl.noaa.gov/gmd/grad/solcalc/
            # (Datetime, lon, lat), (sunrise (hour, min), sza)
            [(datetime(2020, 6, 20, 12), 0, 50), ((3, 56), 90-63.44)],
            [(datetime(2020, 1, 1, 12), 0, 50), ((8, 1), 90-17.03)]]

    def test_SunriseTime(self):
        for dat in self.data:
            stime = astro.sunrise_time(
                dat[0][2], time.datetime_to_ydh(dat[0][0])[1])
            assert int(stime//1) == dat[1][0][0]
            assert int(np.round(60*(stime % 1))) == dat[1][0][1]

    def test_SZA(self):
        for dat in self.data:
            sza = astro.solar_zenith_angle_time(
                dat[0][2], time.datetime_to_ydh(dat[0][0])[1],
                time.datetime_to_lst(dat[0][0], dat[0][1]))
            assert np.abs(sza-dat[1][1]) < 0.25
