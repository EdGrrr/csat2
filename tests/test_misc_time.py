from csat2.misc import time
import unittest
import numpy as np
from datetime import datetime


class TestLSTConversion(unittest.TestCase):
    def test_dataIsGMT(self):
        """
        Data is UTC\n
        Uses UTC value as the data, so should return the UTC \n
        value of lst_time at each longitude"""
        lst_time = 13.5
        gmt_times = np.arange(0, 25, 3)
        lons = 30
        longitudes = np.arange(-180, 180, 360 / lons)
        data = np.zeros((2, lons, len(gmt_times)))
        for i in range(len(gmt_times)):
            data[:, :, i] = gmt_times[i]
        data = data.astype("float")
        output = time.toLocalSolarTime(lst_time, gmt_times, longitudes, data)
        self.assertTrue(
            (time.utc_to_lst(output, longitudes) == lst_time).all(),
            "GMT - LST conversion failed",
        )

    def test_dataIsGMTRev(self):
        """
        Data is UTC (reversed longitudes)\n
        As test_data_is_rev but with longitudes reversed"""
        lst_time = 13.5
        gmt_times = np.arange(0, 25, 3)
        lons = 30
        longitudes = np.arange(180, -180, -360 / lons)
        data = np.zeros((2, lons, len(gmt_times)))
        for i in range(len(gmt_times)):
            data[:, :, i] = gmt_times[i]
        data = data.astype("float")
        output = time.toLocalSolarTime(lst_time, gmt_times, longitudes, data)
        self.assertTrue(
            (time.utc_to_lst(output, longitudes) == lst_time).all(),
            "GMT - LST reverse conversion failed",
        )

    def test_dataIsLST(self):
        """
        Data is LST\n
        Runs a test where the data to be lst'd is the LST at each longitude and UTC value.  Should return the LST, with some possible edge cases near the dateline and greenwich which are not completely sorted"""
        lst_time = 13.5
        gmt_times = np.arange(0, 25, 3)
        lons = 30
        longitudes = np.arange(-180, 180, 360 / lons)
        data = np.zeros((2, lons, len(gmt_times)))
        for i in range(len(gmt_times)):
            data[:, :, i] = time.utc_to_lst(gmt_times[i], longitudes)[
                np.newaxis, :
            ].repeat(2, axis=0)
        data = data.astype("float")
        output = time.toLocalSolarTime(lst_time, gmt_times, longitudes, data)
        total = np.isfinite(output).ravel().sum()
        equal = float((abs(output - lst_time) < 0.001).ravel().sum())
        self.assertTrue(((equal / total) > 0.9), "LST as data conversion failed")

    def test_raiseErrors(self):
        lst_time = 13.5
        gmt_times = np.array([0, 0])
        lons = 30
        longitudes = np.arange(-180, 180, 360 / lons)
        data = np.zeros((2, lons, len(gmt_times)))
        self.assertRaises(
            ValueError, time.toLocalSolarTime, lst_time, gmt_times, longitudes, data
        )
        gmt_times = np.array([0, 12])
        self.assertRaises(
            ValueError, time.toLocalSolarTime, lst_time, gmt_times, longitudes, data
        )


class TestDateConversion(unittest.TestCase):
    def test_doy_date_conversion(self):
        self.assertEqual(time.doy_to_date(2010, 1), (2010, 1, 1))
        self.assertEqual(time.doy_to_date(2010, 100), (2010, 4, 10))
        # Leap year
        self.assertEqual(time.doy_to_date(2008, 366), (2008, 12, 31))
        self.assertRaises(ValueError, time.doy_to_date, 2010, 0)
        self.assertRaises(ValueError, time.doy_to_date, 2010, 366)

        self.assertEqual(time.date_to_doy(2010, 1, 1), (2010, 1))
        self.assertEqual(time.date_to_doy(2010, 4, 10), (2010, 100))
        self.assertEqual(time.date_to_doy(2008, 12, 31), (2008, 366))
        self.assertRaises(ValueError, time.date_to_doy, 2010, 12, 32)
        self.assertRaises(ValueError, time.date_to_doy, 2010, 1, 32)

    def test_doy_step(self):
        self.assertEqual(time.doy_step(2010, 1, 1), (2010, 2))
        self.assertEqual(time.doy_step(2010, 2, -1), (2010, 1))
        self.assertEqual(time.doy_step(2010, 365, 1), (2011, 1))
        self.assertEqual(time.doy_step(2011, 1, -1), (2010, 365))

    def test_doy_exists(self):
        self.assertTrue(time.doy_exists(2010, 1))
        self.assertTrue(time.doy_exists(2010, 365))
        self.assertFalse(time.doy_exists(2010, 366))

    def test_season(self):
        assert time.get_season(1) == 0
        assert time.get_season(100) == 1
        assert time.get_season(200) == 2
        assert time.get_season(300) == 3


class TestTimeConversion(unittest.TestCase):
    def test_tai_utc_conversion(self):
        self.assertEqual(time.tai_to_utc(0), (1993, 1, 0))
        self.assertEqual(time.tai_to_utc(time.utc_to_tai(2001, 1, 0)), (2001, 1, 0))
        self.assertEqual(time.utc_to_tai(*time.tai_to_utc(20000000)), 20000000)

    def test_utc_sat_offset(self):
        self.assertEqual(
            time.utc_to_sat_offset(np.array([13.5]), np.array([0]), 13.5), 0
        )
        self.assertEqual(
            time.utc_to_sat_offset(np.array([13.5]), np.array([-15]), 13.5), -1
        )
        assert time.utc_to_sat_offset(np.array([13.5]), np.array([15]), 13.5) == 1
        assert (
            time.utc_to_sat_offset(np.array([13.5]), np.array([345]), 13.5, col="6")
            == -1
        )

        assert (
            time.utc_to_sat_offset(np.array([13.5]), np.array([195]), 13.5, col="5")
            == 13
        )
        assert (
            time.utc_to_sat_offset(np.array([13.5]), np.array([195]), 13.5, col="6")
            == -11
        )


class TestDatetimeConversion(unittest.TestCase):
    def test_ydh_to_datetime(self):
        assert time.ydh_to_datetime(2020, 1, 0) == datetime(2020, 1, 1, 0)
        assert time.ydh_to_datetime(2020, 1, 1) == datetime(2020, 1, 1, 1)
        assert time.ydh_to_datetime(2020, 10, 1) == datetime(2020, 1, 10, 1)
        assert time.ydh_to_datetime(2020, 40, 12) == datetime(2020, 2, 9, 12)

    def test_datetime_to_ydh(self):
        assert time.datetime_to_ydh(datetime(2020, 1, 1, 0)) == (2020, 1, 0)
        assert time.datetime_to_ydh(datetime(2020, 1, 1, 1)) == (2020, 1, 1)
        assert time.datetime_to_ydh(datetime(2020, 1, 10, 1)) == (2020, 10, 1)
        assert time.datetime_to_ydh(datetime(2020, 2, 9, 12)) == (2020, 40, 12)
