from csat2 import ECMWF, locator, misc
import unittest
import os
import pytest
import numpy as np


@pytest.mark.network
class TestECMWFDownload(unittest.TestCase):
    def test_ECMWFDownload(self):
        ECMWF.download(
            2020, 1,
            ['U-wind-component', 'V-wind-component'],
            level='1000hPa',
            resolution='1grid',
            days=[1, 2],
            lst=False)
        assert ECMWF.check(
            2020, 1,
            ['U-wind-component', 'V-wind-component'],
            level='1000hPa',
            resolution='1grid',
            days=[1, 2])[0]
        ECMWF.download(
            2020, 1,
            ['Temperature'],
            level='700hPa',
            resolution='1grid',
            days=[1, 2],
            lst=False)
        assert ECMWF.check(
            2020, 1,
            ['Temperature'],
            level='700hPa',
            resolution='1grid',
            days=[1, 2])[0]

    def test_ECMWFDownloadLST(self):
        ECMWF.download(
            2020, 1,
            ['Temperature'],
            level='1000hPa',
            resolution='1grid',
            days=[1, 2],
            lst=True)
        assert ECMWF.check(
            2020, 1,
            ['Temperature'],
            level='1000hPa',
            resolution='1grid',
            days=[1, 2])[0]
        assert ECMWF.check(
            2020, 1,
            ['Temperature'],
            level='1000hPa',
            resolution='1grid',
            days=[1, 2],
            time='LST')[0]


class TestECMWFBasic(unittest.TestCase):
    def setUp(self):
        # Requires that the download tests have passed
        # Some data is required!
        if not ECMWF.check(2020, 1,
                           ['Temperature'],
                           level='1000hPa',
                           resolution='1grid',
                           days=[1, 2])[0]:
            pytest.skip('Requires 1000hPa Temperature data')
        if not ECMWF.check(2020, 1,
                           ['Temperature'],
                           level='700hPa',
                           resolution='1grid',
                           time='LST',
                           days=[1, 2])[0]:
            pytest.skip('Requires 700hPa Temperature data')
            
    def rel_err(self, calc, actual):
        return (calc-actual)/actual
        
    def test_ERA5_readin(self):
        data = ECMWF.readin(
            'ERA5', 2020, 1,
            'Temperature', level='1000hPa',
            time='timed', resolution='1grid')
        assert data.shape == (8, 360, 180)
        # Check the Geolocation follows the MODIS L3 grid
        assert np.all(data.lon == np.arange(-179.5, 180, 1))
        assert np.all(data.lat == np.arange(89.5, -90, -1))
        # Check the Arctic is colder (NH hemisphere winter)
        zonal_mean = data.mean(axis=(0, 1))
        assert np.mean(zonal_mean[:5]) < np.mean(zonal_mean[:-5])
        # If the data longitude is correct, Australia should be much hotter that 
        timmean = data.mean(axis=0)
        assert timmean[320, 115] > timmean[140, 115]

    def test_ERA5_readin_LST(self):
        data = ECMWF.readin(
            'ERA5', 2020, 1,
            'Temperature', level='1000hPa',
            time='LST', resolution='1grid')
        assert data.shape == (4, 360, 180)
        # Check the Geolocation follows the MODIS L3 grid
        assert np.all(data.lon == np.arange(-179.5, 180, 1))
        assert np.all(data.lat == np.arange(89.5, -90, -1))
        # Check the Arctic is colder (NH hemisphere winter)
        zonal_mean = data.mean(axis=(0, 1))
        assert np.mean(zonal_mean[:5]) < np.mean(zonal_mean[:-5])
        # If the data longitude is correct, Australia should be much hotter that 
        timmean = data.mean(axis=0)
        assert timmean[320, 115] > timmean[140, 115]

    def test_eis_calculations(self):
        t700 = 270
        t1000 = 290
        rh1000 = 80
        # Less than 1% error in component calculations (email from Rob Wood)
        # LTS
        assert self.rel_err(
            ECMWF.ECMWF._calc_lts(t700, t1000),
            9.78815) < 0.01
        # MALR
        assert self.rel_err(
            ECMWF.ECMWF._sat_adiabatic_potential_gradient(
                0.5*(t700+t1000), 85000, rh1000),
            0.00444213) < 0.01
        # LCL
        assert self.rel_err(
            ECMWF.ECMWF._calc_lcl(t1000, rh1000),
            433.884) < 0.01
        # EIS - This is not Rob's number as I am not sure quite how he gets it
        assert self.rel_err(
            ECMWF.ECMWF._calc_eis(t700, t1000, rh1000)[1],
            -1.8) < 0.01

    def test_ERA5_derived(self):
        # Check doesn't throw error
        data = ECMWF.readin(
            'ERA5', 2020, 1,
            'LTS', level='1000hPa',
            time='timed', resolution='1grid')
        assert data.shape == (8, 360, 180)
        
        data = ECMWF.readin(
            'ERA5', 2020, 1,
            'EIS', level='1000hPa',
            time='timed', resolution='1grid')
        assert data.shape == (8, 360, 180)
        
        data = ECMWF.readin(
            'ERA5', 2020, 1,
            'EIS', level='1000hPa',
            time='timed', resolution='1grid', with_rh=False)
        assert data.shape == (8, 360, 180)


class TestECMWFData(unittest.TestCase):
    def setUp(self):
        # Requires that the download tests have passed
        # Some data is required!
        assert ECMWF.check(
            2020, 1,
            ['U-wind-component', 'V-wind-component', 'Temperature'],
            level='1000hPa',
            resolution='1grid',
            days=[1, 2])[0]
        assert ECMWF.check(
            2020, 1,
            ['Temperature'],
            level='700hPa',
            resolution='1grid',
            days=[1, 2])[0]

    def test_era5data_timeslice(self):
        tdata = ECMWF.ERA5Data('Temperature', '1000hPa', res='1grid')
        temp = tdata.get_data_time(
            misc.time.ydh_to_datetime(2020, 1, 6))
        assert temp.shape == (360, 180)
        assert np.all(temp.lon == np.arange(-179.5, 180, 1))
        assert np.all(temp.lat == np.arange(89.5, -90, -1))
        # Check the Arctic is colder (NH hemisphere winter)
        zonal_mean = temp.mean(axis=0)
        assert np.mean(zonal_mean[:5]) < np.mean(zonal_mean[:-5])
        # If the data longitude is correct, Australia should be much hotter that 
        assert temp[320, 115] > temp[140, 115]

        # Ensure the day update works
        temp2 = tdata.get_data_time(
            misc.time.ydh_to_datetime(2020, 2, 6))
        # Check the data has changed!
        assert np.any(temp != temp2)
        # Now some other checks
        assert temp2.shape == (360, 180)
        assert np.all(temp2.lon == np.arange(-179.5, 180, 1))
        assert np.all(temp2.lat == np.arange(89.5, -90, -1))
        # Check the Arctic is colder (NH hemisphere winter)
        zonal_mean = temp2.mean(axis=0)
        assert np.mean(zonal_mean[:5]) < np.mean(zonal_mean[:-5])
        # If the data longitude is correct, Australia should be much hotter that 
        assert temp2[320, 115] > temp2[140, 115]

    def test_era5data_timeslice_timeinterp(self):
        tdata = ECMWF.ERA5Data('Temperature', '1000hPa', res='1grid', linear_interp=False)
        temp_pre = tdata.get_data_time(
            misc.time.ydh_to_datetime(2020, 1, 0))
        temp_post = tdata.get_data_time(tdata.time[1])
        # Ensure these are not the same data!
        assert np.any(temp_pre != temp_post)
        
        tdatai = ECMWF.ERA5Data('Temperature', '1000hPa', res='1grid', linear_interp='time')
        temp = tdatai.get_data_time(
            tdata.time[0]+0.5*(tdata.time[1]-tdata.time[0]))
        assert temp.shape == (360, 180)
        assert np.all(temp.lon == np.arange(-179.5, 180, 1))
        assert np.all(temp.lat == np.arange(89.5, -90, -1))
        # Check the Arctic is colder (NH hemisphere winter)
        zonal_mean = temp.mean(axis=0)
        assert np.mean(zonal_mean[:5]) < np.mean(zonal_mean[:-5])
        # If the data longitude is correct, Australia should be much hotter that 
        assert temp[320, 115] > temp[140, 115]
        # Check that the interpolation is done correctly
        assert np.max(np.abs(0.5*(temp_pre+temp_post) - temp)) < 0.01
        # A harder example to ensure the interal weighting is correct
        for tw in [0, 0.1, 0.9, 0.99]:
                temp = tdatai.get_data_time(
                    tdata.time[0]+tw*(tdata.time[1]-tdata.time[0]))
                assert np.max(np.abs(tw*(temp_post-temp_pre)+temp_pre - temp)) < 0.01

    def test_era5data_points(self):
        tdata = ECMWF.ERA5Data('Temperature', '1000hPa', res='1grid')
        # Australia is warmer than the Pacific
        temp = tdata.get_data(
            [140, -40], [-25, -25],
            misc.time.ydh_to_datetime(2020, 1, 6)) 
        assert temp[0] > temp[1]

        # North pole is colder than south pole
        temp = tdata.get_data(
            [0, 0], [89, -89],
            misc.time.ydh_to_datetime(2020, 1, 6)) 
        assert temp[0] < temp[1]

    def test_era5data_points_simple(self):
        tdata = ECMWF.ERA5Data('Temperature', '1000hPa', res='1grid')
        # Australia is warmer than the Pacific
        temp = tdata.get_data(
            [140, -40], [-25, -25],
            misc.time.ydh_to_datetime(2020, 1, 6), simple=True) 
        assert temp[0] > temp[1]

        # North pole is colder than south pole
        temp = tdata.get_data(
            [0, 0], [89, -89],
            misc.time.ydh_to_datetime(2020, 1, 6), simple=True) 
        assert temp[0] < temp[1]

    def test_era5data_points_timeinterp(self):
        locs_lon = [140, -40, 140.5, -40.5]
        locs_lat = [-25, -25, -25.5, -25.5]

        tdata = ECMWF.ERA5Data('Temperature', '1000hPa', res='1grid', linear_interp=False)
        _ = tdata.get_data(
            locs_lon, locs_lat,
            misc.time.ydh_to_datetime(2020, 1, 0))
        temp_pre = tdata.get_data(
            locs_lon, locs_lat,
            tdata.time[0])
        temp_post = tdata.get_data(
            locs_lon, locs_lat,
            tdata.time[1]) 
        assert np.any(temp_pre != temp_post)
        
        tdatai = ECMWF.ERA5Data('Temperature', '1000hPa', res='1grid', linear_interp='time')
        for tw in [0, 0.1, 0.5, 0.9, 0.99]:
            temp = tdatai.get_data(
                locs_lon, locs_lat,
                tdata.time[0]+tw*(tdata.time[1]-tdata.time[0]))
            assert np.max(np.abs(tw*(temp_post-temp_pre)+temp_pre - temp)) < 0.01

    def test_era5data_points_spaceinterp(self):
        locs_lon = np.array([140.7, -40.5, 0.5, 0, 359.7, 30.5, 29.8, 30])
        locs_lat = np.array([-25.2, -25.5, 0.5, 0.5, 0.5, -89.9, -89.9, 89.9])

        def rdown(data):
            return np.floor(data-0.5)+0.5
        
        tdata = ECMWF.ERA5Data('Temperature', '1000hPa', res='1grid', linear_interp=False)
        ttime = misc.time.ydh_to_datetime(2020, 1, 0)
        
        temp_ll = tdata.get_data(
            rdown(locs_lon), rdown(locs_lat).clip(-90, 90),
            ttime)
        temp_lr = tdata.get_data(
            rdown(locs_lon+1), rdown(locs_lat).clip(-90, 90),
            ttime)
        temp_ul = tdata.get_data(
            rdown(locs_lon), rdown(locs_lat+1).clip(-90, 90),
            ttime)
        temp_ur = tdata.get_data(
            rdown(locs_lon+1), rdown(locs_lat+1).clip(-90, 90),
            ttime) 
        weights_lon = 1-(locs_lon-rdown(locs_lon))
        weights_lat = 1-(locs_lat-rdown(locs_lat))

        temp_interp = (weights_lat *
                       (weights_lon*temp_ll + (1-weights_lon)*temp_lr) +
                       (1-weights_lat) *
                       (weights_lon*temp_ul + (1-weights_lon)*temp_ur))
        
        tdatai = ECMWF.ERA5Data('Temperature', '1000hPa', res='1grid', linear_interp='space')
        tempi = tdatai.get_data(locs_lon, locs_lat, ttime)

        assert np.all(np.isclose(temp_interp, tempi))

