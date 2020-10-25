from csat2 import GOES, locator, misc
import unittest
import os
import pytest
import numpy as np


class TestGOESBasic(unittest.TestCase):
    def test_granule_fromtext(self):
        gran = GOES.Granule.fromtext('G16.2018002.0000.RadC')
        assert gran.year == 2018
        assert gran.doy == 2
        assert gran.hour == 0
        assert gran.minute == 0
        assert gran.sat == 'G16'

    def test_granule_fromfilename(self):
        gran = GOES.Granule.fromfilename(
            'OR_ABI-L1b-RadC-M3C01_G16_s20180020002199_e20180020004572_c20180020005016.nc')
        assert gran.year == 2018
        assert gran.doy == 2
        assert gran.hour == 0
        assert gran.minute == 0
        assert gran.sat == 'G16'

    def test_granule_next(self):
        gran = GOES.Granule.fromtext('G16.2018002.0000.RadC')
        assert gran.next().minute == 5
        gran = GOES.Granule.fromtext('G16.2018002.0000.RadF')
        assert gran.next().minute == 10
        gran = GOES.Granule.fromtext('G16.2018002.0000.RadM1')
        assert gran.next().minute == 1        


@pytest.mark.network
class TestGOESDownload(unittest.TestCase):
    def test_GOES_download(self):
        gran = GOES.Granule.fromtext('G16.2018002.0000.RadC')
        for ch in [1, 13]:
            # Remove exisiting file
            try:
                fname = gran.get_filename(channel=ch)
                os.system('rm {}'.format(fname))
                try:
                    fname = gran.get_filename(channel=ch)
                    raise OSError('Could not delete file')
                except IndexError:
                    pass
            except IndexError:
                pass
            gran.download(channel=ch)
            assert isinstance(gran.get_filename(channel=ch), str)


class TestGoesGranule(unittest.TestCase):
    def setUp(self):
        # Requires that the download tests have passed
        # Some data is required!
        for ch in [1, 13]:
            newfiles = locator.search(
                'GOES', 'L1b',
                year=2018, doy=2,
                hour=0, minute=2,
                channel=ch,
                sat='G16', area='RadC', mode='*')
            if len(newfiles) != 1:
                pytest.fail('Test requires G16.2018002.000.RadC (ch{ch})- run with --runnetwork to download automatically'.format(ch=ch))
        self.gran = GOES.Granule.fromtext('G16.2018002.0000.RadC')
                
    def test_getfilename(self):
        assert len(self.gran.get_filename(channel=13)) > 0

    def test_readin_filename(self):
        data = GOES.readin_radiances_filename(
            self.gran.get_filename(channel=13))
        assert data.shape == (1500, 2500)
        assert data.max() > 0  # Should always be the case for IR channel
        assert data.min() >= 0

    def test_gran_radiances(self):
        data = self.gran.get_band_radiance(channel=13)
        assert data.shape == (1500, 2500)
        assert data.max() > 0  # Should always be the case for IR channel
        assert data.min() >= 0

    def test_gran_refl(self):
        data = self.gran.get_band_radiance(channel=1, refl=True)
        assert data.shape == (3000, 5000)
        assert data.max()>0   # Should always be the case for IR channel

    def test_gran_bt(self):
        data = self.gran.get_band_bt(channel=13)
        assert data.shape == (1500, 2500)
        assert data.max() > 0  # Should always be the case for IR channel
        assert data.min() >= 0

    def test_get_shape(self):
        assert self.gran.get_shape(channel=13) == (1500, 2500)


class TestGOESLocator(unittest.TestCase):
    def setUp(self):
        # Requires that the download tests have passed
        # Some data is required!
        for ch in [1, 13]:
            newfiles = locator.search(
                'GOES', 'L1b',
                year=2018, doy=2,
                hour=0, minute=2,
                channel=ch,
                sat='G16', area='RadC', mode='*')
            if len(newfiles) != 1:
                pytest.fail('Test requires G16.2018002.000.RadC (ch{ch})- run with --runnetwork to download automatically'.format(ch=ch))
        self.gran = GOES.Granule.fromtext('G16.2018002.0000.RadC')

        self.locs = [
            ['Isla Angel de la Guarda', (-113.1140, 28.9846), (380, 1681)],
            ['Verrazzano Narrows', (-74.0431, 40.6054), (3696, 643)]]

    def test_locator_geometry(self):
        glocator = GOES.granule.GOESLocator(self.gran.get_filename(channel=13))
        lon = -84.690932
        lat = 33.846162
        lonr = np.deg2rad(lon)
        latr = np.deg2rad(lat)
        saty = 0.095340
        satx = -0.024052
        assert np.all(np.isclose(glocator._geodetic_to_sat(lonr, latr),
                                 (satx, saty)))
        assert np.all(np.isclose(
            glocator._geodetic_to_sat(
                np.array([lonr, lonr]), np.array([latr, latr])),
            (np.array([satx, satx]), np.array([saty, saty]))))

        assert np.all(np.isclose(glocator._sat_to_geodetic(satx, saty),
                                 (lonr, latr)))
        assert np.all(np.isclose(
            glocator._sat_to_geodetic(
                np.array([satx, satx]), np.array([saty, saty])),
            (np.array([lonr, lonr]), np.array([latr, latr]))))
        # Check points off Earth fail
        assert np.all(np.isnan(glocator._sat_to_geodetic(-0.2, 0.1)))

    def test_locator_gridding(self):
        # Some identified locations in GOES images
        # Examples in the GOES document appear to be wrong...

        for loc in self.locs:
            ilon, ilat = loc[1]
            # Check 1km locator
            glocator = GOES.granule.GOESLocator(self.gran.get_filename(channel=1))
            slon1km, slat1km = loc[2]
            geoloc = glocator.geolocate(np.array([[slon1km, slat1km]]))[0]
            assert misc.geo.haversine(*geoloc, ilon, ilat) < 1
            sloc = glocator.locate(np.array([[ilon, ilat]]))[0]
            assert np.sqrt((slon1km-sloc[0])**2 + (slat1km-sloc[1])**2) < 2
            
            # Check 2km locator
            glocator = GOES.granule.GOESLocator(self.gran.get_filename(channel=13))
            slon2km, slat2km = slon1km//2, slat1km//2
            geoloc = glocator.geolocate(np.array([[slon2km, slat2km]]))[0]
            assert misc.geo.haversine(*geoloc, ilon, ilat) < 2
            sloc = glocator.locate(np.array([[ilon, ilat]]))[0]
            assert np.sqrt((slon2km-sloc[0])**2 + (slat2km-sloc[1])**2) < 2

    def test_locator_gridding_interpolate(self):
        # Some identified locations in GOES images
        # Examples in the GOES document appear to be wrong...

        for loc in self.locs:
            ilon, ilat = loc[1]
            # Check 1km locator
            glocator = GOES.granule.GOESLocator(self.gran.get_filename(channel=1))
            slon1km, slat1km = loc[2]
            geoloc = glocator.geolocate(np.array([[slon1km, slat1km]]), interp=True)[0]
            # Must be within 1km
            assert misc.geo.haversine(*geoloc, ilon, ilat) < 1
            sloc = glocator.locate(np.array([[ilon, ilat]]), interp=True)[0]
            assert np.sqrt((slon1km-sloc[0])**2 + (slat1km-sloc[1])**2) <= np.sqrt(2)
            
            # Check 2km locator
            glocator = GOES.granule.GOESLocator(self.gran.get_filename(channel=13))
            slon2km, slat2km = slon1km//2, slat1km//2
            geoloc = glocator.geolocate(np.array([[slon2km, slat2km]]), interp=True)[0]
            # Must be within 2km
            assert misc.geo.haversine(*geoloc, ilon, ilat) < 2
            sloc = glocator.locate(np.array([[ilon, ilat]]), interp=True)[0]
            assert np.sqrt((slon2km-sloc[0])**2 + (slat2km-sloc[1])**2) <= np.sqrt(2)

    def test_locator_gridding_failures(self):
        '''Check for the correct returns when the input point is outside
        the domain.'''
        glocator = GOES.granule.GOESLocator(self.gran.get_filename(channel=13))
        data = self.gran.get_band_bt(13)
        fail_locs = [[0, 0],
                     [0, -90],
                     [0, 90],
                     [-170, 0]]
        for fl in fail_locs:
            sloc = glocator.locate(np.array([fl]))[0]
            assert np.all(sloc == -999999)
            with pytest.raises(Exception):
                data[sloc[1], sloc[0]]

    def test_locator_gridding_interpolate_failures(self):
        glocator = GOES.granule.GOESLocator(self.gran.get_filename(channel=13))
        data = self.gran.get_band_bt(13)
        fail_locs = [[0, 0],
                     [0, -90],
                     [0, 90],
                     [-170, 0]]
        for fl in fail_locs:
            sloc = glocator.locate(np.array([fl]), interp=True)[0]
            assert np.all(np.isnan(sloc))
            with pytest.raises(Exception):
                data[sloc[1], sloc[0]]
