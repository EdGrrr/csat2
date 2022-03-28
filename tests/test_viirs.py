from csat2 import VIIRS, misc
import unittest
import os
import pytest
import numpy as np


class TestVIIRSBasic(unittest.TestCase):
    def test_granule_fromtext(self):
        gran = VIIRS.Granule.fromtext('2015080.0806N')
        assert gran.year == 2015
        assert gran.doy == 80
        assert gran.time == 806
        assert gran.timestr() == '0806'
        assert gran.sat == 'N'

    def test_granule_fromfilename(self):
        gran = VIIRS.Granule.fromfilename(
            'VNP02IMG.A2015080.0806.001.2017260093209.nc')
        assert gran.year == 2015
        assert gran.doy == 80
        assert gran.time == 806
        assert gran.sat == 'N'

    def test_granule_next(self):
        gran = VIIRS.Granule.fromtext('2015080.0806N')
        assert gran.increment().time == 812
        assert gran.increment(12).time == 918


@pytest.mark.network
class TestVIIRSDownload(unittest.TestCase):
    def test_VIIRS_download(self):
        gran = VIIRS.Granule.fromtext('2015080.0806N')
        # Remove exisiting file
        try:
            fname = gran.get_filename('02IMG')
            os.system('rm {}'.format(fname))
            try:
                fname = gran.get_filename('02IMG')
                raise OSError('Could not delete file')
            except IndexError:
                pass
        except IndexError:
            pass
        gran.download_product('02IMG')
        assert isinstance(gran.get_filename('02IMG'), str)

    def test_VIIRS_download_geometa(self):
        gran = VIIRS.Granule.fromtext('2015080.0806N')
        # Remove exisiting file
        try:
            fname = gran.get_filename('GEOMETA')
            os.system('rm {}'.format(fname))
            try:
                fname = gran.get_filename('GEOMETA')
                raise OSError('Could not delete file')
            except IndexError:
                pass
        except IndexError:
            pass
        gran.download_geometa()
        assert isinstance(gran.get_filename('GEOMETA'), str)

    @pytest.mark.skip
    def test_VIIRS_download_NRT(self):
        pass

    @pytest.mark.skip
    def test_VIIRS_download_geometa_NRT(self):
        pass


class TestVIIRSGranule(unittest.TestCase):
    def setUp(self):
        # Requires that the download tests have passed
        # Some data is required!
        self.gran = VIIRS.Granule.fromtext('2015080.0806N')
        try:
            self.gran.get_filename('02IMG')
        except IndexError:
            pytest.skip(
                'Test requires 2015080.0806N (VNP02IMG)- run with --runnetwork to download automatically')

    def test_getfilename(self):
        assert len(self.gran.get_filename('02IMG')) > 0

    def test_readin_L2_filename(self):
        data = VIIRS.readfiles.readin_VIIRS_L2_filename(
            self.gran.get_filename('02IMG'), names=['I01', 'I02'])
        assert 'I01' in data.keys()
        assert 'I02' in data.keys()

    def test_gran_radiances(self):
        for band in ['I01', 'I05']:  # Check all the band types
            data = self.gran.get_band_radiance(
                band=band, refl=False, bowtie_corr=False)
            assert np.nanmax(data) > 0
            assert np.nanmin(data) >= 0

    def test_gran_refl(self):
        for band in ['I01']:
            data = self.gran.get_band_radiance(
                band=band, refl=True, bowtie_corr=False)
            assert np.nanmax(data) > 0
            assert np.nanmin(data) >= 0

    def test_gran_bt(self):
        for band in ['I05']:
            data = self.gran.get_band_bt(band, bowtie_corr=False)
            assert np.nanmax(data) > 0
            assert np.nanmin(data) >= 0

    def test_bowtie_corr(self):
        for band in ['I01']:  # Check all the band types
            data = self.gran.get_band_radiance(
                band=band, refl=False, bowtie_corr=True)
            assert data.max() > 0
            assert data.min() >= 0


@pytest.mark.skip
class TestVIIRSLocator(unittest.TestCase):
    def setUp(self):
        # Requires that the download tests have passed
        # Some data is required!
        self.gran = VIIRS.Granule.fromtext('2015001.2150A')

        self.locs = [
            ['Isla Angel de la Guarda', (-113.1140, 28.9846), (380, 1681)],
            ['Verrazzano Narrows', (-74.0431, 40.6054), (3696, 643)]]

    def test_locator_geometry(self):
        mlocator = VIIRS.granule.GOESLocator(
            self.gran.get_filename(channel=13))
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
            glocator = GOES.granule.GOESLocator(
                self.gran.get_filename(channel=1))
            slon1km, slat1km = loc[2]
            geoloc = glocator.geolocate(np.array([[slon1km, slat1km]]))[0]
            assert misc.geo.haversine(*geoloc, ilon, ilat) < 1
            sloc = glocator.locate(np.array([[ilon, ilat]]))[0]
            assert np.sqrt((slon1km-sloc[0])**2 + (slat1km-sloc[1])**2) < 2

            # Check 2km locator
            glocator = GOES.granule.GOESLocator(
                self.gran.get_filename(channel=13))
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
            glocator = GOES.granule.GOESLocator(
                self.gran.get_filename(channel=1))
            slon1km, slat1km = loc[2]
            geoloc = glocator.geolocate(
                np.array([[slon1km, slat1km]]), interp=True)[0]
            # Must be within 1km
            assert misc.geo.haversine(*geoloc, ilon, ilat) < 1
            sloc = glocator.locate(np.array([[ilon, ilat]]), interp=True)[0]
            assert np.sqrt((slon1km-sloc[0])**2 +
                           (slat1km-sloc[1])**2) <= np.sqrt(2)

            # Check 2km locator
            glocator = GOES.granule.GOESLocator(
                self.gran.get_filename(channel=13))
            slon2km, slat2km = slon1km//2, slat1km//2
            geoloc = glocator.geolocate(
                np.array([[slon2km, slat2km]]), interp=True)[0]
            # Must be within 2km
            assert misc.geo.haversine(*geoloc, ilon, ilat) < 2
            sloc = glocator.locate(np.array([[ilon, ilat]]), interp=True)[0]
            assert np.sqrt((slon2km-sloc[0])**2 +
                           (slat2km-sloc[1])**2) <= np.sqrt(2)

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
