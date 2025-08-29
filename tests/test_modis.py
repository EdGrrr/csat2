from csat2 import MODIS, locator, misc
import unittest
import os
import pytest
import numpy as np


class TestMODISBasic(unittest.TestCase):
    def test_granule_fromtext(self):
        gran = MODIS.Granule.fromtext("2015001.0150A")
        assert gran.year == 2015
        assert gran.doy == 1
        assert gran.time == 150
        assert gran.timestr() == "0150"
        assert gran.sat == "A"

    def test_granule_fromfilename(self):
        gran = MODIS.Granule.fromfilename(
            "MOD021KM.A2015001.1855.061.2017318143302.hdf"
        )
        assert gran.year == 2015
        assert gran.doy == 1
        assert gran.time == 1855
        assert gran.sat == "T"

    def test_granule_next(self):
        gran = MODIS.Granule.fromtext("2015001.2150A")
        assert gran.increment().time == 2155
        assert gran.increment(12).time == 2250


@pytest.mark.network
class TestMODISDownload(unittest.TestCase):
    def test_MODIS_download(self):
        gran = MODIS.Granule.fromtext("2015001.1220A")
        # Remove exisiting file
        try:
            fname = gran.get_filename("021KM")
            os.system("rm {}".format(fname))
            try:
                fname = gran.get_filename("021KM")
                raise OSError("Could not delete file")
            except IndexError:
                pass
        except IndexError:
            pass
        gran.download("021KM")
        assert isinstance(gran.get_filename("021KM"), str)

    def test_MODIS_download_geometa(self):
        gran = MODIS.Granule.fromtext("2015001.1220A")
        # Remove exisiting file
        try:
            fname = gran.get_filename("GEOMETA")
            os.system("rm {}".format(fname))
            try:
                fname = gran.get_filename("GEOMETA")
                raise OSError("Could not delete file")
            except IndexError:
                pass
        except IndexError:
            pass
        gran.download_geometa()
        assert isinstance(gran.get_filename("GEOMETA"), str)

    @pytest.mark.skip
    def test_MODIS_download_NRT(self):
        pass

    @pytest.mark.skip
    def test_MODIS_download_geometa_NRT(self):
        pass


class TestMODISGranule(unittest.TestCase):
    def setUp(self):
        # Requires that the download tests have passed
        # Some data is required!
        self.gran = MODIS.Granule.fromtext("2015001.1220A")
        try:
            self.gran.get_filename("021KM")
        except IndexError:
            pytest.skip(
                "Test requires 2015001.1220A (MYD021KM)- run with --runnetwork to download automatically"
            )

    def test_getfilename(self):
        assert len(self.gran.get_filename("021KM")) > 0

    def test_readin_L2_filename(self):
        data = MODIS.readfiles.readin_MODIS_L2_filename(
            self.gran.get_filename("021KM"), names=["Latitude", "Longitude"]
        )
        assert "Latitude" in data.keys()
        assert "Longitude" in data.keys()

    def test_gran_radiances(self):
        for band in [1, 3, 8, 20]:  # Check all the band types
            data = self.gran.get_band_radiance(band=band, refl=False, bowtie_corr=False)
            assert data.max() > 0  # Should always be the case for IR channel
            assert data.min() >= 0

    def test_gran_refl(self):
        for band in [1, 3, 8]:  # Check all the band types
            data = self.gran.get_band_radiance(band=band, refl=True, bowtie_corr=False)
            assert data.max() > 0  # Should always be the case for IR channel
            assert data.min() >= 0

    def test_gran_bt(self):
        for band in [20]:
            data = self.gran.get_band_bt(band, bowtie_corr=False)
            assert data.max() > 0  # Should always be the case for IR channel
            assert data.min() >= 0

    def test_bowtie_corr(self):
        for band in [1, 3, 8, 20]:  # Check all the band types
            data = self.gran.get_band_radiance(band=band, refl=False, bowtie_corr=True)
            assert data.max() > 0  # Should always be the case for IR channel
            assert data.min() >= 0


class BaseMODISLocator(unittest.TestCase):
    locator_type = None  # Subclasses set this
    __test__ = False
    
    def setUp(self):
        # Requires that the download tests have passed
        # Some data is required!
        self.gran = MODIS.Granule.fromtext("2015001.1220A")
        self.polar_gran = MODIS.Granule.fromtext("2015001.1235A")

        assert self.locator_type is not None, "Subclasses must set locator_type"
        
        self.locs = [
            ["Pointe de Grave", (-1.0611, 45.5734), (1774, 1244), (1774, 1947)],
            ["Punta del Marchese", (10.0795, 42.6187), (1227, 719), (1227, 1206)],
            ["Corfu", (19.9210, 39.6243), (782, 144), (782, 462)],
        ]

        # These are locations specifically outside the test granule
        self.fail_locs = [
            ["Isla Angel de la Guarda", (-113.1140, 28.9846)],
            ["Verrazzano Narrows", (-74.0431, 40.6054)],
        ]

    def test_locator_valid(self):
        for loc in self.locs:
            ilon, ilat = loc[1]
            # Check 1km locator
            mlon_1km, mlat_1km = self.gran.locate(
                [[ilon, ilat]],
                locator_type=self.locator_type)[0]
            assert np.sqrt((mlon_1km - loc[2][0]) ** 2 + (mlat_1km - loc[2][1]) ** 2) < 2
            nlon, nlat = self.gran.geolocate([[mlon_1km, mlat_1km]])[0]
            # Allow a 1-pixel mis-registration
            assert misc.geo.haversine(nlon, nlat, ilon, ilat) < 1.41

    def test_locator_valid_rectified(self):
        for loc in self.locs:
            ilon, ilat = loc[1]
            # Check 1km locator
            mlon_1km, mlat_1km = self.gran.locate(
                [[ilon, ilat]],
                locator_type=self.locator_type,
                rectified=True)[0]
            assert np.sqrt((mlon_1km - loc[3][0]) ** 2 + (mlat_1km - loc[3][1]) ** 2) < 2

    def test_locator_list_missing(self):
        # Make sure the correct locations are returned even with missing points
        geo_locs = [l[1] for l in self.locs] + [l[1] for l in self.fail_locs]
        test_mlocs = [l[2] for l in self.locs] + [[-1, -1] for l in self.fail_locs]

        calc_mlocs = self.gran.locate(
            geo_locs,
            locator_type=self.locator_type,
            remove_outside=True)
        for i in range(len(test_mlocs)):
            assert np.sqrt(
                (test_mlocs[i][0] - calc_mlocs[i][0]) ** 2 +
                (test_mlocs[i][1] - calc_mlocs[i][1]) ** 2) < 2

    def test_locator_list_missing_rectified(self):
        # Make sure the correct locations are returned even with missing points
        geo_locs = [l[1] for l in self.locs] + [l[1] for l in self.fail_locs]
        test_mlocs = [l[3] for l in self.locs] + [[-1, -1] for l in self.fail_locs]

        calc_mlocs = self.gran.locate(
            geo_locs,
            locator_type=self.locator_type,
            remove_outside=True, rectified=True)
        for i in range(len(test_mlocs)):
            assert np.sqrt(
                (test_mlocs[i][0] - calc_mlocs[i][0]) ** 2 +
                (test_mlocs[i][1] - calc_mlocs[i][1]) ** 2) < 2

    def test_pattern_locations(self):
        pattern_locs = np.array(
            np.meshgrid(np.arange(5, 2100, 200),
                        np.arange(5, 1350, 130))).T
        geo_locs = self.gran.geolocate(pattern_locs)
        calc_plocs = self.gran.locate(
            geo_locs.reshape((-1, 2)),
            locator_type=self.locator_type
        ).reshape((pattern_locs).shape)
        assert (abs(calc_plocs-pattern_locs) <= 1).all()

    def test_pattern_locations_exact(self):
        pattern_locs = np.array(
            np.meshgrid(np.arange(5, 2100, 200),
                        np.arange(5, 1350, 130))).T
        geo_locs = self.gran.geolocate(pattern_locs)
        calc_plocs = self.gran.locate(
            geo_locs.reshape((-1, 2)),
            locator_type=self.locator_type
        ).reshape((pattern_locs).shape)
        assert (abs(calc_plocs-pattern_locs) ==0 ).all()

    def test_pattern_locations_poles(self):
        pattern_locs = np.array(
            np.meshgrid(np.arange(5, 2100, 200),
                        np.arange(5, 1350, 130))).T
        geo_locs = self.polar_gran.geolocate(pattern_locs)
        calc_plocs = self.polar_gran.locate(
            geo_locs.reshape((-1, 2)),
            locator_type=self.locator_type
        ).reshape((pattern_locs).shape)
        assert (abs(calc_plocs-pattern_locs) <= 1).all()

    def test_pattern_locations_poles_exact(self):
        pattern_locs = np.array(
            np.meshgrid(np.arange(5, 2100, 200),
                        np.arange(5, 1350, 130))).T
        geo_locs = self.polar_gran.geolocate(pattern_locs)
        calc_plocs = self.polar_gran.locate(
            geo_locs.reshape((-1, 2)),
            locator_type=self.locator_type
        ).reshape((pattern_locs).shape)
        assert (abs(calc_plocs-pattern_locs) ==0 ).all()

    def test_locator_invalid(self):
        for loc in self.fail_locs:
            ilon, ilat = loc[1]
            # Check 1km locator
            mlon_1km, mlat_1km = self.gran.locate(
                [[ilon, ilat]],
                locator_type=self.locator_type)[0]
            assert np.isnan(mlon_1km) or (mlon_1km == -1)
            assert np.isnan(mlat_1km) or (mlat_1km == -1)


class TestMODISLocator_BallTree(BaseMODISLocator):
    locator_type = 'BallTree'
    __test__ = True

    @pytest.mark.xfail
    def test_pattern_locations_exact(self):
        super().test_pattern_locations_exact()

    @pytest.mark.xfail
    def test_pattern_locations_poles_exact(self):
        super().test_pattern_locations_poles_exact()


class TestMODISLocator_BallTree1km(BaseMODISLocator):
    locator_type = 'BallTree1km'
    __test__ = True
        

class TestMODISLocator_FullSearch(BaseMODISLocator):
    locator_type = 'FullSearch'
    __test__ = True


class TestMODISLocator_SphereRemap(BaseMODISLocator):
    locator_type = 'SphereRemap'
    __test__ = True


class TestMODISLocator_SphereRemap2km(BaseMODISLocator):
    locator_type = 'SphereRemap2km'
    __test__ = True

    @pytest.mark.xfail
    def test_pattern_locations_exact(self):
        super().test_pattern_locations_exact()

    @pytest.mark.xfail
    def test_pattern_locations_poles_exact(self):
        super().test_pattern_locations_poles_exact()
