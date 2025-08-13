from csat2 import EarthCARE, locator, misc
import unittest
import os
import pytest
import numpy as np
import datetime


class TestEarthCAREBasic(unittest.TestCase):
    def test_granule_fromtext(self):
        gran = EarthCARE.Granule.fromtext("EC.05134B")
        assert gran.orbit == 5134
        assert gran.frame == 'B'

    def test_granule_datetime(self):
        gran = EarthCARE.Granule.fromtext("EC.05134B")
        assert gran.datetime() == datetime.datetime(2025, 4, 23, 21, 45, 7)
        
    def test_granule_fromfilename(self):
        gran = EarthCARE.Granule.fromfilename('ECA_EXAE_ATL_NOM_1B_20250423T214507Z_20250423T220738Z_05134B')
        assert gran.orbit == 5134
        assert gran.frame == 'B'

    def test_granule_next(self):
        gran = EarthCARE.Granule.fromtext("EC.05134B")
        assert gran.next().orbit == 5134
        assert gran.next().frame == 'C'
        assert gran.next(8).orbit == 5135
        assert gran.next(8).frame == 'B'


@pytest.mark.network
class TestEarthCAREDownload(unittest.TestCase):
    def test_GOES_download(self):
        gran = EarthCARE.Granule.fromtext("EC.05234A", baseline='AE')
        try:
            filename = gran.get_filename('ATL_NOM_1B')
            os.remove(filename)
        except FileNotFoundError:
            pass
        gran.download('ATL_NOM_1B')
        assert isinstance(gran.get_filename('ATL_NOM_1B'), str)


class TestEarthCAREGranule(unittest.TestCase):
    def setUp(self):
        # Requires that the download tests have passed
        # Some data is required!
        newfiles = locator.search(
            "EarthCARE",
            "ATL_NOM_1B",
            year=2025,
            doy=120,
            orbit=5234,
            frame='A',
            baseline='AE'
        )
        if len(newfiles) != 1:
            pytest.skip(
                "Test requires ATL_NOM_1B for EC.05234A - run with --runnetwork to download automatically"
            )
        self.gran = EarthCARE.Granule.fromtext('EC.05234A')

    def test_getfilename(self):
        assert len(self.gran.get_filename('ATL_NOM_1B')) > 0

    def test_readin_filename(self):
        data = self.gran.get_variable(
            'ATL_NOM_1B',
            'mie_attenuated_backscatter')['mie_attenuated_backscatter']
        assert data.shape[0] > 10000 # In theory we don't know how large this should be
        assert data.shape[1] == 254
        assert data.max() > 0  # Should always be the case for IR channel
        assert data.min() >= -1e-4 # Some smallish number to account for changes

    def test_lonlat(self):
        lonlat = self.gran.get_lonlat('ATL_NOM_1B')
        assert len(lonlat) == 2
        assert abs(lonlat[0].min() - -93.6) <= 0.1
        assert abs(lonlat[0].max() - -84.7) <= 0.1
        assert abs(lonlat[1].min() - -22.9) <= 0.1
        assert abs(lonlat[1].max() -  23.2) <= 0.1

    def test_locate(self):
        assert self.gran.locate('ATL_NOM_1B', [[-85.2052, -20.4157]]) == 1000

    def test_geolocate(self):
        loc = self.gran.geolocate('ATL_NOM_1B', [1000])
        assert (loc[0] - -85.20529938) <= 1e-6        
        assert (loc[1] - -20.41574287) <= 1e-6
