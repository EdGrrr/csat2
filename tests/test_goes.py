from csat2 import GOES
import unittest
import os


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


class TestGOESDownload(unittest.TestCase):
    def test_GOES_download(self):
        gran = GOES.Granule.fromtext('G16.2018002.0000.RadC')
        # Remove exisiting file
        try:
            fname = gran.get_filename(channel=13)
            os.system('rm {}'.format(fname))
            try:
                fname = gran.get_filename(channel=13)
                raise OSError('Could not delete file')
            except IndexError:
                pass
        except IndexError:
            pass
        gran.download(channel=13)
        assert isinstance(gran.get_filename(channel=13), str)

        
class TestGoesGranule(unittest.TestCase):
    # Requires that the download tests have passed
    # Some data is required!
    pass


class TestGOESLocator(unittest.TestCase):
    pass

