'''This file can be useful to test for download/access issues, as it had debug logging already set'''
import logging
logging.basicConfig(level=logging.DEBUG)
from csat2 import EarthCARE

# Quick test - this takes about a minute on hardin
# This opens a variable from the ACM_CAP product via the streaming api
gran = EarthCARE.Granule(6062, 'E', stream=True)
ds = gran.get_variable('ACM_CAP_2B', 'longitude', baseline='BA')

# Opening a variable via direct download
gran = EarthCARE.Granule(6062, 'E', stream=False)
gran.download('ACM_CAP_2B', baseline='BA', force_redownload=True)
ds = gran.get_variable('ACM_CAP_2B', 'ice_normalized_number_concentration', baseline='BA')
