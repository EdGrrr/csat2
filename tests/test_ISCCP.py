from csat2 import ISCCP
import unittest
import os
import pytest
import numpy as np



class TestISCCP(unittest.TestCase):
    @pytest.mark.network
    def test_ISCCP_download(self):
        """Test downloading an ISCCP file. Check that the local file location matches what is expected"""
        year = 2015
        doy = 1
        time = 0
        product = 'isccp-basic'
        version = 'hgg'
        # initialize granule onject
        gran = ISCCP.Granule(year,doy,time,product,version)
        target_path = gran.get_fileloc()
        #remove instance of the local file if it has previously been downloaded
        if os.path.exists(target_path):
            os.remove(target_path)

        local_path = gran.check(download=True)

        # Check if the file was downloaded
        self.assertTrue(os.path.exists(local_path), f"File not found at {local_path}")