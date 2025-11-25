import os
import numpy as np
import pytest
from csat2.CERES.granule import Granule 


## define constants for testing

TEST_YEAR = 2020
TEST_DOY = 1       # 1 Jan 2020
TEST_TIME = 12     # midday UTC/ GMT
TEST_VAR = "obs_cld_amount"  


# ---------------------------
# FIXTURES
# ---------------------------

@pytest.fixture(scope="module")
def granule(tmp_path_factory):
    """
    Create a Granule object, ensure the file is deleted, then downloaded fresh.
    """
    g = Granule(TEST_YEAR, TEST_DOY, TEST_TIME)

    fileloc = g.get_fileloc()

    # Delete old file so test is clean
    if os.path.exists(fileloc):
        os.remove(fileloc)

    # Fresh download
    g.download(force_redownload=True)

    assert os.path.exists(fileloc), "CERES file did not download correctly."

    return g


# ---------------------------
# BASIC TESTS
# ---------------------------

def test_download(granule):
    """File should exist after download."""
    assert granule.check()


def test_list_variables(granule):
    """Ensure we can list variables and the test var is present."""
    vars = granule.list_variables()
    assert isinstance(vars, list)
    assert TEST_VAR in vars


def test_get_variable(granule):
    """Ensure get_variable returns an xarray.DataArray."""
    da = granule.get_variable(TEST_VAR)
    assert hasattr(da, "shape")
    assert da.ndim >= 2


def test_get_lonlat(granule):
    """Basic check lon/lat load correctly."""
    lon, lat = granule.get_lonlat()

    assert lon.ndim == 1
    assert lat.ndim == 1
    assert lon.size == 360
    assert lat.size == 180


def test_linear_geolocation(granule):
    """Test that linear geolocation returns finite values."""
    lons = np.array([0.0, 45.0])
    lats = np.array([0.0, -10.0])

    vals = granule.geolocate(TEST_VAR, lons, lats, method="linear")

    assert np.isfinite(vals).all()


def test_xarray_geolocation(granule):
    """Same test but with xarray method."""
    lons = np.array([0.0, 45.0])
    lats = np.array([0.0, -10.0])

    vals = granule.geolocate(TEST_VAR, lons, lats, method="xarray")

    assert np.isfinite(vals).all()


def test_linear_matches_xarray(granule):
    """Basic consistency check between methods."""
    lons = np.array([0.0, 179.0, -123.4])
    lats = np.array([10.0, -45.0, 60.0])

    v_lin = granule.geolocate(TEST_VAR, lons, lats, method="linear")
    v_xr  = granule.geolocate(TEST_VAR, lons, lats, method="xarray")

    # They should match because both do nearest-neighbour on a regular grid
    np.testing.assert_allclose(v_lin, v_xr, rtol=1e-5, atol=1e-5)