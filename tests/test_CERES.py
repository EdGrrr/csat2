import os
import numpy as np
import pytest
import csat2.CERES


## constants for the test, this could be made 'random' for future testing

TEST_YEAR = 2020
TEST_DOY = 1       # 1 Jan 2020
TEST_TIME = 12     # midday UTC/ GMT
TEST_VAR = "obs_cld_amount"  


@pytest.fixture(scope="module")
def granule(tmp_path_factory):
    """
    Create a Granule object, ensure the file is deleted, then downloaded fresh.
    """
    granule = csat2.CERES.Granule(TEST_YEAR, TEST_DOY, TEST_TIME)

    fileloc = granule.get_fileloc()

    # Delete old file so test is clean
    if os.path.exists(fileloc):
        os.remove(fileloc)

    # Fresh download
    granule.download(force_redownload=True)

    assert os.path.exists(fileloc), "CERES file did not download correctly."

    return granule

## defie some tests
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

    N = 10 # number of random points to test
    lons = np.random.uniform(-180, 180, N)
    lats = np.random.uniform(-90, 90, N)

    eps = 1e-4 ## this is so that a point exactly at the midpoint between two grid points is nudged either side so that both methods should (in theopry) give the same results, since xarrays tie breaker logic was causing this test to fdails since if two points are equally close xarray chooses the lower index one

    def nudge(coords,grid_spacing=1.0,eps=eps):
        frac_part = np.mod(coords,grid_spacing)
        mask = np.isclose(frac_part,0.0,atol=1e-8)
        coords[mask] += np.random.choice([-1,1],size=np.sum(mask))*eps
        return coords
    lons = nudge(lons)
    lats = nudge(lats)

    v_lin = granule.geolocate(TEST_VAR, lons, lats, method="linear")
    v_xr  = granule.geolocate(TEST_VAR, lons, lats, method="xarray")

    # They should match because both do nearest-neighbour on a regular grid
    np.testing.assert_allclose(v_lin, v_xr, rtol=1e-5, atol=1e-5)

