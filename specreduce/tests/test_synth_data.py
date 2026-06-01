import numpy as np
import pytest
from astropy import units as u
from astropy.modeling import models
from astropy.nddata import CCDData
from astropy.wcs import WCS

from specreduce.utils.synth_data import SynthImage


def test_empty_image_is_zeros():
    arr = SynthImage(nx=50, ny=20).to_array()
    assert arr.shape == (20, 50)
    assert np.all(arr == 0)


def test_add_background_constant():
    arr = SynthImage(nx=50, ny=20).add_background(7).to_array()
    assert arr.shape == (20, 50)
    assert np.allclose(arr, 7)


def test_add_background_stacks():
    arr = SynthImage(nx=10, ny=10).add_background(3).add_background(4).to_array()
    assert np.allclose(arr, 7)


def test_immutability():
    base = SynthImage(nx=10, ny=10).add_background(5)
    derived = base.add_background(2)
    assert np.allclose(base.to_array(), 5)
    assert np.allclose(derived.to_array(), 7)
    assert base is not derived


def test_to_ccddata_no_wcs():
    ccd = SynthImage(nx=30, ny=15).add_background(5).to_ccddata()
    assert isinstance(ccd, CCDData)
    assert ccd.unit == u.count
    assert ccd.data.shape == (15, 30)
    assert ccd.wcs is None


def test_add_source_adds_flux_at_trace():
    base = SynthImage(nx=200, ny=100).add_background(5)
    src = base.add_source(profile=models.Gaussian1D(amplitude=100, stddev=10))
    base_arr = base.to_array()
    src_arr = src.to_array()
    # source adds flux somewhere; total flux strictly increases
    assert src_arr.sum() > base_arr.sum()
    assert src_arr.max() > base_arr.max()


def test_add_source_default_profile():
    arr = SynthImage(nx=200, ny=100).add_source().to_array()
    assert arr.shape == (100, 200)
    assert arr.max() > 0


def test_add_source_stackable():
    one = SynthImage(nx=200, ny=100).add_source(
        profile=models.Gaussian1D(amplitude=100, stddev=10)
    )
    two = one.add_source(profile=models.Gaussian1D(amplitude=100, stddev=10))
    assert two.to_array().sum() > one.to_array().sum()


def _linear_wcs(nx, extent=(3000, 6000), wave_unit=u.Angstrom):
    wcs = WCS(naxis=2)
    wcs.wcs.ctype[0] = "WAVE"
    wcs.wcs.ctype[1] = "PIXEL"
    wcs.wcs.cunit[0] = wave_unit
    wcs.wcs.cunit[1] = u.pixel
    wcs.wcs.crval[0] = extent[0]
    wcs.wcs.cdelt[0] = (extent[1] - extent[0]) / nx
    wcs.wcs.crval[1] = 0
    wcs.wcs.cdelt[1] = 1
    return wcs


# --- local (no network) validation tests ---

def test_arc_requires_wcs_or_extent():
    img = SynthImage(nx=100, ny=50, extent=None).add_arcs()
    with pytest.raises(ValueError, match="Must specify either a wavelength extent or a WCS"):
        img.to_array()


def test_arc_extent_wrong_length():
    img = SynthImage(nx=100, ny=50, extent=[1, 2, 3]).add_arcs()
    with pytest.raises(ValueError, match="Wavelength extent must be of length 2"):
        img.to_array()


def test_arc_bad_wave_unit():
    img = SynthImage(nx=100, ny=50, extent=[100, 300], wave_unit=u.pixel).add_arcs()
    with pytest.raises(ValueError, match="Wavelength unit must be a length unit"):
        img.to_array()


def test_arc_bad_tilt_func():
    wcs = _linear_wcs(100)
    img = SynthImage(nx=100, ny=50, wcs=wcs).add_arcs(tilt_func=models.Gaussian1D)
    with pytest.raises(
        ValueError, match="The only tilt functions currently supported are 1D polynomials"
    ):
        img.to_array()


def test_provided_wcs_must_have_spectral_axis():
    wcs = WCS(naxis=2)
    wcs.wcs.ctype[0] = "PIXEL"
    wcs.wcs.ctype[1] = "PIXEL"
    img = SynthImage(nx=100, ny=50, wcs=wcs).add_arcs()
    with pytest.raises(ValueError, match="Provided WCS must have a spectral axis"):
        img.to_array()


def test_provided_wcs_must_be_2d():
    wcs = WCS(naxis=1)
    wcs.wcs.ctype[0] = "WAVE"
    wcs.wcs.cunit[0] = u.Angstrom
    img = SynthImage(nx=100, ny=50, wcs=wcs).add_arcs()
    with pytest.raises(ValueError, match="WCS must have NAXIS=2 for a 2D image"):
        img.to_array()


# --- network tests (load pypeit line lists) ---

@pytest.mark.remote_data
@pytest.mark.filterwarnings("ignore:No observer defined on WCS")
def test_arc_image_with_extent_builds_wcs():
    ccd = SynthImage(nx=300, ny=100, extent=(3500, 7000)).add_arcs(["HeI"]).to_ccddata()
    assert isinstance(ccd, CCDData)
    assert ccd.data.shape == (100, 300)
    assert ccd.wcs is not None
    assert ccd.data.max() > 0


@pytest.mark.remote_data
@pytest.mark.filterwarnings("ignore:No observer defined on WCS")
def test_arc_image_supplied_wcs_and_tilt():
    wcs = _linear_wcs(300)
    tilt = models.Chebyshev1D(degree=2, c0=50, c1=0, c2=100)
    ccd = SynthImage(nx=300, ny=100, wcs=wcs).add_arcs(["HeI"], tilt_func=tilt).to_ccddata()
    assert ccd.data.shape == (100, 300)


@pytest.mark.remote_data
@pytest.mark.filterwarnings("ignore:No observer defined on WCS")
def test_arcs_stackable():
    wcs = _linear_wcs(300)
    one = SynthImage(nx=300, ny=100, wcs=wcs).add_arcs(["HeI"])
    two = one.add_arcs(["NeI"])
    assert two.to_array().sum() > one.to_array().sum()
