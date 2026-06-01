import numpy as np
import pytest
from astropy import units as u
from astropy.modeling import models
from astropy.nddata import CCDData
from astropy.utils.exceptions import AstropyDeprecationWarning
from astropy.wcs import WCS

from specutils import Spectrum

from specreduce.utils.synth_data import SynthImage
from specreduce.utils.synth_data import (
    make_2d_trace_image,
    make_2d_arc_image,
    make_2d_spec_image,
)


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


def test_add_source_spectrum_modulates_and_normalizes():
    nx, ny = 100, 50
    extent = (3000, 6000)
    wave = np.linspace(extent[0], extent[1], 400) * u.AA
    flux = np.exp(-0.5 * ((wave.value - 4500) / 300) ** 2) * u.count
    spectrum = Spectrum(flux=flux, spectral_axis=wave)
    img = (
        SynthImage(nx=nx, ny=ny, extent=extent)
        .add_source(
            profile=models.Gaussian1D(amplitude=10, stddev=5),
            trace_coeffs={"c0": 0, "c1": 0, "c2": 0},
            spectrum=spectrum,
        )
        .to_array()
    )
    # peak value equals the profile amplitude (normalized flux peak == 1)
    assert np.isclose(img.max(), 10.0)
    # flux varies along the dispersion (X) axis -- not a flat continuum
    col_sums = img.sum(axis=0)
    assert not np.allclose(col_sums, col_sums[0])
    # most flux concentrated near the 4500 A peak
    waves = _linear_wcs(nx, extent=extent).spectral.pixel_to_world(np.arange(nx)).to(u.AA).value
    assert col_sums.argmax() == np.abs(waves - 4500).argmin()


def test_add_source_spectrum_cropped_to_extent():
    nx, ny = 120, 40
    extent = (3000, 6000)
    wcs = _linear_wcs(nx, extent=extent)
    waves = wcs.spectral.pixel_to_world(np.arange(nx)).to(u.AA).value
    wave = np.linspace(4000, 5000, 200) * u.AA
    flux = np.ones(200) * u.count
    spectrum = Spectrum(flux=flux, spectral_axis=wave)
    img = SynthImage(nx=nx, ny=ny, wcs=wcs).add_source(spectrum=spectrum).to_array()
    # flux outside the spectrum's range is zero
    assert np.all(img[:, waves < 4000] == 0)
    assert np.all(img[:, waves > 5000] == 0)
    # flux present within the spectrum's range
    assert img[:, (waves > 4400) & (waves < 4600)].sum() > 0


def test_add_source_spectrum_outside_extent_is_zero():
    nx, ny = 80, 30
    extent = (3000, 6000)
    # spectrum lies entirely outside the image's wavelength range -> no overlap
    wave = np.linspace(8000, 9000, 100) * u.AA
    flux = np.ones(100) * u.count
    spectrum = Spectrum(flux=flux, spectral_axis=wave)
    img = (
        SynthImage(nx=nx, ny=ny, extent=extent)
        .add_source(profile=models.Gaussian1D(amplitude=10, stddev=5), spectrum=spectrum)
        .to_array()
    )
    # nothing overlaps, so the normalized flux (and the source) is all zero
    assert np.all(img == 0)


def test_add_source_spectrum_builds_wcs_from_extent():
    extent = (3000, 6000)
    wave = np.linspace(extent[0], extent[1], 100) * u.AA
    flux = np.ones(100) * u.count
    spectrum = Spectrum(flux=flux, spectral_axis=wave)
    img = SynthImage(nx=80, ny=30, extent=extent).add_source(spectrum=spectrum)
    arr = img.to_array()
    assert arr.shape == (30, 80)
    assert arr.max() > 0
    # the built WCS is carried through to the rendered CCDData
    assert img.to_ccddata().wcs is not None


def test_add_source_spectrum_requires_wcs_or_extent():
    wave = np.linspace(3000, 6000, 50) * u.AA
    flux = np.ones(50) * u.count
    spectrum = Spectrum(flux=flux, spectral_axis=wave)
    img = SynthImage(nx=60, ny=20, extent=None).add_source(spectrum=spectrum)
    with pytest.raises(ValueError, match="Must specify either a wavelength extent or a WCS"):
        img.to_array()


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


@pytest.mark.remote_data
@pytest.mark.filterwarnings("ignore:No observer defined on WCS")
def test_add_skylines_matches_add_arcs():
    wcs = _linear_wcs(300, extent=(6500, 9500))
    sky = SynthImage(nx=300, ny=100, wcs=wcs, seed=1).add_skylines("OH_GMOS")
    arc = SynthImage(nx=300, ny=100, wcs=wcs, seed=1).add_arcs("OH_GMOS")
    assert np.array_equal(sky.to_array(), arc.to_array())


def test_read_noise_changes_image_and_is_reproducible():
    base = SynthImage(nx=40, ny=40, seed=42).add_background(100)
    noisy = base.add_read_noise(5)
    a = noisy.to_array()
    b = noisy.to_array()
    assert np.array_equal(a, b)              # same seed -> reproducible
    assert not np.allclose(a, 100)           # read noise actually applied
    assert abs(a.std() - 5) < 1.0            # ~ sigma


def test_poisson_noise_reproducible_with_seed():
    img = SynthImage(nx=40, ny=40, seed=7).add_background(100).add_poisson_noise()
    assert np.array_equal(img.to_array(), img.to_array())


def test_different_seed_differs():
    a = SynthImage(nx=40, ny=40, seed=1).add_background(100).add_poisson_noise().to_array()
    b = SynthImage(nx=40, ny=40, seed=2).add_background(100).add_poisson_noise().to_array()
    assert not np.array_equal(a, b)


def test_unseeded_runs():
    arr = SynthImage(nx=20, ny=20).add_background(100).add_poisson_noise().to_array()
    assert arr.shape == (20, 20)


def test_rdnoise_alias_equivalent():
    a = SynthImage(nx=20, ny=20, seed=3).add_background(50).add_read_noise(4).to_array()
    b = SynthImage(nx=20, ny=20, seed=3).add_background(50).add_rdnoise(4).to_array()
    assert np.array_equal(a, b)


def test_poisson_then_read_noise_both_applied():
    # both noise stages applied; result differs from poisson-only with same seed
    poisson_only = SynthImage(nx=40, ny=40, seed=5).add_background(100).add_poisson_noise()
    both = poisson_only.add_read_noise(5)
    assert not np.array_equal(poisson_only.to_array(), both.to_array())


def test_to_spectrum_no_wcs():
    sp = SynthImage(nx=60, ny=20).add_background(5).to_spectrum()
    assert isinstance(sp, Spectrum)
    assert sp.flux.shape == (20, 60)
    assert sp.flux.unit == u.count


@pytest.mark.remote_data
@pytest.mark.filterwarnings("ignore:No observer defined on WCS")
def test_to_spectrum_with_wcs():
    wcs = _linear_wcs(300)
    sp = SynthImage(nx=300, ny=100, wcs=wcs).add_arcs(["HeI"]).to_spectrum()
    assert isinstance(sp, Spectrum)
    assert sp.flux.shape == (100, 300)
    assert sp.spectral_axis is not None


def test_make_2d_trace_image_deprecated_and_matches():
    with pytest.warns(AstropyDeprecationWarning):
        ccd = make_2d_trace_image(add_noise=False)
    expected = (
        SynthImage(nx=3000, ny=1000)
        .add_background(5)
        .add_source(profile=models.Moffat1D(amplitude=10, alpha=0.1))
        .to_array()
    )
    assert isinstance(ccd, CCDData)
    assert np.array_equal(ccd.data, expected)


@pytest.mark.remote_data
@pytest.mark.filterwarnings("ignore:No observer defined on WCS")
def test_make_2d_arc_image_deprecated_and_matches():
    with pytest.warns(AstropyDeprecationWarning):
        ccd = make_2d_arc_image(nx=300, ny=100, add_noise=False)
    expected = (
        SynthImage(nx=300, ny=100, extent=(3500, 7000))
        .add_background(5)
        .add_arcs(("HeI",))
        .to_array()
    )
    assert np.array_equal(ccd.data, expected)


@pytest.mark.remote_data
@pytest.mark.filterwarnings("ignore:No observer defined on WCS")
def test_make_2d_spec_image_deprecated_and_matches():
    with pytest.warns(AstropyDeprecationWarning):
        ccd = make_2d_spec_image(nx=300, ny=100, add_noise=False)
    expected = (
        SynthImage(nx=300, ny=100, extent=(6500, 9500))
        .add_background(5)
        .add_arcs(("OH_GMOS",))
        .add_source(profile=models.Moffat1D(amplitude=10, alpha=0.1))
        .to_array()
    )
    assert isinstance(ccd, CCDData)
    assert ccd.data.shape == (100, 300)
    assert np.array_equal(ccd.data, expected)
