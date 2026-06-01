import numpy as np
import pytest
from astropy import units as u
from astropy.modeling import models
from astropy.nddata import CCDData

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
