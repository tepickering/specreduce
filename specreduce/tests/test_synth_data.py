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
