import numpy as np

import astropy.units as u

from ..calibration_data import (
    BaseAtmosphericExtinction,
    ObservatoryExtinction,
    CustomAtmosphericExtinction,
    AtmosphericTransmission,
    SUPPORTED_EXTINCTION_MODELS
)


def test_base_extinction():
    """
    Test base extinction data class
    """
    ext = BaseAtmosphericExtinction()
    assert(len(ext.extinction_curve > 0))
    assert(len(ext.wavelength) > 0)


def test_supported_models():
    """
    Test loading of supported models
    """
    for observatory in SUPPORTED_EXTINCTION_MODELS:
        ext = ObservatoryExtinction(observatory=observatory)
        assert(len(ext.extinction_curve > 0))
        assert(len(ext.wavelength) > 0)


def test_at_airmass():
    """
    Load default observatory model and test it at different airmasses
    """
    ext = ObservatoryExtinction()
    pre_val = ext.transmission(airmass=1.0)[0]
    post_val = ext.transmission(airmass=1.2)[0]
    assert(post_val < pre_val)
    pre_val = ext.extinction(airmass=1.0)[0]
    post_val = ext.extinction(airmass=1.2)[0]
    assert(post_val > pre_val)


def test_custom_mag_model():
    """
    Test creation of custom model from Quantity arrays
    """
    wave = np.linspace(0.3, 2.0, 50) * u.um
    extinction = u.Magnitude(1. / wave.value, u.MagUnit(u.dimensionless_unscaled))
    ext = CustomAtmosphericExtinction(wavelength=wave, extinction_curve=extinction)
    assert(len(ext.extinction()) > 0)
    assert(len(ext.transmission()) > 0)


def test_custom_raw_mag_model():
    """
    Test creation of custom model from Quantity arrays
    """
    wave = np.linspace(0.3, 2.0, 50)
    extinction = (1. / wave) * u.mag
    ext = CustomAtmosphericExtinction(wavelength=wave, extinction_curve=extinction)
    assert(len(ext.extinction()) > 0)
    assert(len(ext.transmission()) > 0)


def test_custom_linear_model():
    """
    Test creation of custom model from Quantity arrays
    """
    wave = np.linspace(0.3, 2.0, 50)
    extinction = (1. / wave) * u.dimensionless_unscaled
    ext = CustomAtmosphericExtinction(wavelength=wave, extinction_curve=extinction)
    assert(len(ext.extinction()) > 0)
    assert(len(ext.transmission()) > 0)


def test_missing_extinction_unit():
    """
    Test creation of custom model from Quantity arrays
    """
    wave = np.linspace(0.3, 2.0, 50)
    extinction = 1. / wave
    ext = CustomAtmosphericExtinction(wavelength=wave, extinction_curve=extinction)
    assert(len(ext.extinction()) > 0)
    assert(len(ext.transmission()) > 0)


def test_transmission_model():
    ext = AtmosphericTransmission()
    assert(len(ext.extinction()) > 0)
    assert(len(ext.transmission()) > 0)
