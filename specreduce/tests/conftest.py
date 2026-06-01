import numpy as np
import pytest
from astropy.modeling import models
from astropy.nddata import StdDevUncertainty
from astropy.wcs import WCS

from specreduce.tilt_correction import TiltCorrection
from specreduce.utils.synth_data import SynthImage


# Arc frame creation code taken from Tim Pickering's example notebook
@pytest.fixture
def mk_arc_frames():
    blue_channel_header = {
        "CTYPE1": "AWAV-GRA",  # Grating dispersion function with air wavelengths
        "CUNIT1": "Angstrom",  # Dispersion units
        "CRPIX1": 1344 // 2,  # Reference pixel [pix]
        "CRVAL1": 5410,  # Reference value [Angstrom]
        "CDELT1": 1.19 * 2,  # Linear dispersion [Angstrom/pix]
        "PV1_0": 5.0e5,  # Grating density [1/m]
        "PV1_1": 1,  # Diffraction order
        "PV1_2": 8.05,  # Incident angle [deg]
        "CTYPE2": "PIXEL",  # Spatial detector coordinates
        "CUNIT2": "pix",  # Spatial units
        "CRPIX2": 1,  # Reference pixel
        "CRVAL2": 0,  # Reference value
        "CDELT2": 1,  # Spatial units per pixel
    }
    blue_channel_wcs = WCS(header=blue_channel_header)
    tilt_mod = models.Legendre1D(degree=2, c0=25, c1=40, c2=25)

    arcs = []
    for ll in (["HeI", "NeI", "XeI"], ["ArI"]):
        arc = (
            SynthImage(nx=512, ny=128, wcs=blue_channel_wcs)
            .add_background(5)
            .add_arcs(
                linelists=ll,
                line_fwhm=3,
                tilt_func=tilt_mod,
                amplitude_scale=1e-2,
            )
            .add_poisson_noise()
            .to_ccddata()
        )
        arc.wcs = None
        arc.uncertainty = StdDevUncertainty(np.full_like(arc.data, 5))
        arcs.append(arc)
    return arcs


@pytest.fixture
def mk_default_tc(mk_arc_frames):
    tc = TiltCorrection(
        arc_frames=mk_arc_frames,
        cdisp_ref_pixel=64,
        n_cdisp_samples=8,
    )
    return tc
