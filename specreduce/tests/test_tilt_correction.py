import matplotlib.pyplot as plt
import numpy as np
import pytest
from astropy.modeling import models
from astropy.nddata import StdDevUncertainty
from astropy.wcs import WCS

from specreduce.tilt_correction import TiltCorrection
from specreduce.tilt_solution import diff_poly2d_x
from specreduce.tracing import FlatTrace
from specreduce.utils.synth_data import make_2d_arc_image


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
        arc = make_2d_arc_image(
            nx=512,
            ny=128,
            linelists=ll,
            wcs=blue_channel_wcs,
            line_fwhm=3,
            tilt_func=tilt_mod,
            amplitude_scale=1e-2,
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


def test_diff_poly2d_x_valid_derivative():
    model = models.Polynomial2D(degree=2, c0_0=1, c1_0=2, c2_0=3, c0_1=4, c1_1=5, c0_2=6)
    derivative = diff_poly2d_x(model)
    assert derivative.degree == 1
    assert derivative.c0_0 == 2
    assert derivative.c1_0 == 6
    assert derivative.c0_1 == 5


def test_init_trace(mk_arc_frames):
    arcs = mk_arc_frames
    trace = FlatTrace(arcs[0], arcs[0].shape[0] // 2)
    tc = TiltCorrection(arc_frames=arcs, trace=trace)
    assert tc.ref_pixel == (arcs[0].shape[0] // 2, arcs[0].shape[1] // 2)


def test_init_default_params(mk_arc_frames):
    arcs = mk_arc_frames
    tc = TiltCorrection(arcs, cdisp_ref_pixel=64, disp_ref_pixel=256)
    assert tc.ref_pixel == (64, 256)
    assert tc.disp_axis == 1
    assert tc.mask_treatment == "apply"
    assert len(tc.arc_frames) == 2

    tc = TiltCorrection(arcs, cdisp_ref_pixel=64)
    assert tc.ref_pixel == (64, arcs[0].shape[1] // 2)

    tc = TiltCorrection(arcs[0], cdisp_ref_pixel=64)
    assert tc.ref_pixel == (64, arcs[0].shape[1] // 2)

    with pytest.raises(ValueError, match="cdisp_ref_position must be provided"):
        TiltCorrection(arcs)

    tc = TiltCorrection(arcs[0], cdisp_ref_pixel=64, cdisp_samples=[10, 20, 30, 40, 50])
    np.testing.assert_array_equal(tc.cd_samples, np.array([10, 20, 30, 40, 50]))


def test_find_lines(mk_default_tc):
    tc = mk_default_tc
    tc.find_arc_lines(3.0, 5.0)
    np.testing.assert_array_equal(tc.cd_samples, np.array([14, 28, 43, 57, 71, 85, 100, 114]))


def test_fit(mk_default_tc):
    tc = mk_default_tc
    tc.find_arc_lines(3.0, 5.0)
    tc.fit(4)


def test_plot_fit_quality(mk_default_tc):
    tc = mk_default_tc
    tc.find_arc_lines(3.0, 5.0)
    tc.fit(4)
    tc.plot_fit_quality()


def test_plot_wavelength_contours(mk_default_tc):
    tc = mk_default_tc
    tc.find_arc_lines(3.0, 5.0)
    tc.fit(4)
    tc.plot_wavelength_contours()


def test_resample(mk_default_tc, mk_arc_frames):
    arcs = mk_arc_frames
    tc = mk_default_tc
    tc.find_arc_lines(3.0, 5.0)
    tc.fit(4)
    tc.solution.resample(arcs[0])


def test_tilt_solution_gwcs(mk_default_tc):
    tc = mk_default_tc
    tc.find_arc_lines(3.0, 5.0)
    tc.fit(4)
    ts = tc.solution
    wcs = ts.gwcs

    # Verify frame names
    assert wcs.input_frame.name == "rectified"
    assert wcs.output_frame.name == "detector"

    # Verify forward transform matches rec_to_det
    disp_arr = np.array([100.0, 200.0, 300.0])
    cdisp_arr = np.array([30.0, 60.0, 90.0])
    det_x_expected, det_y_expected = ts.corr_to_det(disp_arr, cdisp_arr)
    det_x_gwcs, det_y_gwcs = wcs(disp_arr, cdisp_arr)
    np.testing.assert_allclose(det_x_gwcs, det_x_expected)
    np.testing.assert_allclose(det_y_gwcs, det_y_expected)

    # Verify cdisp passes through unchanged
    np.testing.assert_allclose(det_y_gwcs, cdisp_arr)


def test_tilt_solution_gwcs_cache_invalidation(mk_default_tc):
    tc = mk_default_tc
    tc.find_arc_lines(3.0, 5.0)
    tc.fit(4)
    ts = tc.solution

    # Access gwcs to populate the cache
    _ = ts.gwcs
    assert "gwcs" in ts.__dict__

    # Setting cor2det should invalidate the gwcs cache
    ts.c2d = ts.c2d
    assert "gwcs" not in ts.__dict__


def test_det_to_corr_round_trip(mk_default_tc):
    tc = mk_default_tc
    tc.find_arc_lines(3.0, 5.0)
    tc.fit(4)
    ts = tc.solution

    # Test round-trip: corr_to_det then det_to_corr should recover original coordinates
    disp_cor = np.array([100.0, 200.0, 300.0, 400.0])
    cdisp = np.array([20.0, 40.0, 80.0, 100.0])

    disp_det, cdisp_out = ts.corr_to_det(disp_cor, cdisp)
    disp_cor_recovered, cdisp_out2 = ts.det_to_corr(disp_det, cdisp_out)

    np.testing.assert_allclose(disp_cor_recovered, disp_cor, atol=0.01)
    np.testing.assert_allclose(cdisp_out2, cdisp)


def test_d2c_cache_invalidation(mk_default_tc):
    tc = mk_default_tc
    tc.find_arc_lines(3.0, 5.0)
    tc.fit(4)
    ts = tc.solution

    # Access d2c to populate the cache
    _ = ts.d2c
    assert "d2c" in ts.__dict__

    # Setting c2d should invalidate the d2c cache
    ts.c2d = ts.c2d
    assert "d2c" not in ts.__dict__


def test_gwcs_inverse(mk_default_tc):
    tc = mk_default_tc
    tc.find_arc_lines(3.0, 5.0)
    tc.fit(4)
    ts = tc.solution

    disp_cor = np.array([100.0, 200.0, 300.0])
    cdisp = np.array([30.0, 60.0, 90.0])

    # Forward: rectified -> detector via GWCS
    disp_det, cdisp_det = ts.gwcs(disp_cor, cdisp)

    # Inverse: detector -> rectified via GWCS
    disp_cor_inv, cdisp_cor_inv = ts.gwcs.invert(disp_det, cdisp_det)

    # Should match det_to_corr
    disp_rec_direct, cdisp_direct = ts.det_to_corr(disp_det, cdisp_det)
    np.testing.assert_allclose(disp_cor_inv, disp_rec_direct)
    np.testing.assert_allclose(cdisp_cor_inv, cdisp_direct)


def test_refine_fit_before_fit(mk_default_tc):
    tc = mk_default_tc
    tc.find_arc_lines(3.0, 5.0)
    with pytest.raises(ValueError, match="solution must be calculated"):
        tc.refine_fit()


def test_match_lines_before_fit(mk_default_tc):
    tc = mk_default_tc
    tc.find_arc_lines(3.0, 5.0)
    with pytest.raises(ValueError, match="solution must be calculated"):
        tc.match_lines()


def test_plot_wavelength_contours_options(mk_default_tc):
    tc = mk_default_tc
    tc.find_arc_lines(3.0, 5.0)
    tc.fit(4)

    # Test with line_args and a pre-created ax
    fig, ax = plt.subplots()
    tc.plot_wavelength_contours(ax=ax, line_args={"c": "red"})

    # Test with explicit disp_values
    disp_values = np.array([50.0, 100.0, 200.0, 300.0, 400.0])
    tc.plot_wavelength_contours(disp_values=disp_values)
    plt.close("all")


def test_plot_fit_quality_with_rlim(mk_default_tc):
    tc = mk_default_tc
    tc.find_arc_lines(3.0, 5.0)
    tc.fit(4)
    tc.plot_fit_quality(rlim=(-1, 1))
    plt.close("all")
