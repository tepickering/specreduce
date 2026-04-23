import matplotlib.pyplot as plt
import numpy as np
import pytest

from specreduce.tilt_correction import TiltCorrection
from specreduce.tracing import FlatTrace


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
