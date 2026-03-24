import numpy as np
from astropy.modeling import models

from specreduce.tilt_solution import diff_poly2d_x


def test_diff_poly2d_x_valid_derivative():
    model = models.Polynomial2D(degree=2, c0_0=1, c1_0=2, c2_0=3, c0_1=4, c1_1=5, c0_2=6)
    derivative = diff_poly2d_x(model)
    assert derivative.degree == 1
    assert derivative.c0_0 == 2
    assert derivative.c1_0 == 6
    assert derivative.c0_1 == 5


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


def test_resample(mk_default_tc, mk_arc_frames):
    arcs = mk_arc_frames
    tc = mk_default_tc
    tc.find_arc_lines(3.0, 5.0)
    tc.fit(4)
    tc.solution.resample(arcs[0])
