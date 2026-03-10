__all__ = ["TiltSolution"]

from functools import cached_property
from typing import Sequence, Literal

import astropy.units as u
import gwcs
import numpy as np
from astropy.modeling import models, Model
from astropy.modeling.models import Identity, Mapping
from astropy.nddata import NDData
from gwcs import coordinate_frames
from numpy import ndarray

from specreduce.core import _ImageParser


def diff_poly2d_x(model: models.Polynomial2D) -> models.Polynomial2D:
    """Compute the partial derivative of a 2D polynomial model with respect to x.

    Generates a new 2D polynomial model representing the derivative of the input
    model in the x-direction. The coefficients of the resulting model are calculated
    by multiplying the coefficients from the input model by their respective x
    index and reducing the order in the x-dimension.

    Parameters
    ----------
    model
        An `astropy.modeling.models.Polynomial2D` model.

    Returns
    -------
    models.Polynomial2D
        A new 2D polynomial model representing the derivative of the input model
        with respect to x. The degree of the resulting model will be decreased
        by 1 in the x-dimension.
    """
    coeffs = {}
    for n in model.param_names:
        ix, iy = int(n[1]), int(n[3])
        if ix > 0:
            coeffs[f"c{ix-1}_{iy}"] = ix * getattr(model, n).value
    return models.Polynomial2D(model.degree - 1, **coeffs)


class TiltSolution:
    def __init__(self, solution: Model, disp_axis: int = 1):
        """A solution for 2D spectral tilt correction.

        This class encapsulates the polynomial transformation from a tilt-corrected
        (rectified) coordinate space to detector space. It provides methods for
        coordinate transformation, flux-conserving resampling, and export to a GWCS
        object.

        Parameters
        ----------
        solution
            An Astropy compound model representing the transformation from
            tilt-corrected space to detector space along the dispersion axis.
        disp_axis
            The index of the image's dispersion axis, by default 1.
        """
        self._shift: Model = solution[:2]
        self._r2d: Model = solution
        self._r2d_dxdx: None | Model = None
        self.disp_axis = disp_axis

    @property
    def cor2det(self):
        """Transformation from tilt-corrected to detector space along the dispersion axis."""
        return self._r2d

    @cor2det.setter
    def cor2det(self, value):
        self._r2d = value
        if "cor2det_derivative" in self.__dict__:
            del self.cor2det_derivative
        if "gwcs" in self.__dict__:
            del self.gwcs

    @cached_property
    def cor2det_derivative(self):
        """Dispersion-axis derivative of the tilt-corrected to detector space transformation."""
        self._calculate_derivative()
        return self._r2d_dxdx

    @cached_property
    def gwcs(self) -> gwcs.wcs.WCS:
        """GWCS object defining the 2D rectified-to-detector coordinate mapping.

        The forward transform maps ``(disp_rectified, cdisp)`` to
        ``(disp_detector, cdisp)``, where the cross-dispersion coordinate
        passes through unchanged.
        """
        rectified_frame = coordinate_frames.CoordinateFrame(
            2, ("PIXEL", "PIXEL"), (0, 1),
            axes_names=("disp", "cdisp"), unit=[u.pix, u.pix],
            name="rectified",
        )
        detector_frame = coordinate_frames.CoordinateFrame(
            2, ("PIXEL", "PIXEL"), (0, 1),
            axes_names=("disp", "cdisp"), unit=[u.pix, u.pix],
            name="detector",
        )
        # self._r2d maps (disp, cdisp) -> disp_det (2 inputs, 1 output).
        # Build a 2D->2D transform: (disp, cdisp) -> (disp_det, cdisp).
        full_transform = Mapping((0, 1, 1)) | (self._r2d & Identity(1))
        pipeline = [(rectified_frame, full_transform), (detector_frame, None)]
        return gwcs.wcs.WCS(pipeline)

    def _calculate_derivative(self):
        """Calculate the derivative for the tilt-corrected space -> detector space transform."""
        self._r2d_dxdx = self._shift | diff_poly2d_x(self._r2d[-1])

    def rec_to_det(self, disp: ndarray, cdisp: ndarray) -> tuple[ndarray, ndarray]:
        """Transform coordinates from the 2D tilt-corrected space to 2D detector space.

        Parameters
        ----------
        disp
            The dispersion-axis coordinates to be transformed.
        cdisp
            The cross-dispersion coordinates, returned as is.

        Returns
        -------
        tuple of (ndarray, ndarray)
            A tuple containing the transformed dispersion-axis coordinates as the first element
            and the original cross-dispersion-axis coordinates as the second element.
        """
        return self._r2d(disp, cdisp), cdisp

    def resample(
        self,
        flux: NDData,
        nbins: int | None = None,
        bounds: tuple[float, float] | None = None,
        bin_edges: None | Sequence[float] = None,
        mask_treatment: Literal[
            "apply",
            "ignore",
            "propagate",
            "zero_fill",
            "nan_fill",
            "apply_mask_only",
            "apply_nan_only",
        ] = "apply",
    ):
        """Resample a 2D spectrum from the detector space to a tilt-corrected space.

        Resample a 2D spectrum from the detector space to a tilt-corrected space where the
        wavelength is constant along the cross-dispersion axis. The grid edges are based on the
        specified number of bins, bounds, or bin edges. The resampling is exact and conserves
        flux (as long as the tilt-corrected space covers the whole detector space.)

        Parameters
        ----------
        flux
            2D spectrum as an NDData instance. The dispersion and cross-dispersion axis alignment
            should be the same as in the arc frames.
        nbins
            Number of bins in the tilt-corrected space. If None, the number of bins will be set
            to the number of columns in the `flux` input image.
        bounds
            Tuple specifying the start and end coordinates for the tilt-corrected space along the
            x-axis. If None, the bounds default to (0, number of columns in ``flux``).
        bin_edges
            Explicitly provided edges of the bins in the tilt-corrected space. If None, bin
            edges are automatically calculated as a uniform grid based on ``nbins`` and
            ``bounds``.
        mask_treatment
             Specifies how to handle masked or non-finite values in the input image.
             The accepted values are:

             - ``apply``: The image remains unchanged, and any existing mask is combined with a mask
             derived from non-finite values.
             - ``ignore``: The image remains unchanged, and any existing mask is dropped.
             - ``propagate``: The image remains unchanged, and any masked or non-finite pixel
             causes the mask to extend across the entire cross-dispersion axis.
             - ``zero_fill``: Pixels that are either masked or non-finite are replaced with 0.0,
             and the mask is dropped.
             - ``nan_fill``:  Pixels that are either masked or non-finite are replaced with nan,
             and the mask is dropped.
             - ``apply_mask_only``: The  image and mask are left unmodified.
             - ``apply_nan_only``: The  image is left unmodified, the old mask is dropped,
             and a new mask is created based on non-finite values.

        Returns
        -------
        NDData
            NDData instance containing the flux values resampled into the uniform grid
            defined by ``nbins``, ``bounds``, or ``bin_edges``.
        """

        # TODO: In the future, we want to make sure that we don't copy the data unless absolutely
        # necessary.
        ip = _ImageParser()
        im = ip._parse_image(flux, disp_axis=self.disp_axis, mask_treatment=mask_treatment)
        flux = im.flux.value

        ny, nx = flux.data.shape
        ypix = np.arange(ny)
        nbins = nx if nbins is None else nbins
        l1, l2 = bounds if bounds is not None else (0, nx)

        bin_edges_rec = bin_edges if bin_edges is not None else np.linspace(l1, l2, num=nbins + 1)
        bin_edges_det = np.clip(self._r2d(*np.meshgrid(bin_edges_rec, ypix)), 0, nx - 1e-12)
        bin_edge_ix = np.floor(bin_edges_det).astype(int)
        bin_edge_w = bin_edges_det - bin_edge_ix

        resampled_flux = np.zeros((ny, nbins))
        weights = np.zeros((ny, nx))

        # Calculate the derivative of the tilt-corrected space -> detector space transformation  with
        # respect to the detector coordinate (dx_rec / dx_det). This is needed for flux
        # conservation, as it represents how the pixel width changes.
        dtdx = self.cor2det_derivative(*np.meshgrid(np.arange(nx), np.arange(ny)))

        # Calculate a normalization factor 'n' for flux conservation. This factor accounts for the
        # change in pixel size due to the distortion, and ensures that the total flux in each row
        # is conserved after tilt-correction
        n = flux.sum(1) / (dtdx * flux).sum(1)

        ixs = np.tile(np.arange(flux.shape[1]), (flux.shape[0], 1))
        ys = np.arange(flux.shape[0])

        for i in range(nbins):
            # Get the detector pixel indices (left and right edges) for the current tilt-corrected bin.
            i1, i2 = bin_edge_ix[:, i : i + 2].T

            # Create a mask 'm' where the left and right detector pixel edges are the same.
            # This means the entire tilt-corrected bin falls within a single detector pixel.
            m = i1 == i2

            # For rows where the tilt-corrected bin falls within a single detector pixel,
            # the tilt-corrected flux is the detector flux in that pixel, scaled by the width of the
            # tilt-corrected bin in detector coordinates and the derivative dtdx.
            if m.any():
                resampled_flux[:, i] = (
                    (bin_edges_det[:, i + 1] - bin_edges_det[:, i]) * flux[ys, i1] * dtdx[ys, i1]
                )

            # For rows where the tilt-corrected bin spans multiple detector pixels, calculate the
            # tilt-corrected flux as a weighted sum of the detector flux, multiplied by dtdx,
            # within the span [imin, imax].
            if not m.all():
                imin, imax = i1.min(), i2.max() + 1
                ixc = ixs[:, imin:imax]
                w = weights[:, imin:imax]
                w[:] = 0.0
                w[(ixc > i1[:, None]) & (ixc < i2[:, None])] = 1
                w[ys, i1 - imin] = 1.0 - bin_edge_w[:, i]
                w[ys, i2 - imin] = bin_edge_w[:, i + 1]
                resampled_flux[~m, i] = (flux[~m, imin:imax] * dtdx[~m, imin:imax] * w[~m]).sum(1)

        # Apply the normalization factor to conserve flux
        resampled_flux *= n[:, None]
        if self.disp_axis == 0:
            resampled_flux = resampled_flux.T

        return NDData(resampled_flux * im.unit)
