__all__ = ["TiltSolution"]

from functools import cached_property
from typing import Sequence

import numpy as np
from astropy.modeling import models
from astropy.nddata import NDData
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
    def __init__(self, solution: models.Model, disp_axis: int = 1):
        self._shift: models.Model = solution[0]
        self._r2d: models.Model = solution
        self._r2d_dxdx: None | models.Model = None
        self.disp_axis = disp_axis

    @property
    def cor2det(self):
        """Tilt-corrected space to detector space transform."""
        return self._r2d

    @cor2det.setter
    def cor2det(self, value):
        self._r2d = value
        if "cor2det_derivative" in self.__dict__:
            del self.cor2det_derivative

    @cached_property
    def cor2det_derivative(self):
        """Tilt-corrected space to detector space transform derivative along the x-axis."""
        self._calculate_derivative()
        return self._r2d_dxdx

    def _calculate_derivative(self):
        """Calculate the derivative for the tilt-corrected space -> detector space transform."""
        self._r2d_dxdx = self._shift | diff_poly2d_x(self._r2d[-1])

    def rec_to_det(self, disp: ndarray, cdisp: ndarray) -> tuple[ndarray, ndarray]:
        """Transform coordinates from the tilt-corrected space to detector space.

        Parameters
        ----------
        disp : ndarray
            The dispersion-axis coordinates to be transformed.
        cdisp : ndarray
            The cross-dispersion coordinates, returned as is.

        Returns
        -------
        tuple of (ndarray, ndarray)
            A tuple containing the transformed dispersion-axis coordinates as the first element
            and the original cross-dispersion-axis coordinates as the second element..
        """
        return self._r2d(disp, cdisp), cdisp

    def resample(
        self,
        flux: NDData,
        nbins: int | None = None,
        bounds: tuple[float, float] | None = None,
        bin_edges: None | Sequence[float] = None,
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
        bound
            Tuple specifying the start and end coordinates for the tilt-corrected space along the
            x-axis. If None, the bounds default to (0, number of columns in ``flux``).
        bin_edges
            Explicitly provided edges of the bins in the tilt-corrected space. If None, bin
            edges are automatically calculated as a uniform grid based on ``nbins`` and
            ``bounds``.

        Returns
        -------
        resampled_flux : ndarray
            NDData instance containing the flux values resampled into the uniform grid
            defined by ``nbins``, ``bounds``, or ``bin_edges``.
        """

        # TODO: In the future, we want to make sure that we don't copy the data unless absolutely
        # necessary.
        ip = _ImageParser()
        im = ip._parse_image(flux, disp_axis=self.disp_axis, mask_treatment=self.mask_treatment)
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
