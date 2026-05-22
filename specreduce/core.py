# Licensed under a 3-clause BSD style license - see LICENSE.rst

import inspect
from copy import deepcopy
from dataclasses import dataclass
from typing import Literal

import numpy as np
from astropy import units as u
from astropy.nddata import VarianceUncertainty, NDData

from specutils import Spectrum

__all__ = ["SpecreduceOperation", "parse_image"]

MaskingOption = Literal[
    "apply", "ignore", "propagate", "zero_fill", "nan_fill", "apply_mask_only", "apply_nan_only"
]

ImageLike = np.ndarray | NDData | u.Quantity

# The '_valid_mask_treatment_methods' tuples in the Background, Trace, and
# Extract classes are subsets of these implemented methods.
_IMPLEMENTED_MASK_TREATMENTS = (
    "apply",
    "ignore",
    "propagate",
    "zero_fill",
    "nan_fill",
    "apply_mask_only",
    "apply_nan_only",
)


def parse_image(
    image: ImageLike, disp_axis: int = 1, mask_treatment: MaskingOption = "apply"
) -> Spectrum:
    """
    Convert all accepted image types to a consistently formatted Spectrum object.

    Fills any and all of uncertainty, mask, units, and spectral axis that are
    missing in the provided image with generic values. Accepted image types are:

        - `~specutils.Spectrum` (preferred)
        - `~astropy.nddata.ccddata.CCDData`
        - `~astropy.nddata.ndddata.NDDData`
        - `~astropy.units.quantity.Quantity`
        - `~numpy.ndarray`

    Parameters
    ----------
    image
        Input image from which data is extracted. This can be a 2D numpy
        array, Quantity, or an NDData object.
    disp_axis
        The index of the image's dispersion axis. Should not be changed until
        operations can handle variable image orientations.
    mask_treatment
        Specifies how to handle masked or non-finite values in the input image.
        The accepted values are:

        - ``apply``: The image remains unchanged, and any existing mask is combined\
            with a mask derived from non-finite values.
        - ``ignore``: The image remains unchanged, and any existing mask is dropped.
        - ``propagate``: The image remains unchanged, and any masked or non-finite pixel\
            causes the mask to extend across the entire cross-dispersion axis.
        - ``zero_fill``: Pixels that are either masked or non-finite are replaced with 0.0,\
            and the mask is dropped.
        - ``nan_fill``:  Pixels that are either masked or non-finite are replaced with nan,\
            and the mask is dropped.
        - ``apply_mask_only``: The  image and mask are left unmodified.
        - ``apply_nan_only``: The  image is left unmodified, the old mask is dropped, and a\
            new mask is created based on non-finite values.

    Returns
    -------
    Spectrum
    """
    if isinstance(image, u.quantity.Quantity):
        img = image.value
    elif isinstance(image, np.ndarray):
        img = image
    else:  # NDData, including CCDData and Spectrum
        img = image.data

    mask = getattr(image, "mask", None)
    crossdisp_axis = (disp_axis + 1) % 2

    # A mask is built from any non-finite image data, combined with any existing
    # mask. If a fill value is chosen the image data is also modified. The
    # returned Spectrum always carries a mask, even if everything is finite.
    img, mask = _apply_mask_treatment(
        image=img, mask=mask, mask_treatment=mask_treatment, crossdisp_axis=crossdisp_axis
    )

    # mask and uncertainty default to None on a Spectrum, so test for both
    # absence of the attribute and presence-with-None.
    if hasattr(image, "uncertainty"):
        uncertainty = image.uncertainty
    else:
        uncertainty = VarianceUncertainty(np.ones(img.shape))

    unit = getattr(image, "unit", u.Unit("DN"))

    spectral_axis = getattr(image, "spectral_axis", np.arange(img.shape[disp_axis]) * u.pix)

    return Spectrum(
        img * unit, spectral_axis=spectral_axis, uncertainty=uncertainty, mask=mask,
        spectral_axis_index=img.ndim - 1
    )


def _apply_mask_treatment(
    image: np.ndarray,
    mask: np.ndarray | None = None,
    mask_treatment: MaskingOption = "apply",
    crossdisp_axis: int = 1,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Handle the treatment of masked and non-finite data.

    Parameters
    ----------
    image
        The input image data array that may contain non-finite values.
    mask
        An optional Boolean mask array. Non-finite values in the image will be
        added to this mask.
    mask_treatment
        Specifies how to handle masked or non-finite values in the input image.
    crossdisp_axis
        Cross-dispersion axis, used for the ``propagate`` treatment.
    """
    if mask_treatment not in _IMPLEMENTED_MASK_TREATMENTS:
        raise ValueError(
            f"'mask_treatment' must be one of {_IMPLEMENTED_MASK_TREATMENTS}"
        )

    if mask is not None and (mask.dtype not in (bool, int)):
        raise ValueError("'mask' must be a boolean or integer array.")

    match mask_treatment:
        case "apply":
            mask = mask | (~np.isfinite(image)) if mask is not None else ~np.isfinite(image)
        case "ignore":
            mask = np.zeros(image.shape, dtype=bool)
        case "propagate":
            if mask is None:
                mask = ~np.isfinite(image)
            else:
                mask = mask | (~np.isfinite(image))
            mask[:] = mask.any(axis=crossdisp_axis, keepdims=True)
        case "zero_fill" | "nan_fill":
            mask = mask | (~np.isfinite(image)) if mask is not None else ~np.isfinite(image)
            image = deepcopy(image)
            if mask_treatment == "zero_fill":
                image[mask] = 0.0
            else:
                image[mask] = np.nan
            mask[:] = False
        case "apply_nan_only":
            mask = ~np.isfinite(image)
        case "apply_mask_only":
            mask = mask.copy() if mask is not None else np.zeros(image.shape, dtype=bool)

    if mask.all():
        raise ValueError("Image is fully masked. Check for invalid values.")

    return image, mask


@dataclass
class SpecreduceOperation:
    """
    An operation to perform as part of a spectroscopic reduction pipeline.

    This class primarily exists to define the basic API for operations:
    parameters for the operation are provided at object creation,
    and then the operation object is called with the data objects required for
    the operation, which then *return* the data objects resulting from the
    operation.
    """

    def __call__(self):
        raise NotImplementedError("__call__ on a SpecreduceOperation needs to " "be overridden")

    @classmethod
    def as_function(cls, *args, **kwargs):
        """
        Run this operation as a function.  Syntactic sugar for e.g.,
        ``Operation.as_function(arg1, arg2, keyword=value)`` maps to
        ``Operation(arg2, keyword=value)(arg1)`` (if the ``__call__`` of
        ``Operation`` has only one argument)
        """
        argspec = inspect.getargs(cls.__call__.__code__)
        if argspec.varargs:
            raise NotImplementedError(
                "There is not a way to determine the "
                "number of inputs of a *args style "
                "operation"
            )
        ninputs = len(argspec.args) - 1

        callargs = args[:ninputs]
        noncallargs = args[ninputs:]
        op = cls(*noncallargs, **kwargs)
        return op(*callargs)
