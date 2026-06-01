# Licensed under a 3-clause BSD style license - see ../../licenses/LICENSE.rst
import warnings
from dataclasses import dataclass, field

import numpy as np
from astropy import units as u
from astropy.modeling import models, Model
from astropy.nddata import CCDData
from astropy.stats import gaussian_fwhm_to_sigma
from astropy.wcs import WCS

from specutils import Spectrum

from specreduce.calibration_data import load_pypeit_calibration_lines

__all__ = [
    "SynthImage",
    "make_2d_trace_image",
    "make_2d_arc_image",
    "make_2d_spec_image",
]

_ALLOWED_TILT = (
    models.Legendre1D,
    models.Chebyshev1D,
    models.Polynomial1D,
    models.Hermite1D,
)


@dataclass(frozen=True)
class _RenderContext:
    """Shared geometry passed to every layer's ``render`` method."""
    nx: int
    ny: int
    xx: np.ndarray
    yy: np.ndarray
    wcs: WCS | None
    disp_axis: int


@dataclass(frozen=True)
class BackgroundLayer:
    """A constant additive background level in counts."""
    level: float

    def render(self, ctx: _RenderContext) -> np.ndarray:
        return np.full((ctx.ny, ctx.nx), float(self.level))


class SynthImage:
    """
    Immutable, composable builder for synthetic 2D spectroscopic images.

    Build an image by chaining ``add_*`` methods, then render it with one of the
    ``to_*`` terminal methods. Each ``add_*`` returns a *new* ``SynthImage``; the
    original is never mutated, so a base configuration can be safely branched.

    Parameters
    ----------
    nx
        Size of the image along the X (dispersion) axis.
    ny
        Size of the image along the Y (spatial) axis.
    wcs
        Optional 2D WCS with a single spectral axis. If not provided and arc
        layers are present, a linear ``WAVE``/``PIXEL`` WCS is built from
        ``extent`` and ``wave_unit``.
    extent
        Beginning and end wavelengths used to build a default WCS when ``wcs``
        is not supplied and arc layers are present.
    wave_unit
        Wavelength unit for the default WCS.
    seed
        Seed for the random number generator used by the noise layers. If
        ``None``, noise is non-deterministic.
    """

    def __init__(
        self,
        nx: int = 3000,
        ny: int = 1000,
        wcs: WCS | None = None,
        extent=(3500, 7000),
        wave_unit: u.Unit = u.Angstrom,
        seed: int | None = None,
    ):
        self.nx = nx
        self.ny = ny
        self._wcs = wcs
        self._extent = extent
        self._wave_unit = wave_unit
        self._seed = seed
        self._layers = ()
        self._poisson = False
        self._read_noise = None

    def _clone(self, **changes) -> "SynthImage":
        new = SynthImage.__new__(SynthImage)
        new.nx = self.nx
        new.ny = self.ny
        new._wcs = self._wcs
        new._extent = self._extent
        new._wave_unit = self._wave_unit
        new._seed = self._seed
        new._layers = self._layers
        new._poisson = self._poisson
        new._read_noise = self._read_noise
        for key, value in changes.items():
            setattr(new, key, value)
        return new

    def add_background(self, level: float) -> "SynthImage":
        """Add a constant background level (counts)."""
        return self._clone(_layers=self._layers + (BackgroundLayer(level),))

    def add_source(
        self,
        profile: Model | None = None,
        trace_center: float | None = None,
        trace_order: int = 3,
        trace_coeffs: dict | None = None,
    ) -> "SynthImage":
        """Add a continuum source with a Chebyshev-traced spatial profile."""
        if profile is None:
            profile = models.Moffat1D(amplitude=10, alpha=0.1)
        layer = SourceLayer(profile, trace_center, trace_order, trace_coeffs)
        return self._clone(_layers=self._layers + (layer,))

    def _resolve_wcs(self):
        has_arc = any(isinstance(layer, ArcLayer) for layer in self._layers)
        if self._wcs is not None:
            wcs = self._wcs
            if wcs.spectral.naxis != 1:
                raise ValueError("Provided WCS must have a spectral axis.")
            if wcs.naxis != 2:
                raise ValueError("WCS must have NAXIS=2 for a 2D image.")
        elif has_arc:
            if self._extent is None:
                raise ValueError("Must specify either a wavelength extent or a WCS.")
            if len(self._extent) != 2:
                raise ValueError("Wavelength extent must be of length 2.")
            if u.get_physical_type(self._wave_unit) != "length":
                raise ValueError("Wavelength unit must be a length unit.")
            wcs = WCS(naxis=2)
            wcs.wcs.ctype[0] = "WAVE"
            wcs.wcs.ctype[1] = "PIXEL"
            wcs.wcs.cunit[0] = self._wave_unit
            wcs.wcs.cunit[1] = u.pixel
            wcs.wcs.crval[0] = self._extent[0]
            wcs.wcs.cdelt[0] = (self._extent[1] - self._extent[0]) / self.nx
            wcs.wcs.crval[1] = 0
            wcs.wcs.cdelt[1] = 1
        else:
            wcs = None

        if wcs is None:
            disp_axis = 1
        else:
            is_spectral = [a["coordinate_type"] == "spectral" for a in wcs.get_axis_types()]
            disp_axis = 0 if is_spectral[0] else 1
        return wcs, disp_axis

    def _render(self):
        wcs, disp_axis = self._resolve_wcs()
        x = np.arange(self.nx)
        y = np.arange(self.ny)
        xx, yy = np.meshgrid(x, y)
        ctx = _RenderContext(self.nx, self.ny, xx, yy, wcs, disp_axis)

        signal = np.zeros((self.ny, self.nx))
        for layer in self._layers:
            signal = signal + layer.render(ctx)

        rng = np.random.default_rng(self._seed)
        if self._poisson:
            from photutils.datasets import apply_poisson_noise
            signal = apply_poisson_noise(signal, seed=rng)
        if self._read_noise is not None:
            signal = signal + rng.normal(0.0, self._read_noise, size=signal.shape)

        return signal, wcs

    def to_array(self) -> np.ndarray:
        """Render and return the image as a plain ``numpy.ndarray`` (counts)."""
        return self._render()[0]

    def to_ccddata(self) -> CCDData:
        """Render and return the image as a `~astropy.nddata.CCDData`."""
        data, wcs = self._render()
        return CCDData(data, unit=u.count, wcs=wcs)

    def to_spectrum(self) -> Spectrum:
        """Render and return the image as a `~specutils.Spectrum`."""
        data, wcs = self._render()
        if wcs is not None:
            return Spectrum(flux=data * u.count, wcs=wcs)
        return Spectrum(flux=data * u.count, spectral_axis_index=data.ndim - 1)


@dataclass(frozen=True)
class SourceLayer:
    """A continuum source whose spatial profile follows a Chebyshev trace.

    The dispersion axis is the X (column) axis, matching the historical
    ``make_2d_trace_image`` behaviour (the trace ignores any WCS).
    """
    profile: Model
    trace_center: float | None = None
    trace_order: int = 3
    trace_coeffs: dict | None = None

    def render(self, ctx: _RenderContext) -> np.ndarray:
        trace_center = ctx.ny / 2 if self.trace_center is None else self.trace_center
        trace_coeffs = (
            {"c0": 0, "c1": 50, "c2": 100}
            if self.trace_coeffs is None
            else self.trace_coeffs
        )
        trace_mod = models.Chebyshev1D(degree=self.trace_order, **trace_coeffs)
        trace = ctx.yy - trace_center + trace_mod(ctx.xx / ctx.nx)
        return self.profile(trace)


@dataclass(frozen=True)
class ArcLayer:
    pass  # implemented in Task 3


def make_2d_trace_image(
    nx: int = 3000,
    ny: int = 1000,
    background: int | float = 5,
    trace_center: int | float | None = None,
    trace_order: int = 3,
    trace_coeffs: None | dict[str, int | float] = None,
    profile: Model = models.Moffat1D(amplitude=10, alpha=0.1),
    add_noise: bool = True
) -> CCDData:
    """
    Create synthetic 2D spectroscopic image with a single source. The spatial (y-axis) position
    of the source along the dispersion (x-axis) direction is modeled using a Chebyshev polynomial.
    The flux units are counts and the noise is modeled as Poisson.

    Parameters
    ----------
    nx
        Size of image in X axis which is assumed to be the dispersion axis

    ny
        Size of image in Y axis which is assumed to be the spatial axis

    background
        Level of constant background in counts

    trace_center
        Zeropoint of the trace. If None, then use center of Y (spatial) axis.

    trace_order
        Order of the Chebyshev polynomial used to model the source's trace

    trace_coeffs
        Dict containing the Chebyshev polynomial coefficients to use in the trace model

    profile
        Model to use for the source's spatial profile

    add_noise
        If True, add Poisson noise to the image

    Returns
    -------
    ccd_im
        `~astropy.nddata.CCDData` instance containing synthetic 2D spectroscopic image
    """
    if trace_coeffs is None:
        trace_coeffs = {'c0': 0, 'c1': 50, 'c2': 100}
    x = np.arange(nx)
    y = np.arange(ny)
    xx, yy = np.meshgrid(x, y)

    if trace_center is None:
        trace_center = ny / 2

    trace_mod = models.Chebyshev1D(degree=trace_order, **trace_coeffs)
    trace = yy - trace_center + trace_mod(xx/nx)
    z = background + profile(trace)

    if add_noise:
        from photutils.datasets import apply_poisson_noise
        trace_image = apply_poisson_noise(z)
    else:
        trace_image = z

    ccd_im = CCDData(trace_image, unit=u.count)

    return ccd_im


def make_2d_arc_image(
    nx: int = 3000,
    ny: int = 1000,
    wcs: WCS | None = None,
    extent=(3500, 7000),
    wave_unit: u.Unit = u.Angstrom,
    wave_air: bool = False,
    background: int | float = 5,
    line_fwhm: float = 5.,
    linelists: list[str] = ('HeI',),
    amplitude_scale: float = 1.,
    tilt_func: Model = models.Legendre1D(degree=0),
    add_noise: bool = True
) -> CCDData:
    """
    Create synthetic 2D spectroscopic image of reference emission lines, e.g. a calibration
    arc lamp. Currently, linelists from ``pypeit`` are supported and are selected by string or
    list of strings that is passed to `~specreduce.calibration_data.load_pypeit_calibration_lines`.
    If a ``wcs`` is not provided, one is created using ``extent`` and ``wave_unit`` with
    dispersion along the X axis.

    Parameters
    ----------
    nx
        Size of image in X axis which is assumed to be the dispersion axis

    ny
        Size of image in Y axis which is assumed to be the spatial axis

    wcs
        2D WCS to apply to the image. Must have a spectral axis defined along with
        appropriate spectral wavelength units.

    extent
        If ``wcs`` is not provided, this defines the beginning and end wavelengths
        of the dispersion axis.

    wave_unit
        If ``wcs`` is not provided, this defines the wavelength units of the
        dispersion axis.

    wave_air
        If True, convert the vacuum wavelengths used by ``pypeit`` to air wavelengths.

    background
        Level of constant background in counts

    line_fwhm
        Gaussian FWHM of the emission lines in pixels

    linelists
        Specification for linelists to load from ``pypeit``

    amplitude_scale
        Scale factor to apply to amplitudes provided in the linelists

    tilt_func
        The tilt function to apply along the cross-dispersion axis to simulate
        tilted or curved emission lines.

    add_noise
        If True, add Poisson noise to the image; requires ``photutils`` to be installed.

    Returns
    -------
    ccd_im
        `~astropy.nddata.CCDData` instance containing synthetic 2D spectroscopic image

    Examples
    --------
    This is an example of modeling a spectrograph whose output is curved in the
    cross-dispersion direction:

    .. plot::
        :include-source:

        import matplotlib.pyplot as plt
        import numpy as np
        from astropy.modeling import models
        import astropy.units as u
        from specreduce.utils.synth_data import make_2d_arc_image

        model_deg2 = models.Legendre1D(degree=2, c0=50, c1=0, c2=100)
        im = make_2d_arc_image(
            linelists=['HeI', 'ArI', 'ArII'],
            line_fwhm=3,
            tilt_func=model_deg2
        )
        fig = plt.figure(figsize=(10, 6))
        plt.imshow(im)

    The FITS WCS standard implements ideal world coordinate functions based on the physics
    of simple dispersers. This is described in detail by Paper III,
    https://www.aanda.org/articles/aa/pdf/2006/05/aa3818-05.pdf. This can be used to model a
    non-linear dispersion relation based on the properties of a spectrograph. This example
    recreates Figure 5 in that paper using a spectrograph with a 450 lines/mm volume phase
    holographic grism. Standard gratings only use the first three ``PV`` terms:

    .. plot::
        :include-source:

        import numpy as np
        import matplotlib.pyplot as plt
        from astropy.wcs import WCS
        import astropy.units as u
        from specreduce.utils.synth_data import make_2d_arc_image

        non_linear_header = {
            'CTYPE1': 'AWAV-GRA',  # Grating dispersion function with air wavelengths
            'CUNIT1': 'Angstrom',  # Dispersion units
            'CRPIX1': 719.8,       # Reference pixel [pix]
            'CRVAL1': 7245.2,      # Reference value [Angstrom]
            'CDELT1': 2.956,       # Linear dispersion [Angstrom/pix]
            'PV1_0': 4.5e5,        # Grating density [1/m]
            'PV1_1': 1,            # Diffraction order
            'PV1_2': 27.0,         # Incident angle [deg]
            'PV1_3': 1.765,        # Reference refraction
            'PV1_4': -1.077e6,     # Refraction derivative [1/m]
            'CTYPE2': 'PIXEL',     # Spatial detector coordinates
            'CUNIT2': 'pix',       # Spatial units
            'CRPIX2': 1,           # Reference pixel
            'CRVAL2': 0,           # Reference value
            'CDELT2': 1            # Spatial units per pixel
        }

        linear_header = {
            'CTYPE1': 'AWAV',  # Grating dispersion function with air wavelengths
            'CUNIT1': 'Angstrom',  # Dispersion units
            'CRPIX1': 719.8,       # Reference pixel [pix]
            'CRVAL1': 7245.2,      # Reference value [Angstrom]
            'CDELT1': 2.956,       # Linear dispersion [Angstrom/pix]
            'CTYPE2': 'PIXEL',     # Spatial detector coordinates
            'CUNIT2': 'pix',       # Spatial units
            'CRPIX2': 1,           # Reference pixel
            'CRVAL2': 0,           # Reference value
            'CDELT2': 1            # Spatial units per pixel
        }

        non_linear_wcs = WCS(non_linear_header)
        linear_wcs = WCS(linear_header)

        # this re-creates Paper III, Figure 5
        pix_array = 200 + np.arange(1400)
        nlin = non_linear_wcs.spectral.pixel_to_world(pix_array)
        lin = linear_wcs.spectral.pixel_to_world(pix_array)
        resid = (nlin - lin).to(u.Angstrom)
        plt.plot(pix_array, resid)
        plt.xlabel("Pixel")
        plt.ylabel("Correction (Angstrom)")
        plt.show()

        nlin_im = make_2d_arc_image(
            nx=600,
            ny=512,
            linelists=['HeI', 'NeI'],
            line_fwhm=3,
            wave_air=True,
            wcs=non_linear_wcs
        )
        lin_im = make_2d_arc_image(
            nx=600,
            ny=512,
            linelists=['HeI', 'NeI'],
            line_fwhm=3,
            wave_air=True,
            wcs=linear_wcs
        )

        # subtracting the linear simulation from the non-linear one shows how the
        # positions of lines diverge between the two cases
        plt.imshow(nlin_im.data - lin_im.data)
        plt.show()

    """
    if wcs is None:
        if extent is None:
            raise ValueError("Must specify either a wavelength extent or a WCS.")
        if len(extent) != 2:
            raise ValueError("Wavelength extent must be of length 2.")
        if u.get_physical_type(wave_unit) != 'length':
            raise ValueError("Wavelength unit must be a length unit.")
        wcs = WCS(naxis=2)
        wcs.wcs.ctype[0] = 'WAVE'
        wcs.wcs.ctype[1] = 'PIXEL'
        wcs.wcs.cunit[0] = wave_unit
        wcs.wcs.cunit[1] = u.pixel
        wcs.wcs.crval[0] = extent[0]
        wcs.wcs.cdelt[0] = (extent[1] - extent[0]) / nx
        wcs.wcs.crval[1] = 0
        wcs.wcs.cdelt[1] = 1
    else:
        if wcs.spectral.naxis != 1:
            raise ValueError("Provided WCS must have a spectral axis.")
        if wcs.naxis != 2:
            raise ValueError("WCS must have NAXIS=2 for a 2D image.")

    x = np.arange(nx)
    y = np.arange(ny)
    xx, yy = np.meshgrid(x, y)

    is_spectral = [a['coordinate_type'] == "spectral" for a in wcs.get_axis_types()]
    if is_spectral[0]:
        disp_axis = 0
    else:
        disp_axis = 1

    if tilt_func is not None:
        if not isinstance(
            tilt_func,
            (models.Legendre1D, models.Chebyshev1D, models.Polynomial1D, models.Hermite1D)
        ):
            raise ValueError(
                "The only tilt functions currently supported are 1D polynomials "
                "from astropy.models."
            )

        if disp_axis == 0:
            xx = xx + tilt_func((yy - ny/2)/ny)
        else:
            yy = yy + tilt_func((xx - nx/2)/nx)

    z = background + np.zeros((ny, nx))

    linelist = load_pypeit_calibration_lines(linelists, wave_air=wave_air)

    if linelist is not None:
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message="No observer defined on WCS.*")
            line_disp_positions = wcs.spectral.world_to_pixel(linelist['wavelength'])

        line_sigma = gaussian_fwhm_to_sigma * line_fwhm
        for line_pos, ampl in zip(line_disp_positions, linelist['amplitude']):
            line_mod = models.Gaussian1D(
                amplitude=ampl * amplitude_scale,
                mean=line_pos,
                stddev=line_sigma
            )
            if disp_axis == 0:
                z += line_mod(xx)
            else:
                z += line_mod(yy)

    if add_noise:
        from photutils.datasets import apply_poisson_noise
        arc_image = apply_poisson_noise(z)
    else:
        arc_image = z

    ccd_im = CCDData(arc_image, unit=u.count, wcs=wcs)

    return ccd_im


def make_2d_spec_image(
    nx: int = 3000,
    ny: int = 1000,
    wcs: WCS | None = None,
    extent=(6500, 9500),
    wave_unit: u.Unit = u.Angstrom,
    wave_air: bool = False,
    background: int | float = 5,
    line_fwhm: float = 5.,
    linelists: list[str] = ('OH_GMOS',),
    amplitude_scale: float = 1.,
    tilt_func: Model = models.Legendre1D(degree=0),
    trace_center: int | float | None = None,
    trace_order: int = 3,
    trace_coeffs: None | dict[str, int | float] = None,
    source_profile: Model = models.Moffat1D(amplitude=10, alpha=0.1),
    add_noise: bool = True
) -> CCDData:
    """
    Make a synthetic 2D spectrum image containing both emission lines and
    a trace of a continuum source.

    Parameters
    ----------
    nx
        Number of pixels in the dispersion direction.

    ny
        Number of pixels in the spatial direction.

    wcs
        2D WCS to apply to the image. Must have a spectral axis defined along with
        appropriate spectral wavelength units.

    extent
        If ``wcs`` is not provided, this defines the beginning and end wavelengths
        of the dispersion axis.

    wave_unit
        If ``wcs`` is not provided, this defines the wavelength units of the
        dispersion axis.

    wave_air
        If True, convert the vacuum wavelengths used by ``pypeit`` to air wavelengths.

    background
        Constant background level in counts.

    line_fwhm
        Gaussian FWHM of the emission lines in pixels

    linelists
        Specification for linelists to load from ``pypeit``

    amplitude_scale
        Scale factor to apply to amplitudes provided in the linelists

    tilt_func
        The tilt function to apply along the cross-dispersion axis to simulate
        tilted or curved emission lines.

    trace_center
        Zeropoint of the trace. If None, then use center of Y (spatial) axis.

    trace_order
        Order of the Chebyshev polynomial used to model the source's trace

    trace_coeffs
        Dict containing the Chebyshev polynomial coefficients to use in the trace model

    source_profile
        Model to use for the source's spatial profile

    add_noise
        If True, add Poisson noise to the image; requires ``photutils`` to be installed.


    Returns
    -------
    ccd_im
        `~astropy.nddata.CCDData` instance containing synthetic 2D spectroscopic image
    """
    if trace_coeffs is None:
        trace_coeffs = {'c0': 0, 'c1': 50, 'c2': 100}

    arc_image = make_2d_arc_image(
        nx=nx,
        ny=ny,
        wcs=wcs,
        extent=extent,
        wave_unit=wave_unit,
        wave_air=wave_air,
        background=0,
        line_fwhm=line_fwhm,
        linelists=linelists,
        amplitude_scale=amplitude_scale,
        tilt_func=tilt_func,
        add_noise=False
    )

    trace_image = make_2d_trace_image(
        nx=nx,
        ny=ny,
        background=0,
        trace_center=trace_center,
        trace_order=trace_order,
        trace_coeffs=trace_coeffs,
        profile=source_profile,
        add_noise=False
    )

    spec_image = arc_image.data + trace_image.data + background

    if add_noise:
        from photutils.datasets import apply_poisson_noise
        spec_image = apply_poisson_noise(spec_image)

    ccd_im = CCDData(spec_image, unit=u.count, wcs=arc_image.wcs)

    return ccd_im
