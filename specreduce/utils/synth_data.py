# Licensed under a 3-clause BSD style license - see ../../licenses/LICENSE.rst
import warnings
from dataclasses import dataclass, field

import numpy as np
from astropy import units as u
from astropy.modeling import models, Model
from astropy.nddata import CCDData
from astropy.stats import gaussian_fwhm_to_sigma
from astropy.utils.decorators import deprecated
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

    Examples
    --------
    Build a traced continuum source with background and noise, then render it::

        from astropy.modeling import models
        from specreduce.utils.synth_data import SynthImage

        image = (
            SynthImage(nx=1024, ny=400, seed=42)
            .add_background(5)
            .add_source(profile=models.Moffat1D(amplitude=20, alpha=0.1))
            .add_poisson_noise()
            .add_read_noise(3)
            .to_ccddata()
        )

    See the :ref:`synth_data` guide for more examples, including tilted arc
    lines and modeling a non-linear dispersion relation.
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
        profile: Model = None,
        trace_center: float | None = None,
        trace_order: int = 3,
        trace_coeffs: dict | None = None,
    ) -> "SynthImage":
        """Add a continuum source with a Chebyshev-traced spatial profile."""
        if profile is None:
            profile = models.Moffat1D(amplitude=10, alpha=0.1)
        layer = SourceLayer(profile, trace_center, trace_order, trace_coeffs)
        return self._clone(_layers=self._layers + (layer,))

    def add_arcs(
        self,
        linelists=("HeI",),
        line_fwhm: float = 5.0,
        amplitude_scale: float = 1.0,
        wave_air: bool = False,
        tilt_func: Model = None,
    ) -> "SynthImage":
        """Add emission lines from one or more pypeit calibration line lists."""
        if tilt_func is None:
            tilt_func = models.Legendre1D(degree=0)
        if isinstance(linelists, str):
            linelists = (linelists,)
        layer = ArcLayer(
            tuple(linelists), line_fwhm, amplitude_scale, wave_air, tilt_func
        )
        return self._clone(_layers=self._layers + (layer,))

    def add_skylines(self, linelists="OH_GMOS", **kwargs) -> "SynthImage":
        """Add night-sky airglow emission lines (OH lists), wrapping ``add_arcs``."""
        return self.add_arcs(linelists, **kwargs)

    def add_poisson_noise(self) -> "SynthImage":
        """Apply Poisson noise to the rendered signal (requires photutils)."""
        return self._clone(_poisson=True)

    def add_read_noise(self, sigma: float) -> "SynthImage":
        """Add Gaussian read noise of standard deviation ``sigma`` (counts)."""
        return self._clone(_read_noise=sigma)

    add_rdnoise = add_read_noise

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
    """
    A continuum source whose spatial profile follows a Chebyshev trace.

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
    """
    Emission lines from one or more pypeit calibration line lists.

    Requires a resolvable WCS (supplied on ``SynthImage`` or built from its
    ``extent``). ``tilt_func`` applies a cross-dispersion tilt to simulate
    curved lines.
    """
    linelists: tuple
    line_fwhm: float = 5.0
    amplitude_scale: float = 1.0
    wave_air: bool = False
    tilt_func: Model = field(default_factory=lambda: models.Legendre1D(degree=0))

    def render(self, ctx: _RenderContext) -> np.ndarray:
        xx, yy = ctx.xx, ctx.yy
        if self.tilt_func is not None:
            if not isinstance(self.tilt_func, _ALLOWED_TILT):
                raise ValueError(
                    "The only tilt functions currently supported are 1D polynomials "
                    "from astropy.models."
                )
            if ctx.disp_axis == 0:
                xx = xx + self.tilt_func((yy - ctx.ny / 2) / ctx.ny)
            else:
                yy = yy + self.tilt_func((xx - ctx.nx / 2) / ctx.nx)

        z = np.zeros((ctx.ny, ctx.nx))
        linelist = load_pypeit_calibration_lines(list(self.linelists), wave_air=self.wave_air)
        if linelist is not None:
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore", message="No observer defined on WCS.*")
                line_disp_positions = ctx.wcs.spectral.world_to_pixel(linelist["wavelength"])
            line_sigma = gaussian_fwhm_to_sigma * self.line_fwhm
            for line_pos, ampl in zip(line_disp_positions, linelist["amplitude"]):
                line_mod = models.Gaussian1D(
                    amplitude=ampl * self.amplitude_scale,
                    mean=line_pos,
                    stddev=line_sigma,
                )
                if ctx.disp_axis == 0:
                    z += line_mod(xx)
                else:
                    z += line_mod(yy)
        return z


@deprecated("1.10", alternative="SynthImage")
def make_2d_trace_image(
    nx: int = 3000,
    ny: int = 1000,
    background: int | float = 5,
    trace_center: int | float | None = None,
    trace_order: int = 3,
    trace_coeffs: dict | None = None,
    profile: Model = None,
    add_noise: bool = True,
) -> CCDData:
    """
    Deprecated. Use :class:`SynthImage` instead.

    Equivalent to ``SynthImage(nx, ny).add_background(background)
    .add_source(...)`` followed by ``.add_poisson_noise()`` when ``add_noise``,
    then ``.to_ccddata()``.
    """
    if profile is None:
        profile = models.Moffat1D(amplitude=10, alpha=0.1)
    img = (
        SynthImage(nx=nx, ny=ny)
        .add_background(background)
        .add_source(
            profile=profile,
            trace_center=trace_center,
            trace_order=trace_order,
            trace_coeffs=trace_coeffs,
        )
    )
    if add_noise:
        img = img.add_poisson_noise()
    return img.to_ccddata()


@deprecated("1.10", alternative="SynthImage")
def make_2d_arc_image(
    nx: int = 3000,
    ny: int = 1000,
    wcs: WCS | None = None,
    extent=(3500, 7000),
    wave_unit: u.Unit = u.Angstrom,
    wave_air: bool = False,
    background: int | float = 5,
    line_fwhm: float = 5.0,
    linelists=("HeI",),
    amplitude_scale: float = 1.0,
    tilt_func: Model = None,
    add_noise: bool = True,
) -> CCDData:
    """Deprecated. Use :class:`SynthImage` with ``.add_arcs(...)`` instead."""
    if tilt_func is None:
        tilt_func = models.Legendre1D(degree=0)
    img = (
        SynthImage(nx=nx, ny=ny, wcs=wcs, extent=extent, wave_unit=wave_unit)
        .add_background(background)
        .add_arcs(
            linelists=linelists,
            line_fwhm=line_fwhm,
            amplitude_scale=amplitude_scale,
            wave_air=wave_air,
            tilt_func=tilt_func,
        )
    )
    if add_noise:
        img = img.add_poisson_noise()
    return img.to_ccddata()


@deprecated("1.10", alternative="SynthImage")
def make_2d_spec_image(
    nx: int = 3000,
    ny: int = 1000,
    wcs: WCS | None = None,
    extent=(6500, 9500),
    wave_unit: u.Unit = u.Angstrom,
    wave_air: bool = False,
    background: int | float = 5,
    line_fwhm: float = 5.0,
    linelists=("OH_GMOS",),
    amplitude_scale: float = 1.0,
    tilt_func: Model = None,
    trace_center: int | float | None = None,
    trace_order: int = 3,
    trace_coeffs: dict | None = None,
    source_profile: Model = None,
    add_noise: bool = True,
) -> CCDData:
    """Deprecated. Use :class:`SynthImage` with ``.add_arcs(...)`` and
    ``.add_source(...)`` instead."""
    if tilt_func is None:
        tilt_func = models.Legendre1D(degree=0)
    if source_profile is None:
        source_profile = models.Moffat1D(amplitude=10, alpha=0.1)
    img = (
        SynthImage(nx=nx, ny=ny, wcs=wcs, extent=extent, wave_unit=wave_unit)
        .add_background(background)
        .add_arcs(
            linelists=linelists,
            line_fwhm=line_fwhm,
            amplitude_scale=amplitude_scale,
            wave_air=wave_air,
            tilt_func=tilt_func,
        )
        .add_source(
            profile=source_profile,
            trace_center=trace_center,
            trace_order=trace_order,
            trace_coeffs=trace_coeffs,
        )
    )
    if add_noise:
        img = img.add_poisson_noise()
    return img.to_ccddata()
