"""
Utilities for defining, loading, and handling spectroscopic calibration data
"""

import os
import pkg_resources
import warnings
from dataclasses import dataclass

import numpy as np

import astropy.units as u
from astropy.table import Table
from astropy.utils.data import download_file
from astropy.utils.exceptions import AstropyUserWarning

import synphot
from specutils import Spectrum1D
from specutils.manipulation import FluxConservingResampler, SplineInterpolatedResampler

__all__ = [
    'get_reference_file_path',
    'load_MAST_calspec',
    'load_onedstds',
    'BaseAtmosphericExtinction',
    'AtmosphericTransmission',
    'ObservatoryExtinction',
    'CustomAtmosphericExtinction'
]

"""
Make specreduce_data optional. If it's available, great and we can access its data via
pkg_resources. If not, we'll fall back to downloading and optionally caching it using
`~astropy.utils.data`.
"""
LOCAL_DATA = True
try:
    import specreduce_data  # noqa
except ModuleNotFoundError:
    warnings.warn(
        "Can't import specreduce_data package. Falling back to downloading data...",
        AstropyUserWarning
    )
    LOCAL_DATA = False

SUPPORTED_EXTINCTION_MODELS = [
    "kpno",
    "ctio",
    "apo",
    "lapalma",
    "mko",
    "mtham",
    "paranal"
]
"""
List of available extinction curves for individual observatories
"""

SPECPHOT_DATASETS = [
    "bstdscal",
    "ctiocal",
    "ctionewcal",
    "eso",
    "gemini",
    "iidscal",
    "irscal",
    "oke1990",
    "redcal",
    "snfactory",
    "spec16cal",
    "spec50cal",
    "spechayescal"
]
"""
List of spectrophotometric standard star datasets
"""


def get_reference_file_path(path=None, cache=False, show_progress=False):
    """
    Basic function to take a path to a file and load it via pkg_resources if the `specreduce_data`
    package is available and load it via ghithub raw user content if not.

    Parameters
    ----------
    path : str or None (default: None)
        Filename of reference file relative to the reference_data directory within
        `specreduce_data` package.

    cache : bool (default: False)
        Set whether file is cached if file is downloaded.

    show_progress : bool (default: False)
        Set whether download progress bar is shown if file is downloaded.

    Returns
    -------
    file_path : str or None
        Path to reference data file or None if the path cannot be constructed or if the file
        itself is not valid.

    Examples
    --------
    >>> from specreduce.calibration_data import get_reference_file_path
    >>> kpno_extinction_file = get_reference_file_path("extinction/kpnoextinct.dat")
    """
    if path is None:
        return None

    if LOCAL_DATA:
        file_path = pkg_resources.resource_filename(
            "specreduce_data",
            os.path.join("reference_data", path)
        )
    else:
        repo_url = "https://raw.githubusercontent.com/astropy/specreduce-data"
        remote_url = f"{repo_url}/master/specreduce_data/reference_data/{path}"
        try:
            file_path = download_file(
                remote_url,
                cache=cache,
                show_progress=show_progress,
                pkgname='specreduce'
            )
        except Exception as e:
            msg = f"Downloading of {path} failed: {e}"
            warnings.warn(msg, AstropyUserWarning)
            return None

    # final sanity check to make sure file_path is actually a file.
    if os.path.isfile(file_path):
        return file_path
    else:
        warnings.warn(f"Able to construct {file_path}, but it is not a file.")
        return None


def load_MAST_calspec(filename, remote=True, cache=True, show_progress=False):
    """
    Load a standard star spectrum from the `calspec` database at MAST. These spectra are provided in
    FITS format and are described in detail at:

    https://www.stsci.edu/hst/instrumentation/reference-data-for-calibration-and-tools/astronomical-catalogs/calspec  # noqa

    If `remote` is True, the spectrum will be downloaded from MAST. Set `remote` to False to load
    a local file.

    Parameters
    ----------
    filename : str
        FITS filename of the standard star spectrum, e.g. g191b2b_005.fits.

    remote : bool (default = True)
        If True, download the spectrum from MAST. If False, check if `filename` exists and load it.
    cache : bool (default = True)
        Toggle whether downloaded data is cached or not.
    show_progress : bool (default = True)
        Toggle whether download progress bar is shown.

    Returns
    -------
    spectrum : None or `~specutils.Spectrum1D`
        If the spectrum can be loaded, return it as a `~specutils.Spectrum1D`.
        Otherwise return None. The spectral_axis units are Å (`~astropy.units.angstrom`) and
        the flux units are milli-Janskys (`~astropy.units.mJy`).
    """
    if remote:
        url = f"https://archive.stsci.edu/hlsps/reference-atlases/cdbs/calspec/{filename}"
        try:
            file_path = download_file(
                url,
                cache=cache,
                show_progress=show_progress,
                pkgname='specreduce'
            )
        except Exception as e:
            msg = f"Downloading of {filename} failed: {e}"
            warnings.warn(msg, AstropyUserWarning)
            file_path = None
    else:
        if os.path.isfile(filename):
            file_path = filename
        else:
            msg = f"Provided filename, {filename}, does not exist or is not a valid file."
            warnings.warn(msg, AstropyUserWarning)
            file_path = None

    if file_path is None:
        return None
    else:
        hdr, wave, flux = synphot.specio.read_fits_spec(file_path)

        # the calspec data stores flux in synphot's FLAM units. convert to flux units
        # supported directly by astropy.units. mJy is chosen since it's the JWST
        # standard and can easily be converted to/from AB magnitudes.
        flux_mjy = synphot.units.convert_flux(wave, flux, u.mJy)
        spectrum = Spectrum1D(spectral_axis=wave, flux=flux_mjy)
        return spectrum


def load_onedstds(dataset="snfactory", specfile="EG131.dat", cache=True, show_progress=False):
    """
    This is a convenience function for loading a standard star spectrum from the 'onedstds'
    dataset in the `specreduce_data` package. If that package is installed, `~pkg_resources`
    will be used to locate the data files locally. Otherwise they will be downloaded from the
    repository on github.

    Parameters
    ----------
    dataset : str  (default = "snfactory")
        Standard star spectrum database. Valid options are described in :ref:`specphot_standards`.

    specfile : str (default = "EG131.dat")
        Filename of the standard star spectrum.

    cache : bool (default = True)
        Enable caching of downloaded data.

    show_progress : bool (default = False)
        Show download progress bar if data is downloaded.

    Returns
    -------
    spectrum : None or `~specutils.Spectrum1D`
        If the spectrum can be loaded, return it as a `~specutils.Spectrum1D`.
        Otherwise return None. The spectral_axis units are Å (`~astropy.units.angstrom`) and
        the flux units are milli-Janskys (`~astropy.units.mJy`).
    """
    if dataset not in SPECPHOT_DATASETS:
        msg = (f"Specfied dataset, {dataset}, not in list of supported datasets of "
               f"spectrophotometric standard stars: f{SPECPHOT_DATASETS}")
        warnings.warn(msg, AstropyUserWarning)
        return None

    spec_path = get_reference_file_path(
        path=os.path.join("onedstds", dataset, specfile),
        cache=cache,
        show_progress=show_progress
    )
    if spec_path is None:
        msg = f"Can't load {specfile} from {dataset}."
        warnings.warn(msg, AstropyUserWarning)
        return None

    t = Table.read(spec_path, format="ascii", names=['wavelength', 'ABmag', 'binsize'])

    # the specreduce_data standard star spectra all provide wavelengths in angstroms
    spectral_axis = t['wavelength'].data * u.angstrom

    # the specreduce_data standard star spectra all provide fluxes in AB mag
    flux = t['ABmag'].data * u.ABmag
    flux = flux.to(u.mJy)  # convert to linear flux units
    spectrum = Spectrum1D(spectral_axis=spectral_axis, flux=flux)
    return spectrum


@dataclass
class BaseAtmosphericExtinction:
    """
    Base atmospheric extinction class to set up common methods/properties

    Attributes
    ----------
    wavelength : `~astropy.units.Quantity`
        Wavelengths of the extinction curve

    extinction_curve : `~astropy.units.Quantity`
        Base extinction curve normalized to an airmass of 1.0 in units of magnitudes/airmass

    Parameters
    ----------
    object_spectrum : `specutils.Spectrum1D` or spectrum-like
        Input spectrum to correct for atmospheric transmission.
    airmass : float (default = 1.0)
        Airmass that `object_spectrum` was observed at.
    resampler : instance of a `~specutils.manipulation.ResamplerBase` subclass or None
        (Optional) Pass in a custom spectral resampler. If None, then the spectral axis samplings
        of the extinction curve and input spectrum are used to decide which resampler to use.
        `~specutils.manipulation.SplineInterpolatedResampler` is used if the extinction curve is
        more coarsely sampled and `~specutils.manipulation.FluxConservingResampler` is used if
        it's not.

    Returns
    -------
    corrected_spectrum : `specutils.Spectrum1D` or spectrum-like
        Spectrum object with flux corrected for atmospheric transmission.
    """
    def __post_init__(self):
        """
        Set up a default null extinction curve with the right units
        """
        self.wavelength = np.linspace(1000, 50000, 10) * u.angstrom
        self.extinction_curve = u.Magnitude(
            np.zeros_like(self.wavelength.data),
            u.MagUnit(u.dimensionless_unscaled)
        )

    def __call__(self, object_spectrum, airmass=1.0, resampler=None):
        """
        Make instance callable to apply extinction model to an input spectrum
        """
        if resampler is None:
            """
            If the wavelength sampling of the extinction curve is coarser than the input
            spectrum's, we want to use `~specutils.manipulation.SplineInterpolatedResampler`
            to estimate intermediate values. However, if the extinction curve is more finely
            sampled than the input (e.g. high resolution spectra of telluric features), we want
            to use `~specutils.manipulation.FluxConservingResampler` so that the absorption is
            properly integrated within each input wavelength bin.
            """
            object_spacing = np.diff(object_spectrum.spectral_axis).max()
            extinction_spacing = np.diff(self.wavelength).min()
            if extinction_spacing > object_spacing:
                resampler = SplineInterpolatedResampler()
            else:
                resampler = FluxConservingResampler(extrapolation_treatment='nan_fill')

        transmission = self.transmission(airmass=airmass)
        # specutils resamplers want the same wavelength units for both spectra
        wave = self.wavelength.to(object_spectrum.spectral_axis_unit)
        trans_spec = Spectrum1D(flux=transmission, spectral_axis=wave)
        resampled_transmission = resampler(trans_spec, object_spectrum.spectral_axis)
        corrected_spectrum = object_spectrum / resampled_transmission
        return corrected_spectrum

    def extinction(self, airmass=1.0):
        """
        Return extinction in magnitudes at a given airmass

        Parameters
        ----------
        airmass : float (default = 1.0)
            Airmass to scale extinction curve to.

        Returns
        -------
        ext : `~astropy.units.Quantity`
            Airmass-scaled extinction curve in magnitudes.
        """
        if airmass < 1.0:
            msg = f"Airmass, {airmass}, must be >= 1.0."
            raise ValueError(msg)
        # multiplying by a scalar converts Magnitude to Quantity so work around that
        ext = self._mult_airmass(airmass=airmass)
        return ext

    def transmission(self, airmass=1.0):
        """
        Return dimensionless transmission at a given airmass

        Parameters
        ----------
        airmass : float (default = 1.0)
            Airmass to scale extinction curve to.

        Returns
        -------
        ext : `~astropy.units.Quantity`
            Airmass-scaled transmission curve.
        """
        if airmass < 1.0:
            msg = f"Airmass, {airmass}, must be >= 1.0."
            raise ValueError(msg)
        ext = self._mult_airmass(airmass=airmass)
        trans = ext.to(u.dimensionless_unscaled)
        return trans

    def _mult_airmass(self, airmass=1.0):
        """
        Workaround astropy.units "feature" where multiplying by a scalar strips
        magnitudes of their physical units.

        Parameters
        ----------
        airmass : float (default = 1.0)
            Airmass to scale extinction curve to.

        Returns
        -------
        ext : `~astropy.units.Quantity`
            Airmass-scaled extinction curve in magnitudes.
        """
        ext = u.Magnitude(
            self.extinction_curve.value * airmass,
            u.MagUnit(u.dimensionless_unscaled)
        )
        return ext


@dataclass
class ObservatoryExtinction(BaseAtmosphericExtinction):
    """
    Load atmospheric extinction models for specific observatories that are
    provided by the `specreduce_data` package.

    Parameters
    ----------
    observatory : str
        Name of atmospheric extinction model provided by `specreduce_data`. Valid
        options are:

        kpno - Kitt Peak National Observatory (default)
        ctio - Cerro Tololo International Observatory
        apo - Apache Point Observatory
        lapalma - Roque de los Muchachos Observatory, La Palma, Canary Islands
        mko - Mauna Kea Observatories
        mtham - Lick Observatory, Mt. Hamilton station
        paranal - European Southern Observatory, Cerro Paranal station

    cache : bool (default = True)
        Toggle caching of downloaded data.

    show_progress : bool (default = False)
        Toggle showing progress bar while downloading data.
    """
    observatory: str = "kpno"
    cache: bool = True
    show_progress: bool = False

    def __post_init__(self):
        if self.observatory not in SUPPORTED_EXTINCTION_MODELS:
            msg = (
                f"Requested observatory extinction model, {self.observatory}, not in list "
                f"of available models: {SUPPORTED_EXTINCTION_MODELS}"
            )
            raise ValueError(msg)
        model_file = os.path.join("extinction", f"{self.observatory}extinct.dat")
        model_path = get_reference_file_path(
            path=model_file,
            cache=self.cache,
            show_progress=self.show_progress
        )
        t = Table.read(model_path, format="ascii", names=['wavelength', 'extinction'])

        # the specreduce_data models all provide wavelengths in angstroms
        self.wavelength = t['wavelength'].data * u.angstrom

        # the specreduce_data models all provide extinction in magnitudes at an airmass of 1
        self.extinction_curve = u.Magnitude(
            t['extinction'].data,
            u.MagUnit(u.dimensionless_unscaled)
        )


@dataclass
class AtmosphericTransmission(BaseAtmosphericExtinction):
    """
    Load atmospheric transmission vs wavelength from a data file, e.g. output of a
    telluric model.

    Parameters
    ----------
    data_path : str
        Path to file containing atmospheric transmission data. Data is assumed to have
        two columns, wavelength and transmission (unscaled dimensionless). If
        this isn't provided, a model is built from a pre-calculated table of values
        from 0.9 to 5.6 microns. The values were generated by the ATRAN model,
        https://atran.arc.nasa.gov/cgi-bin/atran/atran.cgi (Lord, S. D., 1992, NASA
        Technical Memorandum 103957). The extinction is given as a linear transmission
        fraction at an airmass of 1 and 1 mm of precipitable water.

    wave_unit : `~astropy.units.Unit` (default = u.um)
        Units for spectral axis.

    cache : bool (default = True)
        Toggle caching of downloaded data.

    show_progress : bool (default = False)
        Toggle showing progress bar while downloading data.
    """
    data_path: str = os.path.join("extinction", "atm_trans_am1.0.dat")
    wave_unit: u.Quantity = u.um
    cache: bool = True
    show_progress: bool = False

    def __post_init__(self):
        if not os.path.isfile(self.data_path):
            orig_path = self.data_path
            self.data_path = get_reference_file_path(
                path=self.data_path,
                cache=self.cache,
                show_progress=self.show_progress
            )

        if not os.path.isfile(self.data_path):
            msg = (
                f"Can't load atmospheric transmission data file, {orig_path}, locally "
                "or from specreduce_data."
            )
            raise ValueError(msg)

        t = Table.read(self.data_path, format="ascii", names=['wavelength', 'extinction'])

        # spectral axis is given in microns
        self.wavelength = t['wavelength'].data * self.wave_unit

        # extinction is given in a dimensionless transmission fraction
        transmission_curve = t['extinction'].data * u.dimensionless_unscaled
        self.extinction_curve = transmission_curve.to(u.mag(u.dimensionless_unscaled))


@dataclass
class CustomAtmosphericExtinction(BaseAtmosphericExtinction):
    """
    Custom atmospheric extinction model using input arrays for wavelength
    and extinction_curve.

    Parameters
    ----------
    wavelength : list-like or `~astropy.units.Quantity`
        Input array of wavelengths for extinction curve (default: u.um if not provided)

    extinction_curve : list-like or `~astropy.units.Quantity`
        Input extinction curve (default: magnitudes if not provided)
    """
    wavelength: u.Quantity
    extinction_curve: u.Quantity

    def __post_init__(self):
        if not hasattr(self.wavelength, '__len__'):
            msg = "Input wavelengths must be an array."
            raise ValueError(msg)

        if not hasattr(self.extinction_curve, '__len__'):
            msg = "Input extinction curve must be an array."
            raise ValueError(msg)

        if len(self.wavelength) < 2 or len(self.extinction_curve) < 2:
            msg = "Input extinction curve data must have a wavelength range, i.e. length >= 2."
            raise ValueError(msg)

        if not hasattr(self.wavelength, "unit"):
            self.wavelength *= u.um

        if not hasattr(self.extinction_curve, "unit"):
            self.extinction_curve = u.Magnitude(
                self.extinction_curve,
                u.MagUnit(u.dimensionless_unscaled)
            )

        if self.extinction_curve.unit == u.dimensionless_unscaled:
            # force self.extinction_curve into magnitudes
            self.extinction_curve = self.extinction_curve.to(u.mag(u.dimensionless_unscaled))
        elif self.extinction_curve.unit == u.mag:
            # if raw magnitudes, convert to physical ones
            if not hasattr(self.extinction_curve, "physical"):
                self.extinction_curve = u.Magnitude(
                    self.extinction_curve.value,
                    u.MagUnit(u.dimensionless_unscaled)
                )
        else:
            msg = f"Invalid units, {self.extinction_curve.unit}, for input extinction curve."
            raise ValueError(msg)
