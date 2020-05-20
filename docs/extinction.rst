.. _extinction:

Atmospheric Extinction
======================

Introduction
------------

The Earth's atmosphere is highly chromatic in its transmission of light. The wavelength-dependence
is dominated by scattering in the optical (320-700 nm) and molecular features in the near-infrared
and infrared.

Supported Optical Extinction Models
-----------------------------------

`specreduce` offers support for average optical extinction models for a set of observatories:

.. csv-table::
    :header:  "Model Name", "Observatory", "Lat", "Lon", "Elevation (m)", "Ref"

    "ctio", "Cerro Tololo Int'l Observatory", "-30.165", "70.815", "2215", "1"
    "kpno", "Kitt Peak Nat'l Observatory", "31.963", "111.6", "2120", "2"
    "lapalma", "Roque de Los Muchachos Observatory", "28.764", "17.895", "2396", "3"
    "mko", "Mauna Kea Int'l Observatory", "19.828", "155.48", "4160", "4"
    "mtham", "Lick Observatory, Mt. Hamilton Station", "37.341", "121.643", "1290", "5"
    "paranal", "European Southern Obs., Paranal Station", "-24.627", "70.405", "2635", "6"
    "apo", "Apache Point Observatory", "32.780", "105.82", "2788", "7"



1. The CTIO extinction curve was originally distributed with IRAF and comes from the work of
`Stone & Baldwin (1983 MNRAS 204, 347) <https://ui.adsabs.harvard.edu/abs/1983MNRAS.204..347S/abstract>`_
plus `Baldwin & Stone (1984 MNRAS 206, 241) <https://ui.adsabs.harvard.edu/abs/1984MNRAS.206..241B/abstract>`_.
The first of these papers lists the points from 3200-8370 Å while
the second extended the flux calibration from 6056 to 10870 Å but the
derived extinction curve was not given in the paper.  The IRAF table
follows SB83 out to 6436 Å, the redder points presumably come from BS84
with averages used in the overlap region. More recent CTIO extinction
curves are shown as Figures in Hamuy et al.
(`1992, PASP 104, 533 <https://ui.adsabs.harvard.edu/abs/1992PASP..104..533H/abstract>`_ ;
`1994, PASP 106, 566 <https://ui.adsabs.harvard.edu/abs/1994PASP..106..566H/abstract>`_).

2. The KPNO extinction table was originally distributed with IRAF. The ultimate provenance of this data is unclear,
but it has been used as-is in this form for over 30 years.

3. Extinction table for Roque de Los Muchachos Observatory, La Palma.
Described in https://www.ing.iac.es/Astronomy/observing/manuals/ps/tech_notes/tn031.pdf.

4. Median atmospheric extinction data for Mauna Kea Observatory measured by the Nearby SuperNova
Factory: https://www.aanda.org/articles/aa/pdf/2013/01/aa19834-12.pdf.

5. Extinction table for Lick Observatory on Mt. Hamilton constructed from
https://mthamilton.ucolick.org/techdocs/standards/lick_mean_extinct.html.

6. Updated extinction table for ESO-Paranal taken from
https://www.aanda.org/articles/aa/pdf/2011/03/aa15537-10.pdf.

7. Extinction table for Apache Point Observatory. Based on the extinction table used for SDSS and
available at https://www.apo.nmsu.edu/arc35m/Instruments/DIS/ (https://www.apo.nmsu.edu/arc35m/Instruments/DIS/images/apoextinct.dat).

In each case, the extinction is given in magnitudes per airmass and the wavelengths are in Angstroms. Here is an example that
uses the `~specreduce.calibration_data.ObservatoryExtinction` class to load each supported model and plots the extinction in
magnitudes as well as fractional transmission as a function of wavelength:

.. plot::
    :include-source:

    import matplotlib.pyplot as plt
    from specreduce.calibration_data import ObservatoryExtinction, SUPPORTED_EXTINCTION_MODELS

    fig, ax = plt.subplots(2, 1, sharex=True)
    for observatory in SUPPORTED_EXTINCTION_MODELS:
        ext = ObservatoryExtinction(observatory=observatory)
        ax[0].plot(ext.wavelength, ext.extinction(), label=observatory)
        ax[1].plot(ext.wavelength, ext.transmission())
    ax[0].legend(fancybox=True, shadow=True)
    ax[1].set_xlabel(f"Wavelength ({ext.wavelength.unit})")
    ax[0].set_ylabel("Extinction (mag)")
    ax[1].set_ylabel("Transmission")
    plt.tight_layout()
    fig.show()

A convenience class, `~specreduce.calibration_data.AtmosphericTransmission`, is provided for loading data files containing atmospheric transmission versus wavelength.
The common use case for this would be loading the output of telluric models. By default it loads a telluric model for an airmass of 1 and
1 mm of precipitable water. Some resources for generating more realistic model atmospheric transmission spectra include
https://mwvgroup.github.io/pwv_kpno/1.0.0/documentation/html/index.html and http://www.eso.org/sci/software/pipelines/skytools/molecfit.

.. plot::
    :include-source:

    import matplotlib.pyplot as plt
    from specreduce.calibration_data import AtmosphericTransmission

    fig, ax = plt.subplots()
    ext_default = AtmosphericTransmission()
    ext_custom = AtmosphericTransmission(data_path="atm_transmission_secz1.5_1.6mm.dat")
    ax.plot(ext_default.wavelength, ext_default.transmission(), label=r"sec $z$ = 1; 1 mm H$_{2}$O", linewidth=1)
    ax.plot(ext_custom.wavelength, ext_custom.transmission(), label=r"sec $z$ = 1.5; 1.6 mm H$_{2}$O", linewidth=1)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, 1.12), ncol=2, fancybox=True, shadow=True)
    ax.set_xlabel("Wavelength (microns)")
    ax.set_ylabel("Transmission")
    fig.show()

`~specreduce.calibration_data.CustomAtmosphericExtinction` is provided to facilitate constructing atmospheric extinction models
from user provided arrays of for wavelength and extinction versus wavelength. It will honor any provided wavelength units and assume
wavelengths are given in microns if they are not provided. Likewise, the provided extinction curve is assumed to be in magnitudes per airmass
if no units are provided. Otherwise, the supported units are `~astropy.units.function.logarithmic.Magnitude` or, if it's linear,
`dimensionless_unscaled <https://docs.astropy.org/en/stable/units/standard_units.html#doc-dimensionless-unit>`_.

.. plot::
    :include-source:

    import numpy as np

    import astropy.units as u

    import matplotlib.pyplot as plt
    from specreduce.calibration_data import CustomAtmosphericExtinction

    wave_min, wave_max = 0.3, 2.0
    wave = np.linspace(wave_min, wave_max, 50)

    # These are the exact same extinction curve expressed as magnitudes
    # and then linear transmission
    extinction_mag = np.sqrt(wave_max / wave) * u.mag
    extinction_linear = 10**(-0.4*np.sqrt(wave_max / wave)) * u.dimensionless_unscaled

    # Apply some constant offsets in magnitudes
    extinction_01mag = extinction_mag + 0.1 * u.mag
    extinction_02mag = extinction_mag + 0.2 * u.mag

    ext_lin = CustomAtmosphericExtinction(wavelength=wave, extinction_curve=extinction_linear)
    ext_01 = CustomAtmosphericExtinction(wavelength=wave, extinction_curve=extinction_01mag)
    ext_02 = CustomAtmosphericExtinction(wavelength=wave, extinction_curve=extinction_02mag)

    fig, ax = plt.subplots()
    ax.plot(ext_lin.wavelength, ext_lin.transmission(), label="Linear input")
    ax.plot(ext_01.wavelength, ext_01.transmission(), label="+0.1 mag")
    ax.plot(ext_02.wavelength, ext_02.transmission(), label="+0.2 mag")
    ax.legend(fancybox=True, shadow=True)
    ax.set_xlabel(f"Wavelength ({ext_lin.wavelength.unit})")
    ax.set_ylabel("Transmission")
    fig.show()

These classes are sub-classed from `~specreduce.calibration_data.BaseAtmosphericExtinction` and instances of them are callable to
apply an extinction correction an input spectrum. Here is an example:

.. plot::
    :include-source:

    import numpy as np

    import astropy.units as u

    import matplotlib.pyplot as plt

    from specutils import Spectrum1D
    from specreduce.calibration_data import ObservatoryExtinction

    wave = np.linspace(4000, 8000, 100) * u.angstrom
    flux = np.ones_like(wave.value) * u.mJy
    spec = Spectrum1D(flux=flux, spectral_axis=wave)
    atmos_corr = ObservatoryExtinction(observatory="mko")
    airmass = 1.5
    spec_corr = atmos_corr(spec, airmass=airmass)

    fig, ax = plt.subplots()
    ax.plot(spec.spectral_axis, spec.flux, label="Input")
    ax.plot(spec_corr.spectral_axis, spec_corr.flux, label="Corrected")
    ax.legend(fancybox=True, shadow=True)
    ax.set_xlabel(f"Wavelength ({spec.spectral_axis.unit})")
    ax.set_ylabel(f"Flux ({spec.flux.unit})")
    ax.set_title(f"Apply Mauna Kea Observatory extinction model at airmass={airmass}")
    fig.show()
