Glossary
========

This glossary defines terminology used in spectroscopic data reduction,
with emphasis on long-slit spectroscopy and the specreduce package.

.. glossary::
   :sorted:

   1D Spectrum
      Flux (or flux density) as a function of wavelength, frequency, or pixel
      position. This is the primary output of spectroscopic reduction, extracted
      from a :term:`2D spectrum`. A 1D spectrum may or may not be wavelength
      calibrated or flux calibrated. In specreduce and specutils, 1D spectra are
      represented by the `specutils.Spectrum` class.

   2D Spectrum
      A two-dimensional image from a CCD or similar detector containing dispersed
      light from one :term:`or more sources <MOS>`. One axis corresponds to wavelength (the
      :term:`dispersion axis`) and the other to spatial position along the slit
      (the :term:`cross-dispersion axis`). Also called a spectral image. A 2D spectrum is
      the input to the extraction process that produces a :term:`1D spectrum` , and
      is the starting point for reduction in both :term:`long-slit spectroscopy` and
      :term:`fiber-fed spectroscopy`.

   Absorption Line
      A feature in a spectrum where flux falls below the local continuum level
      at a characteristic wavelength, produced when atoms or molecules absorb
      photons at specific energies. Absorption lines in stellar spectra
      reveal the composition and physical conditions of stellar atmospheres,
      while interstellar and :term:`telluric <telluric correction>` absorption
      lines trace intervening material along the line of sight.

   Airmass
      The path length of light through Earth's atmosphere, expressed as a ratio
      relative to the path length at zenith. In the plane-parallel approximation,
      airmass = sec(z), where z is the zenith angle. At zenith, airmass = 1.0; at
      60\ :math:`^\circ` from zenith, airmass :math:`\approx` 2.0. Airmass affects
      the amount of :term:`atmospheric extinction` and is used in :term:`flux
      calibration`.

   Aperture
      A defined region on a 2D spectral image from which flux is extracted. In
      long-slit spectroscopy, the aperture is typically centered on the
      :term:`trace` and has a specified :term:`aperture width`. The aperture
      defines which pixels contribute to the extracted :term:`1D spectrum`.

   Aperture Width
      The spatial extent (in pixels) over which flux is summed during
      :term:`extraction`. For :term:`boxcar extraction`, all pixels within the
      aperture width receive equal weight. The optimal width balances including
      most of the source flux while minimizing background noise.

   Arc Lamp
      A calibration light source that produces emission lines at known
      wavelengths, used for :term:`wavelength calibration`. Common arc lamps
      include helium (He), argon (Ar), neon (Ne), and mercury-argon (HgAr)
      combinations. The resulting calibration exposure is called an
      :term:`arc spectrum`.

   Arc Spectrum
      A calibration exposure taken with an :term:`arc lamp`, showing discrete
      :term:`emission lines <emission line>` at known wavelengths. Arc spectra
      are used to establish the :term:`wavelength solution` that maps pixel
      positions to wavelengths.

   Atmospheric Extinction
      The wavelength-dependent absorption and scattering of light by Earth's
      atmosphere. Atmospheric extinction increases at shorter wavelengths and
      with higher :term:`airmass`. Correction for atmospheric extinction is
      part of :term:`flux calibration`.

   Background
      Signal in a spectral image that does not originate from the science
      target, including sky emission, scattered light, detector dark current,
      and other instrumental artifacts. Background subtraction removes this
      contaminating signal before or during :term:`extraction`. In specreduce,
      the `~specreduce.background.Background` class handles background estimation
      and subtraction.

   Background Aperture
      The spatial region(s) used to measure the :term:`background` level. In
      :term:`two-sided background` subtraction, apertures are placed on both
      sides of the :term:`trace`; in :term:`one-sided background` subtraction,
      a single aperture is used on one side.

   Barycentric Correction
      A velocity correction that accounts for Earth's orbital motion around
      the Solar System barycenter (center of mass). Applying barycentric
      correction converts observed wavelengths to the reference frame where
      the Solar System barycenter is at rest, enabling precise radial velocity
      measurements. See also :term:`heliocentric correction`.

   Bias Frame
      A zero-second exposure that captures the electronic offset (bias level)
      of the detector. Bias subtraction is typically performed during
      :term:`preprocessing` to remove this fixed offset from science exposures.

   Binning
      Division of an axis into discrete sections. In :term:`trace` fitting,
      binning along the :term:`dispersion axis` allows independent measurement
      of the trace position in each bin, which are then fit with a smooth
      function. In detector readout, binning combines adjacent pixels to
      reduce read noise at the cost of spatial or spectral resolution.

   Boxcar Extraction
      A simple :term:`extraction` method that sums all pixel values within the
      :term:`aperture` with equal weights. Boxcar extraction is straightforward
      but does not account for the :term:`spatial profile` shape or
      pixel-to-pixel variations in noise. In specreduce, this is implemented
      by the `~specreduce.extract.BoxcarExtract` class. Compare with
      :term:`optimal extraction`.

   Calibrated 2D Image
      A :term:`2D spectrum` that has been processed to have a
      :term:`wavelength calibration` applied (so pixel positions correspond to
      known wavelengths) but has not been resampled or rectified. The original
      pixel values are preserved.

   Catalog Lines
      Reference wavelengths of :term:`emission lines <emission line>` from
      standard line lists, used to match against detected lines in an
      :term:`arc spectrum` during :term:`wavelength calibration`. Common
      catalogs include HeI, ArI, NeI, and HgI lines.

   CCDData
      An Astropy data class (`astropy.nddata.CCDData`) for representing CCD
      images with associated uncertainty, mask, and metadata. CCDData is a
      subclass of :term:`NDData` and is commonly used for spectroscopic images
      in specreduce.

   Centroid
      The flux-weighted center position of a feature, such as an
      :term:`emission line` or the :term:`spatial profile` of a source. In
      :term:`trace` fitting, centroiding is one method for determining the
      trace position at each :term:`dispersion axis` bin.

   Classification
      The process of identifying the type of astronomical object a spectrum
      represents (e.g., star, galaxy, quasar). Classification may be performed
      automatically through model fitting or manually through
      :term:`visual inspection`.

   Coadding
      Combining multiple spectra of the same object into a single spectrum,
      typically to improve signal-to-noise ratio or to merge spectra covering
      different wavelength ranges. Unlike :term:`stacking`, coadding
      specifically refers to combining spectra of the same source.

   Cross-Dispersion Axis
      The axis perpendicular to the :term:`dispersion axis` in a
      :term:`2D spectrum`. In long-slit spectroscopy, the cross-dispersion
      axis corresponds to spatial position along the slit. Also called the spatial
      axis.

   Dark Frame
      An exposure taken with the detector shutter closed, capturing the
      thermal signal (dark current) accumulated during an exposure. Dark
      subtraction during :term:`preprocessing` removes this signal from
      science exposures.

   Data Cube
      A three-dimensional data structure with two spatial dimensions and one
      spectral dimension. Data cubes are the standard output format for
      :term:`IFU` observations, where each :term:`spatial pixel (spaxel) <spaxel>`
      contains a complete spectrum. Also called a spectral data cube or hyperspectral cube.

   Dispersion
      The separation of light by wavelength, or quantitatively, the change in
      wavelength per pixel (:math:`d\lambda/d\text{pixel}`) along the
      :term:`dispersion axis`. A spectrograph with higher dispersion spreads the
      spectrum over more pixels, providing higher :term:`spectral resolution`.
      Not to be confused with :term:`resolving power`.

   Dispersion Axis
      The axis along which wavelength varies in a :term:`2D spectrum`. Light at
      different wavelengths falls on different positions along this axis.

   Emission Line
      A feature in a spectrum where flux exceeds the local continuum level at
      a characteristic wavelength, produced when atoms or molecules emit
      photons at specific energies. Emission lines range from spectrally
      unresolved features in arc lamps to broad features in active galactic nuclei.
      Emission lines from :term:`arc lamps <arc lamp>` are used for :term:`wavelength
      calibration`; emission lines from astronomical sources provide
      information about their physical conditions, composition, and
      kinematics.

   Extinction Curve
      The wavelength-dependent function describing :term:`atmospheric extinction`.
      Standard extinction curves (e.g., for specific observatories) provide
      the extinction in magnitudes per unit :term:`airmass` as a function of
      wavelength.

   Extraction
      The process of converting a :term:`2D spectrum` into a :term:`1D spectrum`
      by summing or averaging flux along the :term:`cross-dispersion axis`
      within a defined :term:`aperture`. In specreduce, the two main extraction
      methods are :term:`boxcar extraction` and :term:`optimal extraction`,
      performed by the `~specreduce.extract.BoxcarExtract` and
      `~specreduce.extract.HorneExtract` classes, respectively.

   Fiber-Fed Spectroscopy
      A spectroscopic technique in which optical fibers relay light from the
      telescope focal plane to the spectrograph. Fibers decouple the
      spectrograph from the telescope, allowing it to be mounted in a stable
      environment, and enable flexible positioning of inputs across the focal
      plane. Fiber-fed designs are used in :term:`MOS` instruments and :term:`IFU`
      instruments, where each fiber produces a :term:`2D spectrum` on the detector.

   FITS
      Flexible Image Transport System. A common file format for
      astronomical data, capable of storing images, tables, and metadata
      (headers). Most spectroscopic data, including raw and reduced spectra,
      are stored in FITS format.

   Flat Field
      A calibration exposure of a uniformly illuminated source (e.g., dome
      flat, twilight flat, or lamp flat) used to measure pixel-to-pixel
      sensitivity variations. Flat-field correction during :term:`preprocessing`
      divides science exposures by a normalized flat to remove these variations.

   Flux
      In the strict physical sense, energy per unit time per unit area
      (:math:`\text{W/m}^2`). In spectroscopy, "flux" is often used informally to refer
      to :term:`flux density` or simply to the dependent variable (brightness)
      of a spectrum. The specutils `~specutils.Spectrum` class uses ``flux`` as the
      attribute name for the spectral data values.

   Flux Calibration
      The process of converting spectral data from instrumental units
      (counts, ADU, electrons) to physical flux density units
      (e.g., :math:`\text{erg/s/cm}^2\text{/\AA}` or :math:`\text{W/m}^2\text{/nm}`).
      Flux calibration requires observations of :term:`standard stars <standard star>`
      with known spectra and corrections for :term:`atmospheric extinction`.

   Flux Density
      Flux per unit wavelength or frequency interval, the standard physical
      quantity for spectroscopic measurements. Common units include
      :math:`\text{erg/s/cm}^2\text{/\AA}`, :math:`\text{W/m}^2\text{/nm}`,
      and Jy (jansky). See also :term:`flux`.

   FWHM
      Full Width at Half Maximum. A measure of the width of a peak or line
      profile, defined as the width at which the intensity drops (for
      :term:`emission lines <emission line>`) or rises (for :term:`absorption
      lines <absorption line>`) to half its extreme value relative to the local
      continuum. FWHM is commonly used to characterize :term:`spectral resolution`
      and the :term:`spatial profile` of sources.

   Gaussian Profile
      A bell-shaped curve (normal distribution) often used to model spectral
      lines and :term:`spatial profiles <spatial profile>`. In trace fitting,
      fitting a Gaussian to the cross-dispersion profile is one method for
      determining the trace position.

   GWCS
      `Generalized World Coordinate System <https://gwcs.readthedocs.io/en/latest/>`_.
      A Python package providing a framework for coordinate transformations, including
      non-linear mappings. In specreduce, GWCS is used to represent :term:`wavelength
      solutions <wavelength solution>` produced by wavelength calibration.

   Heliocentric Correction
      A velocity correction that accounts for Earth's motion relative to the
      Sun. Applying heliocentric correction converts observed wavelengths to
      the reference frame where the Sun is at rest. See also
      :term:`barycentric correction`, which provides higher precision.

   IFU
      `Integral Field Unit <https://en.wikipedia.org/wiki/Integral_field_spectrograph>`_.
      An instrument or observing mode that obtains
      spectra simultaneously over a contiguous two-dimensional field of view,
      producing a :term:`data cube`. IFUs use techniques such as fiber
      bundles, image slicers, or lenslet arrays to sample the field. Also
      called integral field spectroscopy (IFS).

   Interpolated Profile
      In :term:`optimal extraction`, an empirical :term:`spatial profile`
      constructed by sampling the actual profile at multiple wavelength bins
      and interpolating between them. This accounts for wavelength-dependent
      variations in the profile shape. In specreduce's `~specreduce.extract.HorneExtract`,
      this is enabled with the ``interpolated_profile`` option.

   Inverse Variance Weighting
      A statistical method for combining measurements that weights each value
      by the inverse of its variance (:math:`1/\sigma^2`), giving more weight to
      higher-precision measurements. Inverse variance weighting is used in
      :term:`optimal extraction` to maximize signal-to-noise ratio.

   Line Matching
      The process of associating detected :term:`emission lines <emission
      line>` in an :term:`arc spectrum` with known wavelengths from
      :term:`catalog lines`. Line matching is a key step in
      :term:`wavelength calibration`, establishing the correspondence between
      pixel positions and wavelengths.

   Long-Slit Spectroscopy
      A spectroscopic observing technique using an elongated slit (typically
      a few arcseconds wide and tens of arcseconds to arcminutes long) to
      obtain spectra of extended objects or to capture spatial information
      along one dimension. The resulting :term:`2D spectrum` contains spectral
      information along one axis and spatial information along the other.

   Mask
      A boolean array indicating which pixels in an image should be excluded
      from analysis (True = masked/bad, False = good). Masks flag cosmic rays,
      bad pixels, saturated pixels, and other defects. Specreduce provides
      various :term:`mask treatment` options for handling masked pixels.

   Mask Treatment
      The strategy for handling masked or non-finite pixels during processing.
      Options in specreduce include: ``apply`` (use mask as-is), ``ignore``
      (discard existing mask), ``propagate`` (extend mask along cross-dispersion),
      ``zero_fill`` (replace masked values with 0), ``nan_fill`` (replace with
      NaN), and others. See the Mask Treatment documentation for details.

   MOS
      Multi-Object Spectroscopy. A technique for obtaining spectra of multiple
      objects simultaneously using multiple slits, fiber positioners, or
      other multiplexing methods. Each object produces its own
      :term:`2D spectrum` on the detector.

   NDData
      N-dimensional data. An Astropy base class (`astropy.nddata.NDData`)
      for representing data with associated uncertainty, mask, unit, and WCS.
      Specreduce accepts NDData objects as input and uses them internally.

   One-Sided Background
      A :term:`background` estimation method using a single aperture on one
      side of the :term:`trace`. This is useful when the background on one
      side is contaminated (e.g., by a nearby source). In specreduce, use
      `specreduce.background.Background.one_sided()` to create a one-sided background.

   Optimal Extraction
      An :term:`extraction` method (also called Horne extraction, [Horne1986]_) that weights
      pixels according to the :term:`spatial profile` and their uncertainties,
      maximizing the signal-to-noise ratio of the extracted spectrum. Optimal
      extraction accounts for the varying contribution of different pixels and
      can reject cosmic rays. In specreduce, this is implemented by
      `~specreduce.extract.HorneExtract`.

   Order Sorting Filter
      An optical filter used to block unwanted orders of diffraction
      (e.g., preventing 2nd-order blue light from overlapping with
      1st-order red light). Failure to use or account for these filters
      results in systematic errors during :term:`flux calibration`.

   Peak Method
      The algorithm used to locate the position of a spectral :term:`trace`
      within each bin along the :term:`dispersion axis`. Options in specreduce's
      `~specreduce.tracing.FitTrace` include ``max`` (pixel with maximum value), ``gaussian``
      (fit Gaussian and use center), and ``centroid`` (flux-weighted center).

   Pipeline
      An automated sequence of processing steps that transforms raw data into
      reduced, calibrated data products. A spectroscopic reduction pipeline
      typically includes :term:`preprocessing`, :term:`extraction`,
      :term:`wavelength calibration`, and :term:`flux calibration`. Specreduce
      provides the building blocks for constructing such pipelines.

   Pixel Scale
      The physical or wavelength interval corresponding to one pixel. In the
      spectral direction, pixel scale refers to the :term:`dispersion`
      (wavelength per pixel). In the spatial direction, it refers to the
      angular size per pixel (arcseconds/pixel).

   Preprocessing
      Initial processing steps applied to raw CCD images before
      :term:`extraction`, including bias subtraction, dark subtraction,
      :term:`flat field` correction, and bad pixel masking. Also called
      instrument signature removal. In the Astropy ecosystem, preprocessing
      is typically handled by the `ccdproc <https://ccdproc.readthedocs.io/en/latest/>`__
      package.

   Rectified Spectrum
      A spectrum that has been resampled so that the :term:`dispersion axis`
      is aligned with image rows or columns, with constant wavelength spacing.
      Rectification involves :term:`resampling`. A rectified
      2D spectrum has one axis that is purely spectral.

   Redshift
      The fractional change in wavelength of spectral features due to the
      motion of the source (Doppler shift) or the expansion of the universe
      (cosmological redshift), defined as
      :math:`z = (\lambda_\text{observed} - \lambda_\text{rest}) / \lambda_\text{rest}`.

   Reduction
      The process of converting raw observational data into calibrated,
      science-ready spectra. For spectroscopy, reduction typically includes
      :term:`preprocessing`, :term:`extraction`, :term:`background`
      subtraction, :term:`wavelength calibration`, and optionally
      :term:`flux calibration`. The term comes from "reducing" the complexity
      of the data, though in practice the data volume often increases.

   Reference Pixel
      An anchor point for polynomial :term:`wavelength solutions <wavelength
      solution>`, typically near the center of the detector. The wavelength
      solution coefficients are defined relative to this reference position.

   Resampling
      The process of binning spectral data onto a new wavelength or
      pixel grid. Resampling is needed when combining spectra with different
      wavelength solutions or when creating :term:`rectified spectra
      <rectified spectrum>`. Resampling should be performed with care to
      preserve flux and properly propagate uncertainties.

   Resolving Power
      The ability of a spectrograph to distinguish closely spaced spectral
      features, defined as :math:`R = \lambda/\Delta\lambda`, where
      :math:`\Delta\lambda` is the minimum resolvable
      wavelength difference (related to :term:`FWHM`). Higher resolving power
      means finer spectral detail can be resolved. Not to be confused with
      :term:`dispersion`.

   Row-Stacked Spectra
      A collection of :term:`1D spectra <1D spectrum>` stored as a 2D array
      with one spectrum per row, sharing a common spectral axis. This format
      is used by specutils `~specutils.Spectrum` when representing multiple spectra.

   Sensitivity Function
      The wavelength-dependent response of the instrument, relating observed
      counts to true flux. The sensitivity function accounts for telescope
      throughput, instrument efficiency, and detector quantum efficiency.
      It is derived during :term:`flux calibration` using
      :term:`standard star` observations.

   Sigma Clipping
      An iterative outlier rejection method that excludes values more than a
      specified number of standard deviations (sigma) from the mean or median.
      In specreduce, sigma clipping can be applied during :term:`background`
      estimation to reject cosmic rays and other outliers.

   Sky
      Emission from Earth's atmosphere, including airglow lines, scattered
      moonlight, and light pollution. Sky emission must be subtracted from
      spectra to isolate the astronomical signal. See :term:`sky subtraction`.

   Sky Subtraction
      The removal of additive :term:`sky` emission from spectra. In long-slit
      spectroscopy, sky is typically measured from regions of the slit away
      from the target and subtracted from the object spectrum. Sky subtraction
      is closely related to :term:`background` subtraction.

   Slit
      The narrow aperture at the focal plane of a spectrograph that defines
      the region of the sky being observed. In :term:`long-slit spectroscopy`,
      the slit is elongated in one direction to capture spatial information.
      The slit width affects :term:`spectral resolution`; the slit length
      determines the spatial coverage.

   Slit Curvature
      The geometric distortion where a straight slit appears curved on
      the :term:`2D spectrum`. This occurs because the off-axis
      rays from the slit meet the grating at different angles. This
      must be corrected during :term:`wavelength calibration` or
      :term:`tilt correction`.

   Spatial Axis
      See :term:`cross-dispersion axis`.

   Spatial Profile
      The distribution of flux across the slit in the :term:`cross-dispersion
      axis` direction at a given wavelength. The spatial profile is typically
      determined by the seeing conditions and telescope optics, often
      approximated as a Gaussian. In :term:`optimal extraction`, accurate
      modeling of the spatial profile is essential for maximizing
      signal-to-noise.

   Spaxel
      A spatial pixel in a :term:`data cube`, corresponding to a single
      spatial element of an :term:`IFU` observation. Each spaxel contains a
      complete spectrum spanning the spectral axis of the data cube. The
      spatial sampling of a spaxel is determined by the IFU design and may not
      be square.

   Spectral Resolution
      The ability to distinguish spectral features at nearby wavelengths,
      often characterized by the :term:`FWHM` of unresolved lines. Can be
      expressed as the resolution element :math:`\Delta\lambda` or as
      :term:`resolving power` :math:`R = \lambda/\Delta\lambda`. Higher
      spectral resolution reveals finer detail in spectra.

   Spectrum
      The distribution of light intensity as a function of wavelength or
      frequency. In specreduce and specutils, `~specutils.Spectrum` specifically
      refers to the specutils class representing spectral data with flux,
      spectral axis, and optional uncertainty and mask.

   Stacking
      Combining multiple spectra or images, typically to increase
      signal-to-noise ratio. In some contexts, stacking refers specifically
      to combining data without alignment (e.g., NumPy array stacking), while
      :term:`coadding` implies proper combination of related spectra.

   Standard Star
      A star with a well-calibrated spectrophotometric spectrum, used as a
      reference for :term:`flux calibration`. Observations of standard stars
      under the same conditions as science targets allow determination of the
      :term:`sensitivity function`. Spectrophotometric standards include
      Vega, BD+17 4708, and various white dwarfs.

   Telluric Correction
      The removal of absorption features caused by molecules in Earth's
      atmosphere (telluric features), such as :math:`\text{O}_2`,
      :math:`\text{H}_2\text{O}`, and :math:`\text{CO}_2` bands.
      Telluric correction addresses the multiplicative effect of atmospheric
      absorption, distinct from the additive :term:`sky` emission. Methods
      include division by telluric standard star spectra or model fitting.

   Tilt Correction
      The process of correcting geometric distortions in a :term:`2D spectrum`
      where lines of constant wavelength are not aligned with detector columns
      (or rows). These distortions arise from :term:`slit curvature` and optical
      aberrations in the spectrograph, causing :term:`emission lines <emission
      line>` to appear tilted or curved across the :term:`cross-dispersion axis`.
      Tilt correction fits a 2D polynomial transformation between tilt-corrected and
      detector coordinate spaces using :term:`arc spectra <arc spectrum>`. In
      specreduce, this is handled by ``TiltCorrection`` (calibration) and
      ``TiltSolution`` (transformation and resampling).

   Trace
      The path of a spectrum across the :term:`2D spectrum` image, defining
      where the spectrum falls on the detector as a function of wavelength.
      Due to optical distortions, the trace may not be a straight line. In
      specreduce, trace classes include `~specreduce.tracing.FlatTrace` (constant
      position), `~specreduce.tracing.ArrayTrace` (arbitrary positions), and
      `~specreduce.tracing.FitTrace` (fit to image features).

   Trace Fitting
      The process of determining the :term:`trace` position as a function of
      wavelength by measuring the spectrum's location in bins along the
      :term:`dispersion axis` and fitting a smooth function (typically a
      polynomial) to these positions. In specreduce, `~specreduce.tracing.FitTrace`
      performs trace fitting with configurable :term:`peak method` and polynomial order.

   Two-Sided Background
      A :term:`background` estimation method using apertures on both sides of
      the :term:`trace`, symmetrically placed at equal distances (separation)
      from the center. The background values from both sides are averaged.
      In specreduce, use `specreduce.background.Background.two_sided()` to
      create a two-sided background.

   Uncertainty
      A measure of the error or noise associated with each data value.
      Uncertainties are typically expressed as standard deviations or
      variances and should be propagated through all reduction steps.
      Specreduce supports uncertainty propagation in background subtraction
      and extraction.

   Variance
      The square of the standard deviation, representing the expected squared
      deviation from the mean. Variances are used in :term:`inverse variance
      weighting` and are propagated through calculations according to error
      propagation rules.

   Visual Inspection
      Human examination of spectra or data products to verify quality, identify
      features, or make scientific judgments. Visual inspection remains important
      for assessing reduction quality and identifying issues that automated
      methods may miss.

   Wavelength Calibration
      The process of determining the relationship between pixel position and
      wavelength, establishing a :term:`wavelength solution`. Wavelength
      calibration typically uses :term:`arc spectra <arc spectrum>` with known
      :term:`emission lines <emission line>`. In specreduce, this is handled
      by the `~specreduce.wavecal1d.WavelengthCalibration1D` class, which produces a
      :term:`GWCS`-based WCS.

   Wavelength Solution
      A mathematical function (typically a polynomial) that maps pixel
      positions to wavelengths. The wavelength solution is determined through
      :term:`wavelength calibration` and is applied to the data as a
      :term:`WCS`.

   WCS
      World Coordinate System. A standard for describing the mapping between
      pixel coordinates and physical coordinates such as wavelength, sky
      position, or other quantities. For spectra, the WCS maps pixel numbers
      to wavelengths. The FITS WCS standard and :term:`GWCS` are common
      implementations. See the `Astropy WCS documentation
      <https://docs.astropy.org/en/stable/wcs/>`_ for details.

   Window
      A restricted region used for analysis, such as a spatial window for
      limiting :term:`trace` fitting to a subset of the
      :term:`cross-dispersion axis` (useful when multiple objects are present),
      or a spectral window limiting analysis to a wavelength range.

   Workflow
      A complete sequence of processing and analysis steps, often broader than
      a :term:`pipeline`, potentially including data organization, job
      orchestration, quality control, and archiving in addition to the core
      reduction steps.

References
----------

.. [Horne1986] `Horne, K. (1986). "An optimal extraction algorithm for
   CCD spectroscopy." Publications of the Astronomical Society of the
   Pacific, 98, 609. <https://ui.adsabs.harvard.edu/abs/1986PASP...98..609H/abstract>`__