# Target strength for a calibration sphere

## acousticTS implementation

Benchmarked Validated

[Overview](https://brandynlucca.github.io/acousticTS/articles/calibration/index.md)
[Theory](https://brandynlucca.github.io/acousticTS/articles/calibration/calibration-theory.md)

SOEMS evaluates solid elastic reference spheres used in standard-target
calibration ([Dragonette et al. 1981](#ref-Dragonette_1981); [Foote
1990](#ref-Foote_1990); [MacLennan 1981](#ref-Maclennan_1981)).

Create a `CAL` object, evaluate it with
[`target_strength()`](https://brandynlucca.github.io/acousticTS/reference/target_strength.md),
and inspect the stored spectrum. The object retains the sphere
dimensions and elastic material properties alongside the result.

### Calibration sphere object generation

A calibration sphere is created with
[`cal_generate()`](https://brandynlucca.github.io/acousticTS/reference/cal_generate.md).
`diameter` is supplied in metres. `material` selects a stored preset, or
the density and longitudinal and transverse sound speeds can be supplied
directly.

See
[`cal_generate()`](https://brandynlucca.github.io/acousticTS/reference/cal_generate.md)
for the current preset names and values. Keeping that list on the
function reference page avoids a second copy that can drift from the
package definitions.

When using the defaults:

``` r

library(acousticTS)

cal_sphere <- cal_generate()
# equivalent to: cal_generate(material = "WC", diameter = 38.1e-3)
```

### Calculating a target-strength spectrum

[`target_strength()`](https://brandynlucca.github.io/acousticTS/reference/target_strength.md)
returns an updated object. The registered model name is `"calibration"`.
`"SOEMS"` is its public alias.

``` r

frequency <- seq(1e3, 600e3, 1e3)

cal_sphere <- target_strength(
  object = cal_sphere,
  frequency = frequency,
  model = "calibration"
)
```

### Inspecting model results

Plot the stored spectrum for a quick check or use
[`extract()`](https://brandynlucca.github.io/acousticTS/reference/extract.md)
for downstream analysis.

#### Plotting results

The [`plot()`](https://rdrr.io/r/graphics/plot.default.html) method can
be used to display either the sphere geometry or the modeled output. For
calibration work, the most common use is `type = "model"`, which plots
the stored target-strength spectrum. The optional `x_units` argument can
also be used to display the horizontal axis in terms of frequency or in
terms of radius-scaled wavenumber.

![Pre-rendered calibration-sphere spectra shown against frequency and
three radius-scaled wavenumber axes for the default tungsten-carbide
sphere.](calibration-spectrum-frequency.png)![Pre-rendered
calibration-sphere spectra shown against frequency and three
radius-scaled wavenumber axes for the default tungsten-carbide
sphere.](calibration-spectrum-k-sw.png)![Pre-rendered calibration-sphere
spectra shown against frequency and three radius-scaled wavenumber axes
for the default tungsten-carbide
sphere.](calibration-spectrum-k-l.png)![Pre-rendered calibration-sphere
spectra shown against frequency and three radius-scaled wavenumber axes
for the default tungsten-carbide sphere.](calibration-spectrum-k-t.png)

Those alternatives are useful for different reasons. Frequency is the
natural axis for practical calibration work, while the radius-scaled
wavenumber views are useful when comparing spheres of different
diameters or different materials on a common nondimensional scale.

#### Accessing results

The model results can also be accessed directly with
[`extract()`](https://brandynlucca.github.io/acousticTS/reference/extract.md).
For the calibration workflow, `feature = "model"` returns a data frame
containing the stored spectral outputs.

``` r

model_results <- extract(cal_sphere, "model")$calibration
head(model_results)
```

    ##   frequency        ka         f_bs     sigma_bs        TS
    ## 1      1000 0.0810226 9.735761e-05 9.478503e-09 -80.23260
    ## 2      2000 0.1620452 3.853270e-04 1.484769e-07 -68.28341
    ## 3      3000 0.2430678 8.518051e-04 7.255720e-07 -61.39319
    ## 4      4000 0.3240904 1.477241e-03 2.182242e-06 -56.61097
    ## 5      5000 0.4051130 2.235373e-03 4.996892e-06 -53.01300
    ## 6      6000 0.4861356 3.093986e-03 9.572752e-06 -50.18963

The extracted data frame includes the working frequency grid, the
ambient acoustic-size variable `ka`, the reported backscattering length
`f_bs`, the backscattering cross-section `sigma_bs`, and target strength
`TS`. In the current calibration workflow, these quantities are related
by:

\sigma\_{\mathrm{bs}} = \|f\_{\mathrm{bs}}\|^2, \qquad \mathit{TS} = 10
\log\_{10}\left(\sigma\_{\mathrm{bs}}\right).

Calibration references do not always use the same amplitude
normalization. See the Theory page for the normalization used by SOEMS.

### Comparison workflows

Diameter and material comparisons isolate geometric scaling from elastic
material effects.

#### Diameter comparisons

![Pre-rendered calibration comparison showing how the stored
calibration-sphere spectrum shifts with sphere
diameter.](calibration-diameter-comparison.png)

Diameter changes alter the acoustic-size scaling directly, so the
resonance structure shifts across the frequency axis even when the
sphere material is unchanged. That is one reason calibration practice is
usually tied to standard diameters rather than to an abstract material
class alone.

#### Material comparisons

![Pre-rendered calibration comparison showing how the stored
calibration-sphere spectrum varies with sphere material at fixed
diameter.](calibration-material-comparison.png)

Material comparisons are especially informative because they isolate the
role of elastic wave speeds and density from the purely geometric role
of sphere size. A tungsten carbide sphere and an aluminum sphere of the
same diameter do not simply differ by a vertical offset. Their resonance
structure can also shift because the interior compressional and shear
wave speeds have changed.

#### External implementation comparison

The external check compares the MacLennan elastic-sphere formulation
with SphereTS, echoSMs, and the NWFSC calibration applet ([MacLennan
1981](#ref-Maclennan_1981); [Macaulay 2025](#ref-SphereTS_software);
[Macaulay and contributors 2024](#ref-echoSMs_software)). With
`adaptive = TRUE` (the default), the solver starts at
\operatorname{round}(ka)+10 partial waves and extends the sum until its
tail falls below 10^{-10}. `adaptive = FALSE` uses only the initial
fixed cutoff. The 38.1 mm tungsten-carbide comparison is limited to
1–360 kHz so the applet remains within its stated ka \lesssim 30 range.

| Comparison | N frequency | Max abs. \Delta TS (dB) | Mean abs. \Delta TS (dB) |
|:---|---:|---:|---:|
| acousticTS vs echoSMs | 360 | 0 | 0 |
| acousticTS vs sphereTS | 360 | 0 | 0 |
| acousticTS vs NOAA applet | 360 | 0 | 0 |
| echoSMs vs sphereTS | 360 | 0 | 0 |
| echoSMs vs NOAA applet | 360 | 0 | 0 |
| sphereTS vs NOAA applet | 360 | 0 | 0 |

![Pre-rendered calibration comparison against echoSMs, sphereTS, and the
NOAA calibration applet for the 38.1 mm tungsten-carbide
sphere.](calibration-external-comparison.png)

For the 38.1 mm tungsten-carbide sphere, `adaptive = TRUE` agrees with
the other implementations to about 10^{-10} dB. The fixed cutoff remains
close, with a maximum difference of about 7.2 \times 10^{-5} dB.

To show that this is not unique to the 38.1 mm tungsten-carbide sphere,
the same comparison was repeated for one smaller tungsten-carbide sphere
and one copper sphere from the calibration-target definitions shipped
with echoSMs ([Macaulay and contributors 2024](#ref-echoSMs_software)),
again including the `SphereTS` implementation ([Macaulay
2025](#ref-SphereTS_software)).

| Target | Diameter (mm) | N frequency | Max frequency (kHz) | Max abs. \Delta adapt = TRUE vs echoSMs (dB) | Max abs. \Delta adapt = FALSE vs echoSMs (dB) | Max abs. \Delta adapt = TRUE vs sphereTS (dB) | Max abs. \Delta adapt = FALSE vs sphereTS (dB) | Max abs. \Delta adapt = TRUE vs NOAA applet (dB) | Max abs. \Delta adapt = FALSE vs NOAA applet (dB) |
|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| WC20 calibration sphere | 20.0 | 360 | 360 | 0 | 1.0e-06 | 0 | 1.0e-06 | 0 | 1.0e-06 |
| WC38.1 calibration sphere | 38.1 | 360 | 360 | 0 | 7.2e-05 | 0 | 7.2e-05 | 0 | 7.2e-05 |
| Cu32.1 calibration sphere | 32.1 | 360 | 360 | 0 | 4.5e-05 | 0 | 4.5e-05 | 0 | 4.5e-05 |

Across the additional targets, the adaptive solver keeps the maximum
absolute differences near 10^{-10} dB. The fixed cutoff remains within
about 10^{-5} to 10^{-4} dB of the other implementations.

### Closing note

SOEMS provides a compact workflow for constructing a reference sphere,
computing its spectrum, and comparing diameter, material, or solver
settings.

## References

Dragonette, Louis R., S. K. Numrich, and Laurence J. Frank. 1981.
“Calibration Technique for Acoustic Scattering Measurements.” *The
Journal of the Acoustical Society of America* 69 (4): 1186–89.
<https://doi.org/10.1121/1.385699>.

Foote, K. G. 1990. “Spheres for Calibrating an Eleven-Frequency Acoustic
Measurement System.” *ICES Journal of Marine Science* 46 (3): 284–86.
<https://doi.org/10.1093/icesjms/46.3.284>.

Macaulay, Gavin J. 2025. *gavinmacaulay/SphereTS: V1.0.8*.
<https://github.com/gavinmacaulay/SphereTS>.

Macaulay, Gavin, and contributors. 2024. “echoSMs: Making Acoustic
Scattering Models Available to Fisheries and Plankton Scientists.” In
*GitHub Repository*.
[Https://github.com/ices-tools-dev/echoSMs](https://github.com/ices-tools-dev/echoSMs);
GitHub.

MacLennan, D. N. 1981. *The Theory of Solid Spheres as Sonar Calibration
Targets*. Scottish Fisheries Research Report 22. Department of
Agriculture; Fisheries for Scotland.
