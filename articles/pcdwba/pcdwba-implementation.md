# Phase-compensated distorted wave Born approximation

## acousticTS implementation

Validated Experimental

[Overview](https://brandynlucca.github.io/acousticTS/articles/pcdwba/index.md)
[Theory](https://brandynlucca.github.io/acousticTS/articles/pcdwba/pcdwba-theory.md)

These pages follow the phase-compensated weak-scattering literature for
broadside elongated bodies and krill-style applications ([Chu and Ye
1999](#ref-Chu_1999); [Chu, Foote, and Stanton 1993](#ref-Chu_1993)).

The phase-compensated distorted wave Born approximation is available
through `target_strength(..., model = "pcdwba")`. The implementation is
intended for weakly scattering fluid-like bodies and uses the same
curved-cylinder bookkeeping whether the target starts as a canonical
bent cylinder or an arbitrary fluid-like profile.

This page checks the implementation against two source-level references:

- the `pcdwba_fbs` routine in the `Python` package Echopop (Lucca and
  Lee ([2026](#ref-Echopop_software))),
- the bent-cylinder DWBA routines in the `R`-package ZooScatR (Gastauer,
  Chu, and Cox ([2019](#ref-ZooScatR_software))) .

`PCDWBA` is validated here against source-level reference
implementations rather than against a separate published benchmark
table. The ZooScatR source agrees exactly on the shared case, while the
remaining Echopop drift is attributable to that implementation’s
interpolated Bessel evaluation.

### Reference case

The comparison uses a single reproducible bent-cylinder case:

- length `15 mm`
- radius `1 mm`
- taper order `10`
- curvature ratio `rho_c / L = 3`
- density contrast `g = 1.02`
- sound-speed contrast `h = 1.02`
- broadside incidence
- `12-200 kHz` in `2 kHz` steps
- `51` integration nodes in all three implementations

In acousticTS, that target is built as:

``` r
library(acousticTS)

pcdwba_object <- fls_generate(
  shape = cylinder(
    length_body = 0.015,
    radius_body = 0.001,
    taper = 10,
    radius_curvature_ratio = 3,
    n_segments = 50
  ),
  g_body = 1.02,
  h_body = 1.02,
  theta_body = pi / 2
)

pcdwba_object <- target_strength(
  object = pcdwba_object,
  frequency = seq(12e3, 200e3, by = 2e3),
  model = "pcdwba",
  sound_speed_sw = 1500,
  density_sw = 1026
)

head(extract(pcdwba_object, "model")$PCDWBA)
```

    ##   frequency         ka                        f_bs     sigma_bs        TS
    ## 1     12000 0.05026548 -6.147156e-07+7.354948e-11i 3.778753e-13 -124.2265
    ## 2     14000 0.05864306 -8.131273e-07+1.148885e-10i 6.611760e-13 -121.7968
    ## 3     16000 0.06702064 -1.027176e-06+1.682542e-10i 1.055090e-12 -119.7671
    ## 4     18000 0.07539822 -1.251059e-06+2.344088e-10i 1.565149e-12 -118.0544
    ## 5     20000 0.08377580 -1.478570e-06+3.137691e-10i 2.186169e-12 -116.6032
    ## 6     22000 0.09215338 -1.703221e-06+4.063858e-10i 2.900962e-12 -115.3746

### Validation outputs

#### Comparison summary

| Comparison                    | Max abs. \Delta TS (dB) | Mean abs. \Delta TS (dB) |
|:------------------------------|------------------------:|-------------------------:|
| acousticTS vs echopop         |                0.073947 |                 0.001123 |
| acousticTS vs ZooScatR-source |                0.000000 |                 0.000000 |
| echopop vs ZooScatR-source    |                0.073947 |                 0.001123 |

The ZooScatR and acousticTS outputs are indistinguishable on this grid.
The Echopop comparison remains close as well, but it is not at machine
precision because that implementation evaluates the cylindrical Bessel
term through interpolation rather than a direct nodewise call. On this
grid, the largest mismatch occurs near `112 kHz`; replacing the
interpolated `J_1(x)/x` evaluation with a direct call collapses that
residual onto the acousticTS / ZooScatR curve. So the remaining drift is
numerical, not geometrical.

#### Timings

| Implementation  | Elapsed (s) | f min (kHz) | f max (kHz) | Step (kHz) | Node count |
|:----------------|------------:|------------:|------------:|-----------:|-----------:|
| acousticTS      |      0.0100 |          12 |         200 |          2 |         51 |
| echopop         |      3.3122 |          12 |         200 |          2 |         51 |
| ZooScatR-source |      0.1200 |          12 |         200 |          2 |         51 |

#### Spectrum overlay

![Pre-rendered PCDWBA comparison showing ZooScatR, acousticTS, and
echopop spectra together with the acousticTS residuals against the two
references.](pcdwba-reference-comparison.png)

### Closing note

This is the kind of implementation check that matters for a
phase-compensated bent-cylinder solver. The comparison is not just
against a benchmark curve. It is against two independently written
source routines that share the same governing model. On this reference
case, acousticTS reproduces the direct ZooScatR-style calculation
exactly and stays very close to the Echopop implementation across the
full frequency band.

## References

Chu, Dezhang, Kenneth G. Foote, and Timothy K. Stanton. 1993. “Further
Analysis of Target Strength Measurements of Antarctic Krill at 38 and
120 kHz: Comparison with Deformed Cylinder Model and Inference of
Orientation Distribution.” *The Journal of the Acoustical Society of
America* 93 (5): 2985–88. <https://doi.org/10.1121/1.405818>.

Chu, Dezhang, and Zhen Ye. 1999. “A Phase-Compensated Distorted Wave
Born Approximation Representation of the Bistatic Scattering by Weakly
Scattering Objects: Application to Zooplankton.” *The Journal of the
Acoustical Society of America* 106 (4): 1732–43.
<https://doi.org/10.1121/1.428036>.

Gastauer, Sven, Dezhang Chu, and Martin J. Cox. 2019. “ZooScatR—An
\<Span Style="font-Variant:small-Caps;"\>r\</Span\> Package for
Modelling the Scattering Properties of Weak Scattering Targets Using the
Distorted Wave Born Approximation.” *The Journal of the Acoustical
Society of America* 145 (1): EL102–8.
<https://doi.org/10.1121/1.5085655>.

Lucca, Brandyn, and Wu-Jung Lee. 2026. “OSOceanAcoustics/Echopop:
V0.6.0.” Zenodo. <https://doi.org/10.5281/ZENODO.18975959>.
