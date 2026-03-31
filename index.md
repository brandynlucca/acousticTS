# acousticTS

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7600659.svg)](https://doi.org/10.5281/zenodo.7600659)
[![Documentation](https://img.shields.io/badge/Latest_Documentation-blue)](https://brandynlucca.github.io/acousticTS/)
[![Build
status](https://github.com/brandynlucca/acousticTS/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/brandynlucca/acousticTS/actions/workflows/R-CMD-check.yaml)
[![Coverage](https://codecov.io/gh/brandynlucca/acousticTS/graph/badge.svg?branch=main)](https://app.codecov.io/gh/brandynlucca/acousticTS?branch=main)
[![License:
GPL-3](https://img.shields.io/badge/License-GPL--3-blue.svg)](https://brandynlucca.github.io/acousticTS/LICENSE)

## Overview

acousticTS is an R package for estimating the acoustic target strength
(TS) of aquatic organisms and objects using physics-based scattering
models. It supports a range of scatterer types: fluid-like bodies (e.g.,
fish, zooplankton), gas-filled bodies (e.g. swimbladders), elastic
shells (e.g. pteropods, euphausiids), and calibration spheres. It
further provides a unified interface for parameterizing, running, and
comparing models across frequencies, orientations, and morphologies.
Acoustic backscatter from a single target is expressed as the
*backscattering cross-section* (\sigma\_\text{bs}, m²). Target strength
(TS, dB re. 1 m²) is its logarithmic form:

TS = 10 \log\_{10}(\sigma\_\text{bs})

TS is used to:

- Convert integrated backscatter (e.g. S\_\mathrm{A}, NASC) or
  volumetric backscatter (S\_\text{V}) into number density or biomass
- Classify backscatter by species or taxon based on multi-frequency
  response
- Parameterize and evaluate physics-based scattering models over
  statistical distributions of organism size and orientation

## Installation

Install the latest release from GitHub:

``` r
# install.packages("devtools")
devtools::install_github("brandynlucca/acousticTS")
```

## Scatterer classes

acousticTS organizes targets into five `S4` classes:

| Class | Description               | Example taxa                            |
|-------|---------------------------|-----------------------------------------|
| `FLS` | Fluid-like scatterer      | Euphausiids, myctophids, decapod shrimp |
| `SBF` | Swimbladder-bearing fish  | Herring, cod, sardine                   |
| `GAS` | Gas-filled body           | Siphonophore pneumatophores             |
| `ESS` | Elastic-shelled scatterer | Pteropods, juvenile bivalves            |
| `CAL` | Calibration sphere        | Tungsten carbide, copper spheres        |

Each class stores a shape (position matrix + morphometrics), body
material properties (density and sound speed contrasts), and model
results in a structured `S4` object.

## Models

| Model                                    | Abbreviation | Scatterer type | Boundary      |
|------------------------------------------|--------------|----------------|---------------|
| Distorted-wave Born approximation        | `DWBA`       | FLS            | Fluid         |
| Stochastic DWBA                          | `SDWBA`      | FLS            | Fluid         |
| Kirchhoff-ray mode                       | `KRM`        | FLS, SBF       | Fluid + gas   |
| High-pass approximation                  | `HPA`        | FLS, GAS       | Fluid / gas   |
| Two-ray cylinder model                   | `TRCM`       | FLS            | Fluid         |
| Modal series solution (sphere)           | `SPHMS`      | CAL, ESS       | Multiple      |
| Modal series solution (prolate spheroid) | `PSMS`       | FLS            | Fluid         |
| Modal series solution (elastic shell)    | `ESSMS`      | ESS            | Elastic       |
| Finite-cylinder modal series             | `FCMS`       | CAL            | Rigid / fluid |
| Resonance model (gas sphere)             | `SOEMS`      | GAS, CAL       | Gas           |

## Quick start

``` r
library(acousticTS)

# Build a fluid-like scatterer (e.g. krill) with a cylinder shape
krill_shape <- cylinder(
  length_body = 0.03,          # 30 mm body length
  radius_body = 0.003          # 3 mm max radius
)
# ---- Create the actual Scatterer-class object (FLS)
krill <- fls_generate(
  shape = krill_shape,
  g = 1.0357,             # density contrast
  h = 1.0279,             # sound speed contrast
  theta = pi / 2          # broadside incidence
)

# Run the DWBA model from 1 kHz to 300 kHz
krill <- target_strength(krill,
  frequency = seq(1e3, 300e3, by = 1e3),
  model = "DWBA"
)

# Plot TS vs frequency
plot(krill, type = "model")
```

## Shapes and reforging

Scatterer shapes are created via dedicated constructors:

``` r
# Cylinder (used by FLS / CAL)
fish_shape <- cylinder(length_body = 0.25, radius_body = 0.02)

# Prolate spheroid (used by FLS / PSMS)
ps_shape <- prolate_spheroid(length_body = 0.02, radius_body = 0.002)

# Arbitrary shape from digitized position matrix
arb_shape <- arbitrary(rpos = my_matrix)

# Dorsal/ventral style inputs
arbitrary(
  x_body = c(0, 0.015),
  w_body = c(0.005, 0.0075),
  zU_body = c(0.001, 0.002),
  zL_body = c(-0.001, -0.002)
)
```

Existing scatterer objects can be resized or re-discretized without
reconstructing them from scratch using
[`reforge()`](https://brandynlucca.github.io/acousticTS/reference/reforge.md):

``` r
# Rescale a krill to 40 mm
krill_40mm <- reforge(krill, body_target = c(length = 0.04))

# Change segment count
krill_fine <- reforge(krill, n_segments_body = 200)
```

## Configuring simulations

[`simulate_ts()`](https://brandynlucca.github.io/acousticTS/reference/simulate_ts.md)
runs repeated model evaluations across distributions of input
parameters, supporting both vectorized and batch modes:

``` r
results <- simulate_ts(
  krill,
  model = "DWBA",
  frequency = 120e3,
  n_realizations = 1000,
  parameters = list(
    theta_body = function() rnorm(1, pi / 2, pi / 18),
    length_body = function() rnorm(1, 0.03, 0.003)
  )
)
```

## Built-in datasets

| Dataset        | Description                                                         |
|----------------|---------------------------------------------------------------------|
| `krill`        | Antarctic krill (*Euphausia superba*) shape and material properties |
| `cod`          | Atlantic cod (*Gadus morhua*) shape and swimbladder                 |
| `sardine`      | Pacific sardine (*Sardinops sagax*) shape and swimbladder           |
| `benchmark_ts` | Benchmark TS values for model validation                            |

``` r
data(krill)
data(cod)
data(sardine)
data(benchmark_ts)
```

## Citation

If you use acousticTS in published work, please cite:

``` r
citation("acousticTS")
```

The Zenodo concept DOI for acousticTS, which resolves to the latest
archived release record, is:

<https://doi.org/10.5281/zenodo.7600659>

## Contributing and bug reports

Bug reports and feature requests are welcome via the [GitHub Issues
page](https://github.com/brandynlucca/acousticTS/issues). Please include
a minimal reproducible example when reporting bugs.

## License

GPL-3. See [LICENSE](https://brandynlucca.github.io/acousticTS/LICENSE)
for details.
