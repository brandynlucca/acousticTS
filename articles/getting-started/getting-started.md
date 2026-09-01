# Quick Start

## What acousticTS does

`acousticTS` estimates the acoustic target strength of aquatic organisms
and other objects with physics-based scattering models. It provides a
common interface for defining target geometry and material properties,
running exact or approximate models, and comparing responses across
frequency, orientation, and morphology.

These pages explain how the package represents and computes its
scattering models. They provide the context needed to use and interpret
`acousticTS`, but they do not replace comprehensive acoustics texts or
the original model literature.

## Installation

`acousticTS` requires R 4.0.0 or later. Install the current GitHub
version from source with:

``` r

install.packages("remotes")

remotes::install_github(
  "brandynlucca/acousticTS",
  build_vignettes = FALSE
)
```

Load the package with:

``` r

library(acousticTS)
```

`acousticTS` contains C++17 and Fortran code, so installation from
GitHub requires a working compilation toolchain.

- **Windows:** install the version of
  [Rtools](https://cran.r-project.org/bin/windows/Rtools/) that matches
  your R version.
- **macOS:** install the Xcode Command Line Tools and the GNU Fortran
  compiler recommended for your R release on the [R for macOS tools
  page](https://mac.r-project.org/tools/). Xcode alone does not provide
  a Fortran compiler.
- **Linux:** install `make`, a C++17 compiler, `gfortran`, and the
  development libraries used by your R installation. Package names vary
  by distribution.

After installing or changing compilers, restart R before retrying the
package installation. Errors mentioning `CXX17`, `gfortran`, `quadmath`,
BLAS, or LAPACK usually indicate a missing or mismatched toolchain. If
the problem continues, include
[`sessionInfo()`](https://rdrr.io/r/utils/sessionInfo.html) and the
complete installation log in a [GitHub
issue](https://github.com/brandynlucca/acousticTS/issues).

## The package workflow

The package is organized around three steps: build a geometric `Shape`,
wrap that shape in a physically defined scatterer, and run one or more
compatible models.

![Getting started
workflow](getting-started-workflow.png)[](https://brandynlucca.github.io/acousticTS/articles/building-shapes/building-shapes.md "Build a shape")[](https://brandynlucca.github.io/acousticTS/articles/building-scatterers/building-scatterers.md "Wrap as scatterer")[](https://brandynlucca.github.io/acousticTS/articles/running-models/running-models.md "Run model(s)")

Select **Build a shape**, **Wrap as scatterer**, or **Run model(s)** in
the figure to open the corresponding guide.

## A minimal example

This example creates a small fluid-like cylinder, assigns its acoustic
contrasts, evaluates a frequency response, and plots the stored result:

``` r

library(acousticTS)

# 1. Build a shape
krill_shape <- cylinder(
  length_body = 0.03,
  radius_body = 0.003
)

# 2. Wrap the shape as a scatterer
krill <- fls_generate(
  shape = krill_shape,
  g_body = 1.0357,
  h_body = 1.0279,
  theta_body = pi / 2
)

# 3. Run a model
krill <- target_strength(
  object = krill,
  frequency = seq(38e3, 120e3, by = 2e3),
  model = "DWBA"
)

plot(krill, type = "model")
```

The `Shape` stores geometry only.
[`fls_generate()`](https://brandynlucca.github.io/acousticTS/reference/fls_generate.md)
adds the material contrasts and orientation needed to interpret that
geometry as a fluid-like target.
[`target_strength()`](https://brandynlucca.github.io/acousticTS/reference/target_strength.md)
then stores the model output on the scatterer for plotting, extraction,
comparison, or simulation.

## Where to go next

Continue with [Building
Shapes](https://brandynlucca.github.io/acousticTS/articles/building-shapes/building-shapes.md),
[Building
Scatterers](https://brandynlucca.github.io/acousticTS/articles/building-scatterers/building-scatterers.md),
or [Running
Models](https://brandynlucca.github.io/acousticTS/articles/running-models/running-models.md)
when one of the three workflow stages needs more detail.
