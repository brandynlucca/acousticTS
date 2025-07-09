
# acousticTS

Acoustic backscatter from a single target or organism is expressed as
the intensity of an echo typically denoted as the *backscattering
cross-section* (σ<sub>bs</sub>, m<sup>2</sup>). Target strength (TS, dB
re. 1 m<sup>2</sup>) is the logarithmic representation of σ<sub>bs</sub>
where: TS = 10 log<sub>10</sub> (σ<sub>bs</sub>). TS can be used to
convert integrated (e.g. nautical area scattering coefficient,
S<sub>A</sub>, dB re. 1(m<sup>2</sup> nmi<sup>-2</sup>) or volumetric
backscatter (e.g. S<sub>v</sub>, dB re. 1 m<sup>-1</sup>) collected from
fisheries acoustic surveys into units of number density, such as the
volumetric density of a fish school (i.e. animals m<sup>-3</sup>). This
parameter can also aid in classifying backscatter based on the
multifrequency response of targets, such as separating likely echoes of
large predatory fish from smaller prey. While there are several
approaches for estimating TS, one common method is to apply
physics-based models to predict theoretical TS that comprise exact and
approximate solutions. The models provided in the `acousticTS` package
can help provide TS estimates parameterized using broad statitsical
distributions of inputs. This package is in a constant state of
development with updates to the available model library, computational
efficiency, and quality-of-life improvements.

*Latest release DOI*
<https://doi.org/10.5281/zenodo.7600659>

*General DOI*
[![DOI](https://zenodo.org/badge/161965429.svg)](https://zenodo.org/badge/latestdoi/161965429)

*Latest document & API*
[![GitHub Pages](https://img.shields.io/badge/github%20pages-121013?style=for-the-badge&logo=github&logoColor=white)](https://brandynlucca.github.io/acousticTS/)

## Installation

You can install the current released version of acousticTS via:

``` r
devtools::install_github("brandynlucca/acousticTS")
```

Or you can install the development version of acousticTS like so:

``` r
devtools::install_github("brandynlucca/acousticTS@test-branch")
```