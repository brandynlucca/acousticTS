# Model library

## Introduction

The family library is easiest to navigate when read against
model-comparison reviews, organism-focused scattering surveys, and
companion open-source modeling ecosystems ([Jech et al.
2015](#ref-jech_etal_2015); [Stanton 1996](#ref-stanton_acoustic_1996);
[Lucca and Lee 2026](#ref-brandyn_lucca_osoceanacousticsechopop_2026)).

The package contains enough model families that the site works better
when the model pages are treated as a library rather than as a flat list
in the top navigation. This page is the main entry point for that
library.

### Model-status policy

- Benchmarked means the family has a documented comparison against a
  canonical benchmark ladder or stored benchmark values.
- Validated means the currently supported public scope has a documented
  external software or independent-comparison check.
- Partially validated means some supported branches are externally
  checked, but the full public scope is not yet closed.
- Unvalidated means the package does not yet claim external validation
  across the current public scope.
- Experimental means the family is available to use, but its interface
  or supported workflow should still be treated as provisional.

These tags are intended to be read in three pieces:

- `Benchmarked` is independent of the validation badge and can appear
  alongside `Validated`, `Partially validated`, or `Unvalidated`.
- The validation badge is always exactly one of `Validated`,
  `Partially validated`, or `Unvalidated`.
- `Experimental` is a separate lifecycle tag and can coexist with any
  benchmark or validation badge.

## Modal-series families

### [SPHMS](https://brandynlucca.github.io/acousticTS/articles/sphms/index.md)

Benchmarked Validated

Spherical modal-series solution for canonical spherical targets.

### [FCMS](https://brandynlucca.github.io/acousticTS/articles/fcms/index.md)

Benchmarked Validated

Finite-cylinder modal-series solution for straight cylindrical targets.

### [PSMS](https://brandynlucca.github.io/acousticTS/articles/psms/index.md)

Benchmarked Validated

Prolate-spheroidal modal-series solution for smooth elongated canonical
bodies.

### [SOEMS](https://brandynlucca.github.io/acousticTS/articles/calibration/index.md)

Benchmarked Validated

Solid elastic spherical model used mainly for calibration spheres.

### [ESSMS](https://brandynlucca.github.io/acousticTS/articles/essms/index.md)

Unvalidated

Elastic-shelled spherical family for layered shell targets.

### [BCMS](https://brandynlucca.github.io/acousticTS/articles/bcms/index.md)

Unvalidated Experimental

Bent-cylinder modal-series family for straight and uniformly bent
cylinders.

### [ECMS](https://brandynlucca.github.io/acousticTS/articles/ecms/index.md)

Unvalidated Experimental

Elastic-cylinder modal-series family for fully elastic solid cylinders.

## Approximation and ray-based families

### [DWBA](https://brandynlucca.github.io/acousticTS/articles/dwba/index.md)

Benchmarked Validated

Weak-scattering elongated-body approximation for fluid-like targets.

### [SDWBA](https://brandynlucca.github.io/acousticTS/articles/sdwba/index.md)

Benchmarked Validated

Stochastic DWBA family for unresolved phase variability.

### [KRM](https://brandynlucca.github.io/acousticTS/articles/krm/index.md)

Benchmarked Validated

Kirchhoff-ray mode model for segmented fish-like body-plus-inclusion
targets.

### [HPA](https://brandynlucca.github.io/acousticTS/articles/hpa/index.md)

Benchmarked Validated

High-pass approximation for compact asymptotic screening.

### [TRCM](https://brandynlucca.github.io/acousticTS/articles/trcm/index.md)

Benchmarked Unvalidated

Two-ray cylindrical family for high-frequency locally cylindrical
targets.

### [PCDWBA](https://brandynlucca.github.io/acousticTS/articles/pcdwba/index.md)

Validated Experimental

Phase-compensated DWBA for bent weakly scattering targets.

## Composite and emerging families

### [BBFM](https://brandynlucca.github.io/acousticTS/articles/bbfm/index.md)

Unvalidated Experimental

Composite flesh-plus-backbone family for swimbladder-less fish-like
targets.

### [VESM](https://brandynlucca.github.io/acousticTS/articles/vesm/index.md)

Validated Experimental

Viscous-elastic layered-sphere family for gas-core, shell, and
viscous-layer targets.

### [TMM](https://brandynlucca.github.io/acousticTS/articles/tmm/index.md)

Benchmarked Partially validated Experimental

Single-target transition-matrix family for retained monostatic and
angle-dependent scattering products across supported canonical shapes.

## How to use this page

Each family page is intended to be the local navigation hub for that
model family. Use the family page to move between overview, theory,
implementation, and related package pages, rather than relying on a long
global theory or implementation menu.

## References

Jech, J. Michael, John K. Horne, Dezhang Chu, David A. Demer, David T.
I. Francis, Natalia Gorska, Benjamin Jones, et al. 2015. “Comparisons
Among Ten Models of Acoustic Backscattering Used in Aquatic Ecosystem
Research.” *The Journal of the Acoustical Society of America* 138 (6):
3742–64. <https://doi.org/10.1121/1.4937607>.

Lucca, Brandyn, and Wu-Jung Lee. 2026. “OSOceanAcoustics/Echopop:
V0.6.0.” Zenodo. <https://doi.org/10.5281/ZENODO.18975959>.

Stanton, T. 1996. “Acoustic Scattering Characteristics of Several
Zooplankton Groups.” *ICES Journal of Marine Science* 53 (2): 289–95.
<https://doi.org/10.1006/jmsc.1996.0037>.
