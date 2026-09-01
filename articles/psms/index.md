# PSMS

## Overview

Benchmarked Validated

[Theory](https://brandynlucca.github.io/acousticTS/articles/psms/psms-theory.md)
[Implementation](https://brandynlucca.github.io/acousticTS/articles/psms/psms-implementation.md)

The prolate-spheroidal modal-series solution (`PSMS`) calculates
scattering from a homogeneous prolate spheroid. It is the
geometry-matched reference for smooth elongated canonical targets.

### Core idea

Separate the Helmholtz equation in prolate spheroidal coordinates,
expand the fields in angular and radial spheroidal wave functions, and
solve the retained boundary system ([Spence and Granger
1951](#ref-Spence_1951); [Furusawa 1988](#ref-Furusawa_1988)).

### Best for

- Fixed-rigid, pressure-release, liquid-filled, or gas-filled prolate
  spheroids
- Canonical elongated-body benchmarks
- Reference comparisons for the prolate branch of `TMM`

### Supports

- `ProlateSpheroid` shapes carried by `FLS` or `GAS` scatterers
- Broadside and oblique monostatic calculations with explicit roll angle
- Configurable truncation, integration, adaptive evaluation, and
  numerical precision

### Main assumptions

- An exact prolate-spheroidal outer boundary
- One homogeneous interior region
- Linear time-harmonic acoustics
- No shell, secondary component, or arbitrary body profile

### Validation status

- Benchmarked against the canonical prolate-spheroid spectra stored in
  benchmark_ts.
- Validated against the external Prol_Spheroid and echoSMs
  implementations on shared prolate cases.

### Family pages

- [Implementation](https://brandynlucca.github.io/acousticTS/articles/psms/psms-implementation.md):
  supported branches, numerical controls, output, and benchmarks
- [Theory](https://brandynlucca.github.io/acousticTS/articles/psms/psms-theory.md):
  spheroidal coordinates, wave functions, and boundary systems

## References

Furusawa, Masahiko. 1988. “Prolate Spheroidal Models for Predicting
General Trends of Fish Target Strength.” *Journal of the Acoustical
Society of Japan (E)* 9 (1): 13–24. <https://doi.org/10.1250/ast.9.13>.

Spence, R. D., and Sara Granger. 1951. “The Scattering of Sound from a
Prolate Spheroid.” *The Journal of the Acoustical Society of America* 23
(6): 701–6. <https://doi.org/10.1121/1.1906827>.
