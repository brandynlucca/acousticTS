# PSMS

## Overview

Benchmarked Validated

[Theory](https://brandynlucca.github.io/acousticTS/articles/psms/psms-theory.md)
[Implementation](https://brandynlucca.github.io/acousticTS/articles/psms/psms-implementation.md)

These pages are rooted in exact spheroidal-coordinate separations and
later fisheries-acoustics use of prolate-spheroid models ([Spence and
Granger 1951](#ref-Spence_1951); [Furusawa 1988](#ref-Furusawa_1988)).

The prolate spheroidal modal series solution (`PSMS`) is the exact
single-target modal family for homogeneous prolate spheroids. It is the
spheroidal analogue of spherical partial-wave theory and the natural
exact reference for elongated canonical bodies whose surface follows a
prolate spheroid.

### Core idea

Separate the Helmholtz equation in prolate spheroidal coordinates,
expand the incident, scattered, and interior fields in spheroidal wave
functions, and solve the retained boundary systems order by order.

### Best for

- Rigid, pressure-release, liquid-filled, and gas-filled prolate
  spheroids
- Canonical elongated-body benchmarks
- Validating prolate branches of more general methods such as `TMM`

### Supports

- `ProlateSpheroid` shapes on `FLS` or `GAS` objects
- Exact monostatic target strength for homogeneous single-region
  spheroids
- Spheroidal modal solutions in which medium `1` is seawater and medium
  `2` is the spheroid interior

### Main assumptions

- Perfect prolate spheroidal geometry
- Homogeneous interior region
- Linear, time-harmonic acoustics
- No shell or internal secondary component

### Validation status

- Benchmarked against the canonical prolate-spheroid spectra stored in
  `benchmark_ts`.
- Validated against the external `Prol_Spheroid` implementation on
  shared prolate cases.

### Family pages

- [Implementation](https://brandynlucca.github.io/acousticTS/articles/psms/psms-implementation.md):
  workflows, comparisons, and timing tables
- [Theory](https://brandynlucca.github.io/acousticTS/articles/psms/psms-theory.md):
  full spheroidal-coordinate derivation and retained modal systems

## References

Furusawa, Masahiko. 1988. “Prolate Spheroidal Models for Predicting
General Trends of Fish Target Strength.” *Journal of the Acoustical
Society of Japan (E)* 9 (1): 13–24. <https://doi.org/10.1250/ast.9.13>.

Spence, R. D., and Sara Granger. 1951. “The Scattering of Sound from a
Prolate Spheroid.” *The Journal of the Acoustical Society of America* 23
(6): 701–6. <https://doi.org/10.1121/1.1906827>.
