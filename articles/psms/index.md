# PSMS

## Overview

Benchmarked Validated

*Model-family pages:*
[Overview](https://brandynlucca.github.io/acousticTS/articles/psms/index.md)
[Implementation](https://brandynlucca.github.io/acousticTS/articles/psms/psms-implementation.md)
[Theory](https://brandynlucca.github.io/acousticTS/articles/psms/psms-theory.md)

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
