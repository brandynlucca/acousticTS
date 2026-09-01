# FCMS

## Overview

Benchmarked Validated

[Theory](https://brandynlucca.github.io/acousticTS/articles/fcms/fcms-theory.md)
[Implementation](https://brandynlucca.github.io/acousticTS/articles/fcms/fcms-implementation.md)

The finite-cylinder modal-series solution (`FCMS`) calculates monostatic
scattering from a straight circular cylinder near broadside.

### Core idea

Solve the circular cross-section in cylindrical partial waves, then
apply the finite-length coherence factor that describes addition along
the straight cylinder axis ([Stanton 1988](#ref-Stanton_1988),
[1989](#ref-Stanton_1989_2)).

### Best for

- Straight, constant-radius finite cylinders near broadside
- Comparing fixed-rigid, pressure-release, liquid-filled, and gas-filled
  cylinder boundaries
- Canonical cylinder benchmarks before adopting an asymptotic model

### Supports

- `Cylinder` shapes carried by `FLS` or `GAS` scatterers
- Homogeneous fluid or gas interiors, plus fixed-rigid and
  pressure-release limits
- Monostatic complex amplitude, backscattering cross-section, and target
  strength

### Main assumptions

- A straight centerline and circular, constant-radius cross-section
- Broadside or near-broadside incidence
- Finite length represented by a coherence factor
- No elastic shear response in the cylinder interior

### Validation status

- Benchmarked against the canonical finite-cylinder spectra stored in
  benchmark_ts.
- Validated against the echoSMs finite-cylinder implementation.

### Family pages

- [Implementation](https://brandynlucca.github.io/acousticTS/articles/fcms/fcms-implementation.md):
  cylinder construction, boundaries, output, and comparisons
- [Theory](https://brandynlucca.github.io/acousticTS/articles/fcms/fcms-theory.md):
  cylindrical modes, interface conditions, and finite-length closure

## References

Stanton, T. K. 1988. “Sound Scattering by Cylinders of Finite Length. I.
Fluid Cylinders.” *The Journal of the Acoustical Society of America* 83
(1): 55–63. <https://doi.org/10.1121/1.396184>.

Stanton, T. K. 1989. “Sound Scattering by Cylinders of Finite Length.
III. Deformed Cylinders.” *The Journal of the Acoustical Society of
America* 86 (2): 691–705. <https://doi.org/10.1121/1.398193>.
