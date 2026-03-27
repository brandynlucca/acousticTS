# FCMS

## Overview

Benchmarked Validated

*Model-family pages:*
[Overview](https://brandynlucca.github.io/acousticTS/articles/fcms/index.md)
[Implementation](https://brandynlucca.github.io/acousticTS/articles/fcms/fcms-implementation.md)
[Theory](https://brandynlucca.github.io/acousticTS/articles/fcms/fcms-theory.md)

These pages follow the finite-cylinder modal-series literature for
straight circular cylinders near broadside ([Stanton
1988](#ref-stanton_sound_1988), [1989](#ref-stanton_sound_1989)).

The finite cylinder modal series solution (`FCMS`) is the package’s
geometry-matched cylinder family for straight circular cylinders. It
keeps the exact cylindrical-harmonic treatment of the cross-section and
closes the finite-length problem with the standard near-broadside
coherence factor.

### Core idea

Solve the circular cross-section exactly in cylindrical partial waves,
then multiply by the finite-length directivity factor that accounts for
the coherent addition along the cylinder axis.

### Best for

- Straight finite cylinders near broadside
- Rigid, pressure-release, liquid-filled, and gas-filled cylinder
  comparisons
- Reference checks for cylinder-like targets before using higher-level
  approximations

### Supports

- `Cylinder` shapes on `FLS` or `GAS` objects
- Near-broadside monostatic target strength
- Fluid and gas interiors parameterized relative to seawater as medium
  `1`

### Main assumptions

- Circular cross-section and straight centerline
- Broadside or near-broadside incidence
- Homogeneous interior material properties for penetrable cases
- No elastic shear support in the interior

### Validation status

- Benchmarked against the canonical finite-cylinder spectra stored in
  `benchmark_ts`.
- Validated against the `echoSMs` finite-cylinder implementation.

### Family pages

- [Implementation](https://brandynlucca.github.io/acousticTS/articles/fcms/fcms-implementation.md):
  spectra, result extraction, and comparison tables
- [Theory](https://brandynlucca.github.io/acousticTS/articles/fcms/fcms-theory.md):
  cylindrical modal reduction, boundary conditions, and finite-length
  factorization

## References

Stanton, T. K. 1988. “Sound Scattering by Cylinders of Finite Length. I.
Fluid Cylinders.” *The Journal of the Acoustical Society of America* 83
(1): 55–63. <https://doi.org/10.1121/1.396184>.

———. 1989. “Sound Scattering by Cylinders of Finite Length. III.
Deformed Cylinders.” *The Journal of the Acoustical Society of America*
86 (2): 691–705. <https://doi.org/10.1121/1.398193>.
