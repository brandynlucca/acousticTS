# TMM

## Overview

Benchmarked Partially validated Experimental

[Theory](https://brandynlucca.github.io/acousticTS/articles/tmm/tmm-theory.md)
[Implementation](https://brandynlucca.github.io/acousticTS/articles/tmm/tmm-implementation.md)

The transition-matrix method (`TMM`) represents a target by the linear
map from incident modal coefficients to scattered modal coefficients.
Its retained state can support more than a monostatic target-strength
curve.

### Core idea

Choose a basis matched to the supported geometry, solve the boundary
problem for the coefficient map, and optionally retain its
frequency-wise blocks for angular scattering, orientation averages, and
diagnostics ([Waterman 1969](#ref-Waterman_1969),
[2009](#ref-Waterman_2009)).

### Best for

- Supported single-target canonical geometries requiring retained
  angular products
- Bistatic slices, scattering grids, and orientation post-processing
  from one solved state
- Cross-checking geometry-specific modal families within a
  transition-matrix representation

### Supports

- Sphere, oblate-spheroid, prolate-spheroid, and guarded finite-cylinder
  branches
- Supported fixed-rigid, pressure-release, penetrable, and
  spherical-shell boundaries
- Optional retained T-matrix blocks for branch-specific post-processing

### Main assumptions

- A single target represented by one of the documented geometry-specific
  bases
- Homogeneous-region or supported concentric spherical-shell material
  structure
- Post-processing uses the documented body-fixed angular convention
- Validation and retained-state support differ by geometry and boundary
  branch

### Validation status

- Benchmarked against `SPHMS`, `PSMS`, and `FCMS` on the currently
  supported canonical shape branches.
- Validated against external BEMPP far-field checks for sphere, oblate,
  and prolate pressure-release cases.
- Retained prolate angular products are also checked against the exact
  general-angle spheroidal solution.
- TMM is partially validated because the sphere, oblate, and prolate
  branches have external checks, but retained general-angle cylinder
  products remain outside the validated public scope.
- TMM is currently marked experimental because the retained-state
  workflow and branch matrix are still guarded while shape-specific
  support continues to be tightened.

### Family pages

- [Implementation](https://brandynlucca.github.io/acousticTS/articles/tmm/tmm-implementation.md):
  branch matrix, retained-state workflows, output, and validation
- [Theory](https://brandynlucca.github.io/acousticTS/articles/tmm/tmm-theory.md):
  coefficient maps, boundary operators, and geometry-matched bases

## References

Waterman, P. C. 1969. “New Formulation of Acoustic Scattering.” *The
Journal of the Acoustical Society of America* 45 (6): 1417–29.
<https://doi.org/10.1121/1.1911619>.

Waterman, P. C. 2009. “T -Matrix Methods in Acoustic Scattering.” *The
Journal of the Acoustical Society of America* 125 (1): 42–51.
<https://doi.org/10.1121/1.3035839>.
