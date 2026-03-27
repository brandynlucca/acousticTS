# TMM

## Overview

Benchmarked Validated Experimental

*Model-family pages:*
[Overview](https://brandynlucca.github.io/acousticTS/articles/tmm/index.md)
[Implementation](https://brandynlucca.github.io/acousticTS/articles/tmm/tmm-implementation.md)
[Theory](https://brandynlucca.github.io/acousticTS/articles/tmm/tmm-theory.md)

These pages follow the coefficient-map view of scattering and later
numerical implementations for axisymmetric bodies ([Waterman
1969](#ref-waterman_new_1969), [2009](#ref-waterman_t_2009); [Ganesh and
Hawkins 2022](#ref-ganesh_numerically_2022)).

The transition matrix method (`TMM`) is the package’s current
single-target bridge between exact modal-series solvers and broader
angle-dependent scattering products.

### Core idea

Represent the incident and scattered fields in complete modal bases and
solve for the linear map between those coefficient vectors. In the
current package scope, that gives a reusable single-target retained
state for monostatic target strength and, where externally constrained,
for general-angle or orientation-averaged post-processing.

### Best for

- Single-target axisymmetric scattering problems that need more than one
  post-processed product from one solve
- Sphere, oblate spheroid, and prolate spheroid problems with retained
  scattering products
- Geometry-specific comparisons between exact modal families and a
  T-matrix viewpoint

### Supports

- `Sphere`, `OblateSpheroid`, `ProlateSpheroid`, and guarded `Cylinder`
  branches
- Single homogeneous interiors under rigid, pressure-release,
  liquid-filled, and gas-filled boundaries
- Stored retained state for scattering slices, grids, diagnostics, and
  orientation averages where validated

### Main assumptions

- Single-target scope only
- Geometry-specific basis choice rather than one universal retained
  operator for every shape
- Cylinder branch currently has narrower validated scope than sphere,
  oblate, and prolate branches

### Validation status

- Benchmarked against `SPHMS`, `PSMS`, and `FCMS` on the currently
  supported canonical shape branches.
- Validated against external BEMPP far-field checks for sphere, oblate,
  and prolate pressure-release cases.
- Retained prolate angular products are also checked against the exact
  general-angle spheroidal solution.
- The cylinder branch is benchmark-matched only for the exact monostatic
  workflow; retained general-angle cylinder products remain outside the
  validated public scope.

### Family pages

- [Implementation](https://brandynlucca.github.io/acousticTS/articles/tmm/tmm-implementation.md):
  stored-state workflows, plots, benchmarks, and current scope
- [Theory](https://brandynlucca.github.io/acousticTS/articles/tmm/tmm-theory.md):
  T-matrix interpretation, boundary operators, and geometry-matched
  bases

![Single-target TMM workflow, showing the spherical-coordinate branch
used for spheres and oblates, the spheroidal branch used for prolates,
and the cylindrical monostatic branch used for finite
cylinders.](tmm-branch-schematic.svg)

Single-target TMM workflow, showing the spherical-coordinate branch used
for spheres and oblates, the spheroidal branch used for prolates, and
the cylindrical monostatic branch used for finite cylinders.

## References

Ganesh, M., and Stuart C. Hawkins. 2022. “A Numerically Stable T-Matrix
Method for Acoustic Scattering by Nonspherical Particles with Large
Aspect Ratios and Size Parameters.” *The Journal of the Acoustical
Society of America* 151 (3): 1978–88.
<https://doi.org/10.1121/10.0009679>.

Waterman, P. C. 1969. “New Formulation of Acoustic Scattering.” *The
Journal of the Acoustical Society of America* 45 (6): 1417–29.
<https://doi.org/10.1121/1.1911619>.

———. 2009. “T -Matrix Methods in Acoustic Scattering.” *The Journal of
the Acoustical Society of America* 125 (1): 42–51.
<https://doi.org/10.1121/1.3035839>.
