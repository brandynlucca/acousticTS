# BCMS

## Overview

Unvalidated Experimental

[Theory](https://brandynlucca.github.io/acousticTS/articles/bcms/bcms-theory.md)
[Implementation](https://brandynlucca.github.io/acousticTS/articles/bcms/bcms-implementation.md)

The bent-cylinder modal-series model (`BCMS`) calculates near-broadside
scattering from straight or uniformly bent circular cylinders while
retaining the `FCMS` cross-sectional modal solution.

### Core idea

Use the straight-cylinder partial-wave coefficients locally, then
replace the straight axial coherence factor with the Fresnel-integral
correction for uniform curvature ([Stanton 1989](#ref-Stanton_1989_2)).

### Best for

- Uniformly bent, constant-radius cylinders near broadside
- Isolating curvature-driven coherence changes from cross-sectional
  boundary physics
- Comparisons with the straight `FCMS` limit

### Supports

- `Cylinder` shapes carried by `FLS` or `GAS` scatterers
- Curvature supplied through the shape’s radius-of-curvature ratio
- Liquid-filled, gas-filled, fixed-rigid, and pressure-release
  boundaries

### Main assumptions

- Uniform centerline curvature and constant circular cross-section
- Broadside or near-broadside incidence
- Curvature changes axial coherence but not local modal coefficients
- No arbitrary bending, torsion, taper, or fully three-dimensional
  boundary solve

### Validation status

- BCMS is currently marked experimental because the documented checks
  are internal coherence reconstructions rather than an external
  benchmark or software-comparison ladder.

### Family pages

- [Implementation](https://brandynlucca.github.io/acousticTS/articles/bcms/bcms-implementation.md):
  curvature input, boundaries, output, and internal checks
- [Theory](https://brandynlucca.github.io/acousticTS/articles/bcms/bcms-theory.md):
  Fresnel coherence correction and the straight-cylinder limit

## References

Stanton, T. K. 1989. “Sound Scattering by Cylinders of Finite Length.
III. Deformed Cylinders.” *The Journal of the Acoustical Society of
America* 86 (2): 691–705. <https://doi.org/10.1121/1.398193>.
