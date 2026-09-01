# ECMS

## Overview

Unvalidated Experimental

[Theory](https://brandynlucca.github.io/acousticTS/articles/ecms/ecms-theory.md)
[Implementation](https://brandynlucca.github.io/acousticTS/articles/ecms/ecms-implementation.md)

The elastic-cylinder modal-series model (`ECMS`) calculates
near-broadside scattering from a finite solid elastic cylinder. It
retains both longitudinal and shear waves in the solid.

### Core idea

Solve the infinite circular elastic-cylinder interface through modal
phase shifts, then apply a finite-length coherence factor. Uniform
curvature can be represented by the corresponding Fresnel correction
([Faran 1951](#ref-Faran_1951); [Stanton 1988](#ref-Stanton_1988)).

### Best for

- Straight or uniformly bent solid elastic cylinders near broadside
- Backbone-like canonical components where shear-wave support matters
- Elastic-cylinder reference and sensitivity calculations

### Supports

- Cylindrical `ESS` objects, with the elastic component stored in the
  shell slot
- Legacy cylindrical `FLS` geometry carriers for backward compatibility
- User-supplied density and longitudinal and transverse elastic wave
  speeds

### Main assumptions

- A homogeneous isotropic solid with circular cross-section
- Broadside or near-broadside incidence
- Finite length and uniform curvature represented through coherence
  factors
- No surrounding flesh-body coupling or arbitrary elastic geometry

### Validation status

- ECMS is currently marked experimental because the documented checks
  are independent algebra reconstructions rather than an external
  benchmark or software-comparison ladder.

### Family pages

- [Implementation](https://brandynlucca.github.io/acousticTS/articles/ecms/ecms-implementation.md):
  supported objects, material inputs, output, and current checks
- [Theory](https://brandynlucca.github.io/acousticTS/articles/ecms/ecms-theory.md):
  elastic potentials, phase shifts, and finite-length closure

## References

Faran, James J. 1951. “Sound Scattering by Solid Cylinders and Spheres.”
*The Journal of the Acoustical Society of America* 23 (4): 405–18.
<https://doi.org/10.1121/1.1906780>.

Stanton, T. K. 1988. “Sound Scattering by Cylinders of Finite Length. I.
Fluid Cylinders.” *The Journal of the Acoustical Society of America* 83
(1): 55–63. <https://doi.org/10.1121/1.396184>.
