# ESSMS

## Overview

Unvalidated

[Theory](https://brandynlucca.github.io/acousticTS/articles/essms/essms-theory.md)
[Implementation](https://brandynlucca.github.io/acousticTS/articles/essms/essms-implementation.md)

The elastic-shelled spherical modal-series model (`ESSMS`) calculates
scattering from an elastic spherical shell around an inner fluid or
cavity.

### Core idea

Represent the exterior and interior acoustic fields and the shell’s
longitudinal and transverse elastic fields in spherical modes. Match
pressure, displacement, and traction at both interfaces ([Goodman and
Stern 1962](#ref-Goodman_1962); [Faran 1951](#ref-Faran_1951)).

### Best for

- Spherical elastic shells with fluid or gas-like interiors
- Shell-resonance and thickness-sensitivity studies
- Canonical layered-sphere calculations where shell elasticity must be
  retained

### Supports

- Spherical `ESS` objects with an explicit shell and inner fluid
- Shell density, elastic moduli, longitudinal speed, and shear speed
- User-controlled modal truncation

### Main assumptions

- Concentric spherical interfaces
- Homogeneous isotropic shell and homogeneous inner fluid
- Linear elasticity and linear acoustics
- No viscosity, nonspherical deformation, or eccentric core

### Validation status

- The package does not yet claim external validation across the current
  public scope.

### Family pages

- [Implementation](https://brandynlucca.github.io/acousticTS/articles/essms/essms-implementation.md):
  object construction, output, and current validation scope
- [Theory](https://brandynlucca.github.io/acousticTS/articles/essms/essms-theory.md):
  layered fields, interface tractions, and modal systems

## References

Faran, James J. 1951. “Sound Scattering by Solid Cylinders and Spheres.”
*The Journal of the Acoustical Society of America* 23 (4): 405–18.
<https://doi.org/10.1121/1.1906780>.

Goodman, Ralph R., and Raya Stern. 1962. “Reflection and Transmission of
Sound by Elastic Spherical Shells.” *The Journal of the Acoustical
Society of America* 34 (3): 338–44. <https://doi.org/10.1121/1.1928120>.
