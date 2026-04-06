# ESSMS

## Overview

Unvalidated

[Theory](https://brandynlucca.github.io/acousticTS/articles/essms/essms-theory.md)
[Implementation](https://brandynlucca.github.io/acousticTS/articles/essms/essms-implementation.md)

These pages are grounded in the classical elastic-shell scattering
literature for fluid-filled spherical shells ([Goodman and Stern
1962](#ref-Goodman_1962); [Faran 1951](#ref-Faran_1951); [Stanton
1990](#ref-Stanton_1990)).

The elastic shelled spherical modal series (`ESSMS`) is the package’s
layered spherical family for an elastic shell surrounding either a fluid
interior or a pressure-release cavity.

### Core idea

Represent the exterior acoustic field, the shell’s longitudinal and
transverse elastic fields, and the interior acoustic field in spherical
modes, then match pressure, displacement, and traction conditions at
both shell interfaces.

### Best for

- Spherical elastic shells with fluid or gas interiors
- Layered shell problems where shell resonances matter explicitly
- Theoretical reference work on spherical shell scattering

### Supports

- `ESS` spherical geometries
- Exterior seawater (medium `1`), shell (medium `2`), and interior fluid
  or cavity (medium `3`)
- Pressure-release and fluid-filled shell interiors

### Main assumptions

- Perfectly spherical shell geometry
- Homogeneous isotropic shell material
- Linear elasticity and linear acoustics
- Modal truncation of an exact spherical layered system

### Validation status

- The package does not yet claim external validation across the current
  public scope.

### Family pages

- [Implementation](https://brandynlucca.github.io/acousticTS/articles/essms/essms-implementation.md):
  current usage and validation status
- [Theory](https://brandynlucca.github.io/acousticTS/articles/essms/essms-theory.md):
  full layered shell derivation, tractions, and mode-wise boundary
  systems

## References

Faran, James J. 1951. “Sound Scattering by Solid Cylinders and Spheres.”
*The Journal of the Acoustical Society of America* 23 (4): 405–18.
<https://doi.org/10.1121/1.1906780>.

Goodman, Ralph R., and Raya Stern. 1962. “Reflection and Transmission of
Sound by Elastic Spherical Shells.” *The Journal of the Acoustical
Society of America* 34 (3): 338–44. <https://doi.org/10.1121/1.1928120>.

Stanton, T. K. 1990. “Sound Scattering by Spherical and Elongated
Shelled Bodies.” *The Journal of the Acoustical Society of America* 88
(3): 1619–33. <https://doi.org/10.1121/1.400321>.
