# ESSMS

## Overview

Unvalidated

*Model-family pages:*
[Overview](https://brandynlucca.github.io/acousticTS/articles/essms/index.md)
[Implementation](https://brandynlucca.github.io/acousticTS/articles/essms/essms-implementation.md)
[Theory](https://brandynlucca.github.io/acousticTS/articles/essms/essms-theory.md)

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

- A direct shell-sphere benchmark family exists, but the current ESSMS
  implementation still does not return finite full-grid benchmark
  spectra.

### Family pages

- [Implementation](https://brandynlucca.github.io/acousticTS/articles/essms/essms-implementation.md):
  current usage and validation status
- [Theory](https://brandynlucca.github.io/acousticTS/articles/essms/essms-theory.md):
  full layered shell derivation, tractions, and mode-wise boundary
  systems
