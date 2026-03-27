# PCDWBA

## Overview

Validated Experimental

*Model-family pages:*
[Overview](https://brandynlucca.github.io/acousticTS/articles/pcdwba/index.md)
[Implementation](https://brandynlucca.github.io/acousticTS/articles/pcdwba/pcdwba-implementation.md)
[Theory](https://brandynlucca.github.io/acousticTS/articles/pcdwba/pcdwba-theory.md)

The phase-compensated distorted wave Born approximation (`PCDWBA`) is
the curved-body extension of the weak-scattering `DWBA`. It keeps the
same local fluid-like cylindrical kernel but corrects the along-body
phase accumulation for a bent centerline.

### Core idea

Write the weak-scattering backscatter as a sum over local cross-sections
along a bent body, and retain the position- and tilt-dependent phase
that would be lost by treating the target as straight.

### Best for

- Weakly scattering elongated targets with meaningful curvature
- Bent zooplankton-like bodies parameterized by a centerline and local
  radius profile
- Curved-body problems where a straight `DWBA` is too restrictive

### Supports

- Canonical bent cylinders and arbitrary fluid-like centerline profiles
- Weak-contrast notation relative to seawater as g\_{21} and h\_{21}
- Monostatic and phase-sensitive curved-body backscatter calculations

### Main assumptions

- Born-type weak-scattering regime
- Curvature modifies phase bookkeeping but not the local weak-fluid
  kernel
- Single scattering from an elongated fluid-like target

### Validation status

- Validated against source-level `ZooScatR` and `echopop` PCDWBA
  workflows.

### Family pages

- [Implementation](https://brandynlucca.github.io/acousticTS/articles/pcdwba/pcdwba-implementation.md):
  object workflow and comparison results
- [Theory](https://brandynlucca.github.io/acousticTS/articles/pcdwba/pcdwba-theory.md):
  bent-centerline parameterization and phase-compensated weak-scattering
  sum
