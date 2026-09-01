# PCDWBA

## Overview

Validated Experimental

[Theory](https://brandynlucca.github.io/acousticTS/articles/pcdwba/pcdwba-theory.md)
[Implementation](https://brandynlucca.github.io/acousticTS/articles/pcdwba/pcdwba-implementation.md)

The phase-compensated distorted-wave Born approximation (`PCDWBA`)
calculates weak-fluid scattering from an elongated body while retaining
the phase associated with its curved centerline.

### Core idea

Use the local cross-sectional kernel of `DWBA`, but calculate position,
local tilt, spacing, and two-way phase from the actual centerline rather
than from a straight-axis reduction ([Chu and Ye 1999](#ref-Chu_1999)).

### Best for

- Weakly contrasting elongated bodies with physically meaningful
  curvature
- Uniformly bent cylinders or curved arbitrary profiles
- Assessing the phase effect of curvature relative to a straight
  weak-scattering calculation

### Supports

- `FLS` scatterers converted to a DWBA-style radius and centerline
  profile
- Curvature already present in the profile or supplied as a radius or
  radius-to-length ratio
- Nodewise density and sound-speed contrast bookkeeping

### Main assumptions

- The same first-order weak-scattering regime as `DWBA`
- Curvature changes position and phase but not the local scattering
  kernel
- Single scattering from a fluid-like target
- The centerline and local radius discretization resolve the target
  geometry

### Validation status

- Validated against source-level ZooScatR and Echopop PCDWBA workflows.
- PCDWBA is currently marked experimental because the public package
  workflow is still being tightened even though the current source-
  level comparison cases are documented.

### Family pages

- [Implementation](https://brandynlucca.github.io/acousticTS/articles/pcdwba/pcdwba-implementation.md):
  accepted profiles, curvature inputs, output, and comparisons
- [Theory](https://brandynlucca.github.io/acousticTS/articles/pcdwba/pcdwba-theory.md):
  centerline geometry and phase-compensated sum

## References

Chu, Dezhang, and Zhen Ye. 1999. “A Phase-Compensated Distorted Wave
Born Approximation Representation of the Bistatic Scattering by Weakly
Scattering Objects: Application to Zooplankton.” *The Journal of the
Acoustical Society of America* 106 (4): 1732–43.
<https://doi.org/10.1121/1.428036>.
