# PCDWBA

## Overview

Validated Experimental

*Model-family pages:*
[Overview](https://brandynlucca.github.io/acousticTS/articles/pcdwba/index.md)
[Implementation](https://brandynlucca.github.io/acousticTS/articles/pcdwba/pcdwba-implementation.md)
[Theory](https://brandynlucca.github.io/acousticTS/articles/pcdwba/pcdwba-theory.md)

These pages follow the phase-compensated weak-scattering literature for
broadside elongated bodies and krill-style applications ([Chu and Ye
1999](#ref-chu_phase-compensated_1999); [Chu, Foote, and Stanton
1993](#ref-chu_further_1993); [Stanton 1989](#ref-stanton_sound_1989)).

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

## References

Chu, Dezhang, Kenneth G. Foote, and Timothy K. Stanton. 1993. “Further
Analysis of Target Strength Measurements of Antarctic Krill at 38 and
120 kHz: Comparison with Deformed Cylinder Model and Inference of
Orientation Distribution.” *The Journal of the Acoustical Society of
America* 93 (5): 2985–88. <https://doi.org/10.1121/1.405818>.

Chu, Dezhang, and Zhen Ye. 1999. “A Phase-Compensated Distorted Wave
Born Approximation Representation of the Bistatic Scattering by Weakly
Scattering Objects: Application to Zooplankton.” *The Journal of the
Acoustical Society of America* 106 (4): 1732–43.
<https://doi.org/10.1121/1.428036>.

Stanton, T. K. 1989. “Sound Scattering by Cylinders of Finite Length.
III. Deformed Cylinders.” *The Journal of the Acoustical Society of
America* 86 (2): 691–705. <https://doi.org/10.1121/1.398193>.
