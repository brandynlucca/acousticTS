# SOEMS

## Overview

Benchmarked Validated

*Model-family pages:*
[Overview](https://brandynlucca.github.io/acousticTS/articles/calibration/index.md)
[Implementation](https://brandynlucca.github.io/acousticTS/articles/calibration/calibration-implementation.md)
[Theory](https://brandynlucca.github.io/acousticTS/articles/calibration/calibration-theory.md)

`SOEMS` is the package’s elastic solid-sphere family for standard-target
calibration work. It follows the classic solid elastic sphere literature
and is intended for metal calibration spheres whose compressional and
shear wave physics matter directly.

### Core idea

Expand the exterior acoustic field and the interior elastic displacement
field in spherical modes, enforce pressure, normal-velocity, and
traction conditions at the sphere surface, and reconstruct the far-field
backscatter from the resulting phase-shifted partial waves.

### Best for

- Calibration spheres made of tungsten carbide, copper, aluminum, steel,
  or similar elastic solids
- Frequency-response calculations for standard-target calibration
  workflows
- Benchmarking elastic spherical scattering against literature values

### Supports

- `CAL` objects representing homogeneous elastic spheres
- Solid elastic interiors with both longitudinal and shear wave support
- Medium `1` as seawater and medium `2` as the elastic solid sphere

### Main assumptions

- Perfectly spherical homogeneous solid target
- Linear elasticity in the solid interior
- Inviscid exterior fluid
- No shell or internal cavity

### Validation status

- Benchmarked against published calibration-sphere targets used
  throughout the package documentation.
- Validated against `echoSMs`, `sphereTS`, and the NOAA calibration
  applet.

### Family pages

- [Implementation](https://brandynlucca.github.io/acousticTS/articles/calibration/calibration-implementation.md):
  calibration workflows and benchmark comparisons
- [Theory](https://brandynlucca.github.io/acousticTS/articles/calibration/calibration-theory.md):
  elastic-sphere modal theory and phase-shift representation
