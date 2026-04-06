# SOEMS

## Overview

Benchmarked Validated

[Theory](https://brandynlucca.github.io/acousticTS/articles/calibration/calibration-theory.md)
[Implementation](https://brandynlucca.github.io/acousticTS/articles/calibration/calibration-implementation.md)

These pages are grounded in the standard-target calibration literature
for elastic reference spheres ([Dragonette, Numrich, and Frank
1981](#ref-Dragonette_1981); [Foote 1990](#ref-Foote_1990); [MacLennan
1981](#ref-Maclennan_1981)).

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

## References

Dragonette, Louis R., S. K. Numrich, and Laurence J. Frank. 1981.
“Calibration Technique for Acoustic Scattering Measurements.” *The
Journal of the Acoustical Society of America* 69 (4): 1186–89.
<https://doi.org/10.1121/1.385699>.

Foote, K. G. 1990. “Spheres for Calibrating an Eleven-Frequency Acoustic
Measurement System.” *ICES Journal of Marine Science* 46 (3): 284–86.
<https://doi.org/10.1093/icesjms/46.3.284>.

MacLennan, D. N. 1981. “The Theory of Solid Spheres as Sonar Calibration
Targets.” Scottish Fisheries Research Report 22. Department of
Agriculture; Fisheries for Scotland.
