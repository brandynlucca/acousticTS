# HPA

## Overview

Benchmarked Validated

[Theory](https://brandynlucca.github.io/acousticTS/articles/hpa/hpa-theory.md)
[Implementation](https://brandynlucca.github.io/acousticTS/articles/hpa/hpa-implementation.md)

The high-pass approximation (`HPA`) provides fast broadband estimates
for simple fluid-like canonical bodies by joining low-frequency and
reflection-controlled high-frequency behaviour.

### Core idea

Use a rational form that recovers the Rayleigh-scale response at small
acoustic size while limiting its high-frequency growth. The Johnson
method is sphere-specific, while the Stanton method also covers prolate
spheroids and cylinders ([Johnson 1977](#ref-Johnson_1977); [Stanton
1989](#ref-Stanton_1989_1)).

### Best for

- Rapid broadband trend calculations for spheres, prolate spheroids, and
  cylinders
- Screening shape, size, orientation, or material-contrast effects
- Comparison with a geometry-matched modal solution when fine resonances
  are not required

### Supports

- The Johnson sphere and Stanton sphere, prolate-spheroid, and cylinder
  formulations
- Fluid-like density and sound-speed contrasts relative to the exterior
  medium
- Optional null-position and deviation controls

### Main assumptions

- An interpolation formula rather than an exact boundary-value solution
- A supported canonical shape and the source formulation’s geometric
  prefactors
- No complete modal resonance structure
- Results are interpreted over the approximation’s acoustic-size regime

### Validation status

- Benchmarked against the canonical spherical spectra stored in
  benchmark_ts.
- Validated against the spherical echoSMs implementation.

### Family pages

- [Implementation](https://brandynlucca.github.io/acousticTS/articles/hpa/hpa-implementation.md):
  method selection, supported shapes, output, and comparisons
- [Theory](https://brandynlucca.github.io/acousticTS/articles/hpa/hpa-theory.md):
  low-frequency term, reflection limit, and shape-specific forms

## References

Johnson, Richard K. 1977. “Sound Scattering from a Fluid Sphere
Revisited.” *The Journal of the Acoustical Society of America* 61 (2):
375–77. <https://doi.org/10.1121/1.381326>.

Stanton, Timothy K. 1989. “Simple Approximate Formulas for
Backscattering of Sound by Spherical and Elongated Objects.” *The
Journal of the Acoustical Society of America* 86 (4): 1499–510.
<https://doi.org/10.1121/1.398711>.
