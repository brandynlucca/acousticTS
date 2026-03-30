# HPA

## Overview

Benchmarked Validated

[Theory](https://brandynlucca.github.io/acousticTS/articles/hpa/hpa-theory.md)
[Implementation](https://brandynlucca.github.io/acousticTS/articles/hpa/hpa-implementation.md)

These pages follow Johnson’s asymptotic fluid-sphere formulation and
Stanton’s generalized approximate backscatter formulas ([Johnson
1977](#ref-johnson_sound_1977); [Stanton 1989](#ref-stanton_1989)).

The high-pass approximation (`HPA`) is a compact asymptotic backscatter
family that interpolates between a Rayleigh-style low-frequency limit
and a reflection-controlled high-frequency limit.

### Core idea

Build a rational approximation whose numerator reproduces the low-`ka`
scattering strength and whose denominator suppresses unphysical growth
as frequency increases, then adapt the geometric prefactors to spheres,
spheroids, and cylinders.

### Best for

- Fast approximate spectra for simple canonical bodies
- Broad trend studies when an exact modal solve is unnecessary
- Weakly contrasting or moderately reflecting targets represented by
  simple shapes

### Supports

- Sphere, prolate spheroid, straight cylinder, and bent-cylinder
  branches
- Contrast bookkeeping relative to seawater as medium `1`
- Very fast monostatic backscatter estimates

### Main assumptions

- Asymptotic interpolation rather than an exact boundary-value solution
- Shape-specific prefactors carried from the source literature
- Best interpreted as a broadband approximation, not a
  resonance-resolving solver

### Validation status

- Benchmarked against canonical asymptotic target families rather than
  as an exact modal solver.
- Validated against the spherical `echoSMs::HPModel` branch and the
  published Johnson/Stanton algebra.

### Family pages

- [Implementation](https://brandynlucca.github.io/acousticTS/articles/hpa/hpa-implementation.md):
  quick workflows and validation tables
- [Theory](https://brandynlucca.github.io/acousticTS/articles/hpa/hpa-theory.md):
  Rayleigh term, reflection limit, and shape-specific completions

## References

Johnson, Richard K. 1977. “Sound Scattering from a Fluid Sphere
Revisited.” *The Journal of the Acoustical Society of America* 61 (2):
375–77. <https://doi.org/10.1121/1.381326>.

Stanton, Timothy K. 1989. “Simple Approximate Formulas for
Backscattering of Sound by Spherical and Elongated Objects.” *The
Journal of the Acoustical Society of America* 86 (4): 1499–1510.
<https://doi.org/10.1121/1.398711>.
