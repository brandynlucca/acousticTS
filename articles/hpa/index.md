# HPA

## Overview

Benchmarked Validated

*Model-family pages:*
[Overview](https://brandynlucca.github.io/acousticTS/articles/hpa/index.md)
[Implementation](https://brandynlucca.github.io/acousticTS/articles/hpa/hpa-implementation.md)
[Theory](https://brandynlucca.github.io/acousticTS/articles/hpa/hpa-theory.md)

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
