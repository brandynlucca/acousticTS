# SPHMS

## Overview

Benchmarked Validated

[Theory](https://brandynlucca.github.io/acousticTS/articles/sphms/sphms-theory.md)
[Implementation](https://brandynlucca.github.io/acousticTS/articles/sphms/sphms-implementation.md)

The spherical modal-series solution (`SPHMS`) calculates scattering from
an unshelled, homogeneous sphere in an exterior fluid. It is the
package’s primary spherical reference model.

### Core idea

Expand the incident, scattered, and interior fields in spherical partial
waves, apply the selected interface condition independently at each
angular order, and sum the far-field coefficients. This separation is
exact for a spherical interface ([Anderson 1950](#ref-Anderson_1950);
[Hickling 1962](#ref-Hickling_1962)).

### Best for

- Fixed-rigid, pressure-release, liquid-filled, or gas-filled spheres
- Canonical checks of boundary and material-property effects
- Spherical reference cases for approximate or more general solvers

### Supports

- `Sphere` geometry carried by `FLS` or `GAS` scatterers
- Homogeneous interiors with seawater as exterior medium `1`
- Monostatic complex amplitude, backscattering cross-section, and target
  strength

### Main assumptions

- A perfectly spherical, unshelled interface
- Homogeneous properties within each region
- Linear time-harmonic acoustics
- No elastic shear waves, viscous layer, or nonspherical orientation
  effects

### Validation status

- Benchmarked against the canonical spherical spectra stored in
  benchmark_ts.
- Validated against `KRMr` and `echoSMs` on shared penetrable-sphere
  cases.

### Family pages

- [Implementation](https://brandynlucca.github.io/acousticTS/articles/sphms/sphms-implementation.md):
  object construction, boundaries, output, and comparisons
- [Theory](https://brandynlucca.github.io/acousticTS/articles/sphms/sphms-theory.md):
  field expansions, boundary systems, and modal coefficients

## References

Anderson, Victor C. 1950. “Sound Scattering from a Fluid Sphere.” *The
Journal of the Acoustical Society of America* 22 (4): 426–31.
<https://doi.org/10.1121/1.1906621>.

Hickling, Robert. 1962. “Analysis of Echoes from a Solid Elastic Sphere
in Water.” *The Journal of the Acoustical Society of America* 34 (10):
1582–92. <https://doi.org/10.1121/1.1909055>.
