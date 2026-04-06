# SPHMS

## Overview

Benchmarked Validated

[Theory](https://brandynlucca.github.io/acousticTS/articles/sphms/sphms-theory.md)
[Implementation](https://brandynlucca.github.io/acousticTS/articles/sphms/sphms-implementation.md)

These pages follow the classical exact sphere literature for fluid,
elastic, and approximate shell limits ([Anderson
1950](#ref-Anderson_1950); [Faran 1951](#ref-Faran_1951); [Hickling
1962](#ref-Hickling_1962)).

The spherical modal series solution (`SPHMS`) is the package’s exact
canonical solution for unshelled homogeneous spheres in an exterior
fluid. Because the Helmholtz equation separates exactly in spherical
coordinates, each angular order remains algebraically local after the
boundary conditions are imposed.

### Core idea

Represent the incident, scattered, and interior fields in spherical
partial waves, enforce the chosen boundary condition at the spherical
interface, and sum the retained orders into the far-field backscatter.
`SPHMS` is the natural spherical benchmark for the rest of the package.

### Best for

- Rigid, pressure-release, liquid-filled, and gas-filled spheres
- Canonical benchmark comparisons and sanity checks
- Problems where exact spherical geometry is a reasonable physical
  idealization

### Supports

- `Sphere` shapes on `FLS` or `GAS` objects
- Medium-indexed spherical theory with medium `1` as seawater and medium
  `2` as the sphere interior
- Direct target-strength calculations from exact spherical coefficients

### Main assumptions

- Perfectly spherical interface
- Homogeneous material properties in each region
- Linear, time-harmonic acoustics
- No elastic shell or viscous intermediate layer

### Validation status

- Benchmarked against the canonical spherical spectra stored in
  `benchmark_ts`.
- Validated against `KRMr` and `echoSMs` on shared penetrable-sphere
  cases.

### Family pages

- [Implementation](https://brandynlucca.github.io/acousticTS/articles/sphms/sphms-implementation.md):
  object construction, spectra, and software comparisons
- [Theory](https://brandynlucca.github.io/acousticTS/articles/sphms/sphms-theory.md):
  full spherical separation, boundary conditions, and coefficient
  derivations

## References

Anderson, Victor C. 1950. “Sound Scattering from a Fluid Sphere.” *The
Journal of the Acoustical Society of America* 22 (4): 426–31.
<https://doi.org/10.1121/1.1906621>.

Faran, James J. 1951. “Sound Scattering by Solid Cylinders and Spheres.”
*The Journal of the Acoustical Society of America* 23 (4): 405–18.
<https://doi.org/10.1121/1.1906780>.

Hickling, Robert. 1962. “Analysis of Echoes from a Solid Elastic Sphere
in Water.” *The Journal of the Acoustical Society of America* 34 (10):
1582–92. <https://doi.org/10.1121/1.1909055>.
