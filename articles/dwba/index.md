# DWBA

## Overview

Benchmarked Validated

[Theory](https://brandynlucca.github.io/acousticTS/articles/dwba/dwba-theory.md)
[Implementation](https://brandynlucca.github.io/acousticTS/articles/dwba/dwba-implementation.md)

The distorted-wave Born approximation (`DWBA`) calculates monostatic
backscatter from weakly contrasting fluid-like bodies represented by a
segmented centerline and local radius profile.

### Core idea

Linearize the scattering response in density and compressibility
contrast, calculate each local cross-section, and integrate the complex
contributions with their two-way phase ([Chu et al.
1993](#ref-Chu_1993); [Stanton et al. 1998](#ref-Stanton_1998_2)).

### Best for

- Weakly scattering zooplankton and other fluid-like elongated bodies
- Arbitrary profiles that are not well represented by one sphere,
  cylinder, or spheroid
- Deterministic calculations for a specified geometry and orientation

### Supports

- `FLS` scatterers with canonical or arbitrary body profiles
- Density and sound-speed contrasts relative to the surrounding water
- Frequency-dependent complex amplitude, cross-section, and target
  strength

### Main assumptions

- Small density and sound-speed contrasts
- First-order Born-type scattering with no strong internal reverberation
- A fluid-like interior without elastic shear response
- Sufficient body discretization for the local radius and phase
  variation

### Validation status

- Benchmarked against the canonical spectra stored in benchmark_ts.
- Validated against the published McGehee et al (1998) and echoSMs
  workflows.

### Family pages

- [Implementation](https://brandynlucca.github.io/acousticTS/articles/dwba/dwba-implementation.md):
  object preparation, discretization, output, and comparisons
- [Theory](https://brandynlucca.github.io/acousticTS/articles/dwba/dwba-theory.md):
  weak-scattering contrast and centerline integral

## References

Chu, Dezhang, Kenneth G. Foote, and Timothy K. Stanton. 1993. “Further
Analysis of Target Strength Measurements of Antarctic Krill at 38 and
120 kHz: Comparison with Deformed Cylinder Model and Inference of
Orientation Distribution.” *The Journal of the Acoustical Society of
America* 93 (5): 2985–88. <https://doi.org/10.1121/1.405818>.

Stanton, Timothy K., Dezhang Chu, and Peter H. Wiebe. 1998. “Sound
Scattering by Several Zooplankton Groups. II. Scattering Models.” *The
Journal of the Acoustical Society of America* 103 (1): 236–53.
<https://doi.org/10.1121/1.421110>.
