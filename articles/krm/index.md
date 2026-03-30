# KRM

## Overview

Benchmarked Validated

[Theory](https://brandynlucca.github.io/acousticTS/articles/krm/krm-theory.md)
[Implementation](https://brandynlucca.github.io/acousticTS/articles/krm/krm-implementation.md)

These pages follow the composite body-plus-swimbladder fish modeling
literature initiated for cod and later generalized in open software
implementations ([C. S. Clay 1991](#ref-clay_1991); [Clarence S. Clay
and Horne 1994](#ref-clay_horne_1994); [Gastauer
2025](#ref-sven_gastauer_svengastauerkrmr_2025)).

The Kirchhoff-ray mode model (`KRM`) is the package’s composite fish
family for targets whose body and swimbladder occupy different acoustic
regimes. It keeps a weakly contrasting ray-style body treatment and a
separate swimbladder treatment that switches between low-mode and
high-frequency behavior.

### Core idea

Treat the fish body with a Kirchhoff-style short-segment approximation,
treat the swimbladder with a simplified cylinder-based modal or
high-frequency branch depending on acoustic size, and combine the
complex component amplitudes coherently.

### Best for

- Gas-bearing fish whose body and swimbladder should not be modeled by
  the same exact family
- Fish-like profile data represented by segmented body and bladder
  outlines
- Practical fisheries workflows where composite body-plus-swimbladder
  structure matters

### Supports

- `SBF` composite scatterers
- Body contrasts relative to seawater as medium `1` and swimbladder
  contrasts relative to the body in the local bladder subproblem
- Monostatic composite target strength

### Main assumptions

- Short-segment ray-style treatment of the body
- Simplified swimbladder physics chosen by acoustic size regime
- Coherent component combination without a full coupled body-bladder
  boundary-value solve

### Validation status

- Benchmarked against canonical modal-family targets used for isolated
  gas-filled and weakly scattering cases.
- Validated against `KRMr`, `echoSMs`, and the NOAA KRM applet on
  bundled fish objects and shared workflows.

### Family pages

- [Implementation](https://brandynlucca.github.io/acousticTS/articles/krm/krm-implementation.md):
  scatterer setup, spectra, and validation tables
- [Theory](https://brandynlucca.github.io/acousticTS/articles/krm/krm-theory.md):
  body Kirchhoff reduction, swimbladder mode/ray branches, and coherent
  composite sum

## References

Clay, C. S. 1991. “Low-Resolution Acoustic Scattering Models:
Fluid-Filled Cylinders and Fish with Swim Bladders.” *The Journal of the
Acoustical Society of America* 89 (5): 2168–79.
<https://doi.org/10.1121/1.400910>.

Clay, Clarence S., and John K. Horne. 1994. “Acoustic Models of Fish:
The Atlantic Cod (*Gadus Morhua*).” *The Journal of the Acoustical
Society of America* 96 (3): 1661–68. <https://doi.org/10.1121/1.410245>.

Gastauer, Sven. 2025. “SvenGastauer/KRMr: V0.4.8.” Zenodo.
<https://doi.org/10.5281/ZENODO.15838374>.
