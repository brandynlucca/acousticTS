# KRM

## Overview

Benchmarked Validated

*Model-family pages:*
[Overview](https://brandynlucca.github.io/acousticTS/articles/krm/index.md)
[Implementation](https://brandynlucca.github.io/acousticTS/articles/krm/krm-implementation.md)
[Theory](https://brandynlucca.github.io/acousticTS/articles/krm/krm-theory.md)

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
