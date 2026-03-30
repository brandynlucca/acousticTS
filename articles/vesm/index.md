# VESM

## Overview

Validated Experimental

[Theory](https://brandynlucca.github.io/acousticTS/articles/vesm/vesm-theory.md)
[Implementation](https://brandynlucca.github.io/acousticTS/articles/vesm/vesm-implementation.md)

These pages are motivated by layered gas-bearing fish scattering models
and viscous resonance broadening ([Khodabandeloo et al.
2021](#ref-khodabandeloo_estimating_2021); [Feuillade and Nero
1998](#ref-feuillade_nero_1998)).

The viscous-elastic spherical model (`VESM`) is a layered spherical
resonance family for a gas core surrounded by an elastic shell and an
outer viscous biological layer.

### Core idea

Solve the spherical layered problem mode by mode, with acoustic waves in
the exterior and gas core, viscous compressional/shear branches in the
soft outer layer, and elastic longitudinal/transverse waves in the
shell.

### Best for

- Gas-bearing mesopelagic-fish style layered spheres
- Problems where gas resonance, shell elasticity, and viscous damping
  all matter simultaneously
- Wideband spherical resonance studies beyond a simple bubble or
  shell-only idealization

### Supports

- `ESS` spherical objects with an explicit shell and inner gas core plus
  viscous-layer model arguments
- Four-region spherical bookkeeping with seawater as medium `1`
- Monostatic layered resonance calculations

### Main assumptions

- Spherical symmetry
- Homogeneous material properties within each layer
- Linear acoustics and linear elasticity
- Layered spherical interfaces with no nonspherical posture or geometry
  effects

### Validation status

- Validated against the reference Python VESM implementation on the
  documented layered-sphere case.
- VESM is currently marked experimental because the documented public
  workflow is still limited to the current layered-sphere scope.

### Family pages

- [Implementation](https://brandynlucca.github.io/acousticTS/articles/vesm/vesm-implementation.md):
  layered-object workflow and benchmark comparisons
- [Theory](https://brandynlucca.github.io/acousticTS/articles/vesm/vesm-theory.md):
  full layered spherical field structure, interface conditions, and
  mode-wise systems

## References

Feuillade, C., and R. W. Nero. 1998. “A Viscous-Elastic Swimbladder
Model for Describing Enhanced-Frequency Resonance Scattering from Fish.”
*The Journal of the Acoustical Society of America* 103 (6): 3245–55.
<https://doi.org/10.1121/1.423076>.

Khodabandeloo, Babak, Mette Dalgaard Agersted, Thor Klevjer, Gavin J.
Macaulay, and Webjørn Melle. 2021. “Estimating Target Strength and
Physical Characteristics of Gas-Bearing Mesopelagic Fish from Wideband
*in Situ* Echoes Using a Viscous-Elastic Scattering Model.” *The Journal
of the Acoustical Society of America* 149 (1): 673–91.
<https://doi.org/10.1121/10.0003341>.
