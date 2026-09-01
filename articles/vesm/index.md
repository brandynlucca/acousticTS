# VESM

## Overview

Validated Experimental

[Theory](https://brandynlucca.github.io/acousticTS/articles/vesm/vesm-theory.md)
[Implementation](https://brandynlucca.github.io/acousticTS/articles/vesm/vesm-implementation.md)

The viscous-elastic spherical model (`VESM`) calculates resonant
scattering from a gas core, elastic shell, and outer viscous layer. The
package interface uses `model = "vesms"` and stores output under
`VESMS`.

### Core idea

Solve the concentric layered sphere mode by mode, including acoustic
waves in the exterior and gas, elastic waves in the shell, and complex
compressional and shear response in the viscous layer ([Khodabandeloo et
al. 2021](#ref-Khodabandeloo_2021); [Feuillade and Nero
1998](#ref-Feuillade_1998)).

### Best for

- Gas-bearing layered targets represented by a concentric spherical
  idealization
- Resonance studies in which shell elasticity and viscous damping both
  matter
- Comparison with simpler bubble, fluid-sphere, or elastic-shell models

### Supports

- Spherical `ESS` objects containing the elastic shell and inner
  gas-like fluid
- Run-time viscous-layer density, sound speed, shear viscosity, bulk
  viscosity, and thickness or radius
- Configurable spherical modal truncation

### Main assumptions

- Concentric spherical interfaces and homogeneous properties in each
  layer
- Linear acoustics, elasticity, and the documented viscous constitutive
  model
- No nonspherical posture, eccentric core, or spatial material gradients
- Validation is limited to the documented layered-sphere cases

### Validation status

- Validated against the reference Python VESM implementation on the
  documented layered-sphere case.
- VESM is currently marked experimental because the documented public
  workflow is still limited to the current layered-sphere scope.

### Family pages

- [Implementation](https://brandynlucca.github.io/acousticTS/articles/vesm/vesm-implementation.md):
  building the layered object, required arguments, output, and
  comparison
- [Theory](https://brandynlucca.github.io/acousticTS/articles/vesm/vesm-theory.md):
  fields in all four regions and interface systems

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
