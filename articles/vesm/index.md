# VESM

## Overview

Validated Experimental

*Model-family pages:*
[Overview](https://brandynlucca.github.io/acousticTS/articles/vesm/index.md)
[Implementation](https://brandynlucca.github.io/acousticTS/articles/vesm/vesm-implementation.md)
[Theory](https://brandynlucca.github.io/acousticTS/articles/vesm/vesm-theory.md)

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

### Family pages

- [Implementation](https://brandynlucca.github.io/acousticTS/articles/vesm/vesm-implementation.md):
  layered-object workflow and benchmark comparisons
- [Theory](https://brandynlucca.github.io/acousticTS/articles/vesm/vesm-theory.md):
  full layered spherical field structure, interface conditions, and
  mode-wise systems
