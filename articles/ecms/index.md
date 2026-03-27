# ECMS

## Overview

Experimental Unvalidated

*Model-family pages:*
[Overview](https://brandynlucca.github.io/acousticTS/articles/ecms/index.md)
[Implementation](https://brandynlucca.github.io/acousticTS/articles/ecms/ecms-implementation.md)
[Theory](https://brandynlucca.github.io/acousticTS/articles/ecms/ecms-theory.md)

The elastic cylinder modal series solution (`ECMS`) is the package’s
solid elastic-cylinder family. It combines the phase-shift treatment of
an infinite elastic circular cylinder with the usual finite-length
coherence factor used near broadside.

### Core idea

Represent the elastic interior by longitudinal and transverse
cylindrical waves, determine the phase shifts induced by the boundary
conditions, and apply the near-broadside finite-length coherence factor
to recover monostatic backscatter.

### Best for

- Solid elastic finite cylinders near broadside
- Backbone-like cylindrical structures where shear-wave support matters
- Canonical elastic-cylinder benchmarks and component models

### Supports

- `ESS` or compatible elastic-cylinder scatterers
- Longitudinal and transverse elastic wave speeds in the cylinder
  interior
- Straight and current bent broadside branches

### Main assumptions

- Circular cylinder with homogeneous elastic properties
- Broadside or near-broadside incidence
- Finite-length treatment through a coherence factor rather than a fully
  general 3D elastic solve

### Validation status

- Current ECMS checks are independent algebra reconstructions rather
  than a documented external benchmark ladder.

### Family pages

- [Implementation](https://brandynlucca.github.io/acousticTS/articles/ecms/ecms-implementation.md):
  object setup and current comparison workflow
- [Theory](https://brandynlucca.github.io/acousticTS/articles/ecms/ecms-theory.md):
  elastic cylindrical potentials, phase shifts, and finite-length
  closure
