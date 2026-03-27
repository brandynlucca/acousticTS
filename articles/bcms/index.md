# BCMS

## Overview

Experimental Unvalidated

*Model-family pages:*
[Overview](https://brandynlucca.github.io/acousticTS/articles/bcms/index.md)
[Implementation](https://brandynlucca.github.io/acousticTS/articles/bcms/bcms-implementation.md)
[Theory](https://brandynlucca.github.io/acousticTS/articles/bcms/bcms-theory.md)

The bent cylinder modal series solution (`BCMS`) extends the straight
finite-cylinder modal family to uniformly bent cylinders by keeping the
straight cross-sectional modal physics and modifying only the along-axis
coherence.

### Core idea

Start from the finite-cylinder modal backscatter of a straight cylinder,
then replace the straight-axis coherent length by the curved equivalent
coherent length derived for a uniformly bent axis.

### Best for

- Uniformly bent fluid-like cylinders near broadside
- Curvature studies where the straight cylinder modal content remains
  the right local kernel
- Problems where a full curved-body rederivation is unnecessary or
  unavailable

### Supports

- Bent `Cylinder` geometries represented through a curvature ratio
- Monostatic near-broadside target strength
- The same fluid boundary families carried by the straight-cylinder
  modal kernel

### Main assumptions

- Uniform curvature
- Near-broadside incidence
- Curvature modifies coherence along the axis but not the local
  cross-sectional modal physics

### Validation status

- Current BCMS checks are internal coherence reconstructions; the family
  does not yet have an external benchmark or software ladder.

### Family pages

- [Implementation](https://brandynlucca.github.io/acousticTS/articles/bcms/bcms-implementation.md):
  usage and current reference checks
- [Theory](https://brandynlucca.github.io/acousticTS/articles/bcms/bcms-theory.md):
  bent-axis coherence integrals and relation to the straight-cylinder
  modal kernel
