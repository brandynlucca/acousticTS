# BCMS

## Overview

Unvalidated Experimental

[Theory](https://brandynlucca.github.io/acousticTS/articles/bcms/bcms-theory.md)
[Implementation](https://brandynlucca.github.io/acousticTS/articles/bcms/bcms-implementation.md)

This family follows the deformed-cylinder and coherence-corrected
cylinder literature for weakly scattering elongated bodies ([Stanton
1989](#ref-Stanton_1989_2), [1988](#ref-Stanton_1988)).

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

- BCMS is currently marked experimental because the documented checks
  are internal coherence reconstructions rather than an external
  benchmark or software-comparison ladder.

### Family pages

- [Implementation](https://brandynlucca.github.io/acousticTS/articles/bcms/bcms-implementation.md):
  usage and current reference checks
- [Theory](https://brandynlucca.github.io/acousticTS/articles/bcms/bcms-theory.md):
  bent-axis coherence integrals and relation to the straight-cylinder
  modal kernel

## References

Stanton, T. K. 1988. “Sound Scattering by Cylinders of Finite Length. I.
Fluid Cylinders.” *The Journal of the Acoustical Society of America* 83
(1): 55–63. <https://doi.org/10.1121/1.396184>.

———. 1989. “Sound Scattering by Cylinders of Finite Length. III.
Deformed Cylinders.” *The Journal of the Acoustical Society of America*
86 (2): 691–705. <https://doi.org/10.1121/1.398193>.
