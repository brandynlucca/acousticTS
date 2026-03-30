# TRCM

## Overview

Benchmarked Unvalidated

[Theory](https://brandynlucca.github.io/acousticTS/articles/trcm/trcm-theory.md)
[Implementation](https://brandynlucca.github.io/acousticTS/articles/trcm/trcm-implementation.md)

These pages come from the high-frequency elongated-body literature and
later fish and zooplankton applications ([Stanton et al.
1993](#ref-stanton_etal_1993), [1998](#ref-stanton_sound_1998);
[Stanton, Chu, and Wiebe 1998](#ref-stanton_sound_1998-1)).

The two-ray cylinder model (`TRCM`) is a high-frequency asymptotic
family for elongated fluid-like bodies. It retains only two dominant
coherent paths: a prompt near-side reflection and a through-body path
that returns after an internal reflection.

### Core idea

Replace the full internal reverberation problem by the first two
physically important ray paths and compute the target strength from
their interference, combined with a finite-length directivity factor.

### Best for

- High-frequency elongated fluid-like targets
- Cylinder-like bodies whose scattering is dominated by specular path
  interference
- Fast asymptotic estimates when a full modal or perturbative solve is
  unnecessary

### Supports

- Straight and bent cylindrical-style branches
- Contrast notation relative to seawater as medium `1`
- Monostatic high-frequency backscatter estimates

### Main assumptions

- High-frequency regime
- Only two dominant coherent internal/external paths are retained
- No low-order resonances or full internal multiple scattering series

### Validation status

- Benchmarked within the package validation workflow against the
  straight-cylinder and FCMS-derived bent-cylinder reference
  constructions.

### Family pages

- [Implementation](https://brandynlucca.github.io/acousticTS/articles/trcm/trcm-implementation.md):
  usage and comparison workflows
- [Theory](https://brandynlucca.github.io/acousticTS/articles/trcm/trcm-theory.md):
  two-ray geometry, interference factor, and finite-length directivity

## References

Stanton, Timothy K., Dezhang Chu, and Peter H. Wiebe. 1998. “Sound
Scattering by Several Zooplankton Groups. II. Scattering Models.” *The
Journal of the Acoustical Society of America* 103 (1): 236–53.
<https://doi.org/10.1121/1.421110>.

Stanton, Timothy K., Dezhang Chu, Peter H. Wiebe, and Clarence S. Clay.
1993. “Average Echoes from Randomly Oriented Random‐length Finite
Cylinders: Zooplankton Models.” *The Journal of the Acoustical Society
of America* 94 (6): 3463–72. <https://doi.org/10.1121/1.407200>.

Stanton, Timothy K., Dezhang Chu, Peter H. Wiebe, Linda V. Martin, and
Robert L. Eastwood. 1998. “Sound Scattering by Several Zooplankton
Groups. I. Experimental Determination of Dominant Scattering
Mechanisms.” *The Journal of the Acoustical Society of America* 103 (1):
225–35. <https://doi.org/10.1121/1.421469>.
