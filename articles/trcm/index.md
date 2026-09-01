# TRCM

## Overview

Benchmarked Unvalidated

[Theory](https://brandynlucca.github.io/acousticTS/articles/trcm/trcm-theory.md)
[Implementation](https://brandynlucca.github.io/acousticTS/articles/trcm/trcm-implementation.md)

The two-ray cylinder model (`TRCM`) is a high-frequency approximation
for fluid-like cylindrical targets whose response is dominated by two
coherent cross-sectional ray paths.

### Core idea

Retain a prompt near-side reflection and a transmitted path that returns
after one internal reflection. Their interference is multiplied by the
finite-length directivity term ([Stanton et al.
1993](#ref-Stanton_1993), [1998](#ref-Stanton_1998_2)).

### Best for

- Acoustically large, locally cylindrical fluid-like targets
- Interpreting high-frequency interference between reflected and
  transmitted paths
- Fast comparison with finite-cylinder modal calculations

### Supports

- Cylindrical `FLS` scatterers with density and sound-speed contrasts
- Straight or uniformly curved centerline geometry represented by the
  current `TRCM` interface
- Monostatic complex amplitude, cross-section, and target strength

### Main assumptions

- A high-frequency ray regime and locally circular cross-section
- Only the two leading coherent ray paths are retained
- No low-order modal resonances or full internal reverberation series
- Finite-length and curvature effects enter through approximate
  coherence factors

### Validation status

- Benchmarked within the package validation workflow against the
  canonical spectra stored in benchmark_ts. Further compared to the
  straight-cylinder and FCMS-derived bent-cylinder reference
  constructions.

### Family pages

- [Implementation](https://brandynlucca.github.io/acousticTS/articles/trcm/trcm-implementation.md):
  geometry, output, and reference comparisons
- [Theory](https://brandynlucca.github.io/acousticTS/articles/trcm/trcm-theory.md):
  ray paths, interference, and axial directivity

## References

Stanton, Timothy K., Dezhang Chu, and Peter H. Wiebe. 1998. “Sound
Scattering by Several Zooplankton Groups. II. Scattering Models.” *The
Journal of the Acoustical Society of America* 103 (1): 236–53.
<https://doi.org/10.1121/1.421110>.

Stanton, Timothy K., Dezhang Chu, Peter H. Wiebe, and Clarence S. Clay.
1993. “Average Echoes from Randomly Oriented Random-Length Finite
Cylinders: Zooplankton Models.” *The Journal of the Acoustical Society
of America* 94 (6): 3463–72. <https://doi.org/10.1121/1.407200>.
