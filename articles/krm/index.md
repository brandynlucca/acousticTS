# KRM

## Overview

Benchmarked Validated

[Theory](https://brandynlucca.github.io/acousticTS/articles/krm/krm-theory.md)
[Implementation](https://brandynlucca.github.io/acousticTS/articles/krm/krm-implementation.md)

The Kirchhoff-ray mode model (`KRM`) calculates composite monostatic
backscatter from a fish body and an explicitly represented gas-filled
swimbladder.

### Core idea

Treat the segmented flesh body with a Kirchhoff-style approximation, use
an acoustic-size-dependent mode or ray description for the swimbladder,
and combine the component amplitudes coherently ([Clay
1991](#ref-Clay_1991); [Clay and Horne 1994](#ref-Clay_1994)).

### Best for

- Gas-bearing fish with separate body and swimbladder outlines
- Frequency and orientation studies where the internal gas component
  dominates or interferes with the body return
- Comparison of supported KRM material and embedding variants

### Supports

- `SBF` objects with body and swimbladder shapes, placement, and
  orientation
- Separate body and bladder material descriptions
- Stored body, swimbladder, and coherently combined outputs

### Main assumptions

- Segment-wise body and simplified swimbladder approximations
- The selected KRM variant matches the intended material interpretation
- Coherent first-order component combination
- No fully coupled body-swimbladder boundary solve or repeated
  rescattering

### Validation status

- Benchmarked against the canonical spectra stored in benchmark_ts.
- Validated against KRMr, echoSMs, and the NOAA KRM applet on bundled
  fish objects and shared workflows.

### Family pages

- [Implementation](https://brandynlucca.github.io/acousticTS/articles/krm/krm-implementation.md):
  building fish objects, variants, component output, and comparisons
- [Theory](https://brandynlucca.github.io/acousticTS/articles/krm/krm-theory.md):
  body term, swimbladder branches, and coherent sum

## References

Clay, C. S. 1991. “Low-Resolution Acoustic Scattering Models:
Fluid-Filled Cylinders and Fish with Swim Bladders.” *The Journal of the
Acoustical Society of America* 89 (5): 2168–79.
<https://doi.org/10.1121/1.400910>.

Clay, Clarence S., and John K. Horne. 1994. “Acoustic Models of Fish:
The Atlantic Cod (*Gadus Morhua*).” *The Journal of the Acoustical
Society of America* 96 (3): 1661–68. <https://doi.org/10.1121/1.410245>.
