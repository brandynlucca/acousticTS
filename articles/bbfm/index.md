# BBFM

## Overview

Unvalidated Experimental

[Theory](https://brandynlucca.github.io/acousticTS/articles/bbfm/bbfm-theory.md)
[Implementation](https://brandynlucca.github.io/acousticTS/articles/bbfm/bbfm-implementation.md)

The body-backbone fish model (`BBFM`) is a composite model for
swimbladder-less fish in which the flesh body and an elastic backbone
are retained as separate acoustic components.

### Core idea

Calculate the flesh-body contribution with `DWBA`, calculate the
backbone with `ECMS`, translate both amplitudes into one body-fixed
phase reference, and add them coherently. This separates mechanisms
emphasized in swimbladder-less fish models ([Gorska et al.
2005](#ref-Gorska_2005); [Stanton et al. 1998](#ref-Stanton_1998_1)).

### Best for

- Swimbladder-less fish with explicit body and backbone geometry
- Studying interference between weak-fluid flesh and an elastic internal
  component
- A transparent first-order composite approximation before a coupled
  solver is available

### Supports

- `BBF` objects constructed from separate body and backbone `Shape`
  objects
- Fluid-like body properties and elastic backbone density and wave
  speeds
- Stored body, backbone, phase-shifted, and combined outputs

### Main assumptions

- The flesh remains within the `DWBA` weak-scattering regime
- The backbone is represented by the canonical `ECMS` cylinder solution
- Centroid translation is sufficient to place the component phases
- No shadowing, repeated rescattering, or coupled embedded-cylinder
  boundary solve

### Validation status

- BBFM is currently marked experimental because it has documented
  internal reconstruction checks but no external benchmark ladder or
  independent public implementation comparison.

### Family pages

- [Implementation](https://brandynlucca.github.io/acousticTS/articles/bbfm/bbfm-implementation.md):
  constructing `BBF` objects, component output, and reconstruction
  checks
- [Theory](https://brandynlucca.github.io/acousticTS/articles/bbfm/bbfm-theory.md):
  component amplitudes, phase translation, and interference
- [Combining
  components](https://brandynlucca.github.io/acousticTS/articles/combining-components/combining-components.md):
  requirements for coherent and incoherent sums

## References

Gorska, Natalia, Egil Ona, and Rolf Korneliussen. 2005. “Acoustic
Backscattering by Atlantic Mackerel as Being Representative of Fish That
Lack a Swimbladder. Backscattering by Individual Fish.” *ICES Journal of
Marine Science* 62 (5): 984–95.
<https://doi.org/10.1016/j.icesjms.2005.03.010>.

Stanton, Timothy K., Dezhang Chu, Peter H. Wiebe, Linda V. Martin, and
Robert L. Eastwood. 1998. “Sound Scattering by Several Zooplankton
Groups. I. Experimental Determination of Dominant Scattering
Mechanisms.” *The Journal of the Acoustical Society of America* 103 (1):
225–35. <https://doi.org/10.1121/1.421469>.
