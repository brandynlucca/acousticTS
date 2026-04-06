# BBFM

## Overview

Unvalidated Experimental

[Theory](https://brandynlucca.github.io/acousticTS/articles/bbfm/bbfm-theory.md)
[Implementation](https://brandynlucca.github.io/acousticTS/articles/bbfm/bbfm-implementation.md)

This family is best read alongside the swimbladder-less fish and
composite-scatterer literature that motivates explicit flesh-body and
backbone terms ([Gorska, Ona, and Korneliussen 2005](#ref-Gorska_2005);
[Stanton et al. 1998](#ref-Stanton_1998_1); [Clay and Horne
1994](#ref-Clay_1994)).

The body-backbone fish model (`BBFM`) is the package’s composite
swimbladder-less family for targets whose flesh body and backbone should
remain explicit acoustic components.

### Core idea

Compute a weak-fluid flesh-body term with `DWBA`, compute an elastic
backbone term with `ECMS`, place the backbone inside the same body-fixed
frame with a two-way phase factor, and then sum the two complex
amplitudes coherently.

### Best for

- Swimbladder-less fish where the backbone should remain acoustically
  explicit
- Composite body-plus-backbone studies that are too structured for a
  body-only approximation
- Intermediate modeling between weak-fluid body models and future fully
  coupled composite solvers

### Supports

- `BBF` scatterers built from explicit `body_shape` and `backbone_shape`
  inputs
- Flesh-body material properties as contrasts or absolute density/sound
  speed
- Backbone density plus longitudinal and transverse elastic wave speeds

### Main assumptions

- Flesh body treated in the weak-fluid `DWBA` regime
- Backbone treated as an elastic cylinder through `ECMS`
- Components combined coherently after centroid-based phase placement
- No repeated rescattering or fully coupled embedded elastic-cylinder
  solve

### Validation status

- BBFM is currently marked experimental because it has documented
  internal reconstruction checks but no external benchmark ladder or
  independent public implementation comparison.

### Family pages

- [Implementation](https://brandynlucca.github.io/acousticTS/articles/bbfm/bbfm-implementation.md):
  object setup, stored outputs, and internal reconstruction checks
- [Theory](https://brandynlucca.github.io/acousticTS/articles/bbfm/bbfm-theory.md):
  composite amplitude, phase placement, and interference structure

`BBFM` is already useful as a transparent composite family, but it
should still be treated as experimental. The main open work is external
validation and, in the longer term, a more tightly coupled treatment of
the backbone as an embedded elastic structure rather than a positioned
component surrogate.

### Why it exists

`BBFM` fills the gap between:

- body-only weak-scattering models that ignore the backbone,
- and future fully coupled body-plus-backbone solvers that are not yet
  in the package.

That is why the family belongs in acousticTS even before a more complete
composite solver exists: it gives a physically interpretable way to keep
the two dominant anatomical components explicit in one coherent
target-strength calculation.

## References

Clay, Clarence S., and John K. Horne. 1994. “Acoustic Models of Fish:
The Atlantic Cod (*Gadus Morhua*).” *The Journal of the Acoustical
Society of America* 96 (3): 1661–68. <https://doi.org/10.1121/1.410245>.

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
