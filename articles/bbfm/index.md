# BBFM

## Overview

Experimental Unvalidated

*Model-family pages:*
[Overview](https://brandynlucca.github.io/acousticTS/articles/bbfm/index.md)
[Implementation](https://brandynlucca.github.io/acousticTS/articles/bbfm/bbfm-implementation.md)
[Theory](https://brandynlucca.github.io/acousticTS/articles/bbfm/bbfm-theory.md)

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

- BBFM currently has documented internal reconstruction checks but no
  external benchmark ladder or independent public implementation
  comparison.

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

That is why the family belongs in `acousticTS` even before a more
complete composite solver exists: it gives a physically interpretable
way to keep the two dominant anatomical components explicit in one
coherent target-strength calculation.
