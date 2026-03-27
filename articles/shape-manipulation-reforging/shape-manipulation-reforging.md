# Shape Manipulation and Reforging

## Introduction

Shape manipulation matters because many published target representations
are built from segmented or reformatted coordinate sets rather than from
canonical closed forms ([Clay and Horne 1994](#ref-clay_horne_1994);
[Gastauer
2019](#ref-gastauer_australianantarcticdivisionzooscatr_2019)).

Real workflows often begin with a target description that is close to
useful but not quite in the form required for the next model or
comparison step. The package includes tools for reshaping, bending, or
otherwise transforming existing objects rather than rebuilding them from
scratch.

These tools are not cosmetic. They are part of the modeling layer
because curvature, segmentation, and relative body-swimbladder scaling
can all change the resulting backscatter in physically meaningful ways.

That is the key idea behind this vignette. Shape manipulation is not
only about producing a different outline. It is about deciding how a
target should be represented for a specific scientific or modeling
question. A bent body, an inflated swimbladder, or a resampled
segmentation may all be legitimate transformations, but each of them
changes what the object now means physically.

## Main ideas

This part of the package is centered on utilities such as
[`brake()`](https://brandynlucca.github.io/acousticTS/reference/brake.md)
and
[`reforge()`](https://brandynlucca.github.io/acousticTS/reference/reforge.md).
Conceptually, these tools let the user keep the same target identity
while changing representation, curvature, or model-facing structure.

The scratch notes make the distinction especially clear:

1.  [`brake()`](https://brandynlucca.github.io/acousticTS/reference/brake.md)
    changes curvature while preserving the target identity,
2.  [`reforge()`](https://brandynlucca.github.io/acousticTS/reference/reforge.md)
    changes scaling, target dimensions, or resolution of an existing
    representation.

The practical distinction is simple:

1.  use
    [`brake()`](https://brandynlucca.github.io/acousticTS/reference/brake.md)
    when the main question is body curvature or pose,
2.  use
    [`reforge()`](https://brandynlucca.github.io/acousticTS/reference/reforge.md)
    when the question is dimensions, proportions, or axial resolution.

## Bending with `brake()`

[`brake()`](https://brandynlucca.github.io/acousticTS/reference/brake.md)
applies a smooth curvature transformation to an existing body
description. It can work on a body-style coordinate list or on a full
scatterer object, depending on the workflow.

The key arguments are:

1.  `radius_curvature`, the amount of curvature,
2.  `mode = "ratio"`, which interprets curvature relative to body
    length,
3.  `mode = "measurement"`, which interprets curvature in physical
    units.

Ratio mode is often the more robust choice for comparative studies
because it keeps curvature tied to target size. Measurement mode is more
natural when curvature is known from direct observation.

This distinction matters because the same nominal curvature can mean
very different things for differently sized targets. Ratio mode
preserves geometric interpretation across a size series. Measurement
mode preserves literal physical curvature. Neither is universally
better, but they answer different questions and should not be treated as
interchangeable.

``` r
library(acousticTS)

shape_obj <- cylinder(
  length_body = 0.05,
  radius_body = 0.003,
  n_segments = 120
)

obj <- fls_generate(
  shape = shape_obj,
  density_body = 1045,
  sound_speed_body = 1520
)

bent_ratio <- brake(obj, radius_curvature = 5, mode = "ratio")
bent_measured <- brake(obj, radius_curvature = 0.35, mode = "measurement")
```

After bending, it is good practice to re-plot the geometry before
running a model.

``` r
plot(bent_ratio, type = "shape")
```

That check is worth doing because a mathematically valid transformation
can still create a geometry that is too coarse for the intended model if
the segmentation is sparse.

It is also worth checking because curvature can alter more than the
silhouette. Depending on the downstream model, bending can change
projected length, local orientation, and the phase relationships among
different parts of the body. A curvature transformation is therefore
best interpreted as a new geometric state of the same target, not just
as a cosmetic deformation of a drawing.

## Re-parameterizing with `reforge()`

[`reforge()`](https://brandynlucca.github.io/acousticTS/reference/reforge.md)
is a generic interface for reshaping or resizing an existing object. The
most developed method in the current package is the `SBF` method, which
handles body and swimbladder geometry jointly.

Its most important arguments are:

1.  `body_scale` and `swimbladder_scale` for proportional resizing,
2.  `body_target` and `swimbladder_target` for direct target dimensions,
3.  `maintain_ratio` for preserving relative component proportions when
    appropriate,
4.  `n_segments_body` and `n_segments_swimbladder` for resampling,
5.  `swimbladder_inflation_factor` for controlled bladder-size changes.

There is a critical interpretation difference between scale arguments
and target-dimension arguments. Scale arguments ask, “How much larger or
smaller should this component become relative to its current state?”
Target-dimension arguments ask, “What final dimensions should this
component have?”

``` r
obj_scaled <- reforge(
  obj_sbf,
  body_scale = c(length = 1.2),
  swimbladder_scale = c(height = 0.8),
  n_segments_body = 60,
  n_segments_swimbladder = 40
)

obj_target <- reforge(
  obj_sbf,
  body_target = c(length = 0.12),
  swimbladder_target = c(length = 0.07, height = 0.0025),
  maintain_ratio = FALSE
)
```

Named vectors matter here. For anisotropic changes, dimensions should be
supplied using names such as `length`, `width`, and `height`. That keeps
the transformation explicit and avoids accidental axis mismatches.

The deeper point is that
[`reforge()`](https://brandynlucca.github.io/acousticTS/reference/reforge.md)
often changes representation more directly than
[`brake()`](https://brandynlucca.github.io/acousticTS/reference/brake.md).
Bending usually preserves the same broad anatomical proportions while
changing pose. Reforging may change size, aspect ratio, internal
proportions, and numerical resolution all at once. That makes it
especially powerful for simulation and sensitivity studies, but it also
means the transformed object should be treated as a deliberate new
parameterization rather than as a trivial variant of the original.

## Why this matters

These tools are especially useful when:

1.  an arbitrary measured geometry must be adapted for a model family,
2.  curvature or pose needs to be explored systematically,
3.  one object representation needs to be converted into another for
    comparison.

They are also useful when the scientific question itself is
morphological. For example, one may want to ask whether target strength
is more sensitive to curvature, overall length, or swimbladder
inflation. Those are exactly the kinds of questions that reshaping
utilities make tractable.

In that sense, these utilities are often not just preprocessing steps.
They can be the mechanism by which the actual hypothesis is encoded. If
a workflow is designed to ask what aspect of morphology matters most
acoustically, then the transformation function is part of the
experimental design.

## Resolution and physicality checks

Two geometry checks are especially important after
[`brake()`](https://brandynlucca.github.io/acousticTS/reference/brake.md)
or
[`reforge()`](https://brandynlucca.github.io/acousticTS/reference/reforge.md):

1.  whether the object still has adequate segment resolution for the
    intended model,
2.  whether the transformed swimbladder remains physically plausible
    inside the body.

The second point is particularly important for `SBF` workflows. If the
swimbladder is inflated or rescaled aggressively, geometric containment
can become questionable even before the code throws a warning.

For that reason, a practical reshaping workflow is usually:

1.  transform the object,
2.  plot the transformed geometry,
3.  inspect component dimensions,
4.  then re-run
    [`target_strength()`](https://brandynlucca.github.io/acousticTS/reference/target_strength.md).

## Choosing between `brake()` and `reforge()`

Although both functions manipulate geometry, they answer different
questions.

Use
[`brake()`](https://brandynlucca.github.io/acousticTS/reference/brake.md)
when:

1.  the target identity should remain the same but posture changes,
2.  curvature is the parameter of interest,
3.  the original segmentation is already suitable.

Use
[`reforge()`](https://brandynlucca.github.io/acousticTS/reference/reforge.md)
when:

1.  component lengths or widths need to be rescaled,
2.  body and swimbladder proportions must be altered separately,
3.  segment counts must be regularized for downstream modeling,
4.  simulation workflows need dimension changes as explicit parameters.

## Recommended use

Geometry transformation is best treated as a modeling decision rather
than a purely cosmetic step. Any reshaping or reforging should be
interpreted in light of the model assumptions it is intended to support.

For
[`brake()`](https://brandynlucca.github.io/acousticTS/reference/brake.md),
the practical distinction is between ratio-based curvature and curvature
specified in measured units. For
[`reforge()`](https://brandynlucca.github.io/acousticTS/reference/reforge.md),
the practical distinction is between proportional scaling and
target-dimension-based resizing. Both are useful, but they answer
different modeling questions.

One useful habit is to preserve the original object and store
transformed variants as separate named objects. That makes it easier to
compare before/after model outputs without losing the original
parameterization.

Another useful habit is to document what physical interpretation the
transformation is meant to represent. A bent object may correspond to
posture. A rescaled swimbladder may correspond to inflation state. A
resampled geometry may correspond to numerical regularization rather
than a biological change. Making that meaning explicit helps keep later
comparisons scientifically interpretable instead of turning into an
untracked series of geometric edits.

## Related reading

- [Building
  shapes](https://brandynlucca.github.io/acousticTS/articles/building-shapes/building-shapes.md)
- [Building
  scatterers](https://brandynlucca.github.io/acousticTS/articles/building-scatterers/building-scatterers.md)
- [Comparing models on the same
  target](https://brandynlucca.github.io/acousticTS/articles/comparing-models/comparing-models.md)

## References

Clay, Clarence S., and John K. Horne. 1994. “Acoustic Models of Fish:
The Atlantic Cod (*Gadus Morhua*).” *The Journal of the Acoustical
Society of America* 96 (3): 1661–68. <https://doi.org/10.1121/1.410245>.

Gastauer, Sven. 2019. “AustralianAntarcticDivision/ZooScatR: First
Package Release.” Zenodo. <https://doi.org/10.5281/ZENODO.2536635>.
