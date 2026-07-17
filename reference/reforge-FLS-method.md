# Reforge FLS-class object

Resize a fluid-like scatterer's single body by scaling its `length`
and/or `radius`, supplied either as a direct `body_scale` factor or a
`body_target` dimension (m). The cross-section is circular, so `radius`
drives both lateral and vertical extent together.

## Usage

``` r
# S4 method for class 'FLS'
reforge(
  object,
  body_scale = NULL,
  body_target = NULL,
  isometric_body = TRUE,
  n_segments_body = NULL,
  length = NULL,
  radius = NULL,
  length_radius_ratio_constant = TRUE,
  n_segments = NULL
)
```

## Arguments

- object:

  FLS-class object.

- body_scale:

  Proportional scaling to the body length and radius. When a single
  value is supplied, both dimensions are scaled together. Otherwise,
  this input must be a named numeric vector using `length` and/or
  `radius`.

- body_target:

  Target dimensions (m) for the body length and/or radius. This input
  must be a named numeric vector.

- isometric_body:

  Logical; maintain isometric scaling for body.

- n_segments_body:

  New number of segments along the body profile.

- length:

  Legacy alias for a new body length resize.

- radius:

  Legacy alias for a new maximum body radius.

- length_radius_ratio_constant:

  Legacy toggle controlling whether a length-only resize also rescales
  radius.

- n_segments:

  Legacy alias for `n_segments_body`.

## Details

For bent bodies (see
[`brake()`](https://brandynlucca.github.io/acousticTS/reference/brake.md))
a `length` resize follows the true centerline arc length and rescales
the curved path, while a `radius` resize changes only the tube thickness
and leaves the centerline in place. The legacy `length`, `radius`,
`length_radius_ratio_constant`, and `n_segments` arguments are retained
as thin wrappers over the `body_scale`/`body_target` pathway.
