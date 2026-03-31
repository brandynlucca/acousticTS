# Reforge FLS-class object.

Reforge FLS-class object.

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
