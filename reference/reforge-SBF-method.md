# Resizing function for swimbladdered targets

Resizing function for swimbladdered targets

## Usage

``` r
# S4 method for class 'SBF'
reforge(
  object,
  body_scale = NULL,
  body_target = NULL,
  swimbladder_scale = NULL,
  swimbladder_target = NULL,
  isometric_body = TRUE,
  isometric_swimbladder = TRUE,
  maintain_ratio = TRUE,
  swimbladder_inflation_factor = 1,
  n_segments_body = NULL,
  n_segments_swimbladder = NULL
)
```

## Arguments

- object:

  SBF-class object.

- body_scale:

  Proportional scaling to the body length, width, and height dimensions.
  When a single value is supplied, all dimensions are scaled using the
  same scaling factor. Otherwise, this input must be a named numeric
  vector.

- body_target:

  Target dimensions (m) for the body length, width, and height
  dimensions. This input must be a named numeric vector.

- swimbladder_scale:

  Proportional scaling to the swimbladder length, width, and height
  dimensions. When a single value is supplied, all dimensions are scaled
  using the same scaling factor. Otherwise, this input must be a named
  numeric vector.

- swimbladder_target:

  Target dimensions (m) for the swimbladder length, width, and height
  dimensions. This input must be a named numeric vector.

- isometric_body:

  Logical; maintain isometric scaling for body.

- isometric_swimbladder:

  Logical; maintain isometric scaling for bladder.

- maintain_ratio:

  Maintain size ratio between body and swimbladder.

- swimbladder_inflation_factor:

  Proportional swimbladder volume where the swimbladder x-axis origin
  and terminus are both held constant.

- n_segments_body:

  Number of segments along the body.

- n_segments_swimbladder:

  Number of segments along the bladder.
