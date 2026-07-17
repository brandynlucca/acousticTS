# Reforge SBF-class object

Resize a swimbladdered fish scatterer by rescaling the flesh body and
the gas-filled swimbladder together or independently. Body and
swimbladder geometry each support `length`, `width`, and `height`
scaling supplied as a direct `*_scale` factor or a `*_target` dimension
(m).

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
  n_segments_swimbladder = NULL,
  containment = c("warn", "error", "ignore")
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

- containment:

  Containment policy for internal geometry checks. Use `"warn"` to keep
  the current warning behavior, `"error"` to fail fast for invalid
  internal geometries, or `"ignore"` to skip containment checks.

## Details

Vertical scaling is applied about each component's own centerline so
that curved bodies keep their proportions when resized rather than
warping - a resized body remains geometrically similar whether it is
enlarged or shrunk. After scaling, the swimbladder is re-nested inside
the body: its relative along-body start and its relative vertical offset
within the body envelope are restored, so the bladder tracks the body
instead of drifting out of it. The `swimbladder_inflation_factor` scales
only the bladder width and height (its length endpoints are held fixed)
about the bladder's own center. The `containment` policy governs what
happens when the swimbladder exceeds the body envelope.
