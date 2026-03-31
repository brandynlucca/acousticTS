# Reforge BBF-class object.

Resize a backboned-fish scatterer by rescaling the flesh body and
elastic backbone independently or together. The body follows the same
profile-based length/width/height scaling used by `SBF`, while the
backbone retains the same cylinder-style component bookkeeping and
preserves its relative axial start within the body when body length
changes.

## Usage

``` r
# S4 method for class 'BBF'
reforge(
  object,
  body_scale = NULL,
  body_target = NULL,
  backbone_scale = NULL,
  backbone_target = NULL,
  isometric_body = TRUE,
  isometric_backbone = TRUE,
  maintain_ratio = TRUE,
  n_segments_body = NULL,
  n_segments_backbone = NULL,
  containment = c("warn", "error", "ignore")
)
```

## Arguments

- object:

  BBF-class object.

- body_scale:

  Proportional scaling to the body length, width, and height dimensions.
  When a single value is supplied, all dimensions are scaled uniformly.
  Otherwise, this input must be a named numeric vector.

- body_target:

  Target dimensions (m) for the body length, width, and height
  dimensions. This input must be a named numeric vector.

- backbone_scale:

  Proportional scaling to the backbone length, width, and height
  dimensions. When a single value is supplied, all dimensions are scaled
  uniformly. Otherwise, this input must be a named numeric vector.

- backbone_target:

  Target dimensions (m) for the backbone length, width, and height
  dimensions. This input must be a named numeric vector.

- isometric_body:

  Logical; maintain isometric scaling for body.

- isometric_backbone:

  Logical; maintain isometric scaling for backbone.

- maintain_ratio:

  Maintain size ratio between body and backbone.

- n_segments_body:

  Number of points along the body profile.

- n_segments_backbone:

  Number of points along the backbone profile.

- containment:

  Containment policy for internal geometry checks. Use `"warn"` to keep
  the current warning behavior, `"error"` to fail fast for invalid
  internal geometries, or `"ignore"` to skip containment checks.

## Value

Modified BBF-class object.
