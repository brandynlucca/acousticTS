# Canonicalize one shape into a canonical surrogate

Build an explicit canonical surrogate shape from an existing `Shape`
object. This is intended for workflows where a measured or segmented
body should be approximated deliberately by one of the package's
canonical geometries before being passed to a shape-specific model
family such as `SPHMS`, `PSMS`, or `FCMS`.

## Usage

``` r
canonicalize_shape(
  shape,
  to = c("ProlateSpheroid", "Cylinder", "Sphere", "OblateSpheroid"),
  method = c("auto", "volume", "length_volume", "profile_l2"),
  n_segments = NULL,
  diagnostics = FALSE
)
```

## Arguments

- shape:

  A `Shape` object to approximate.

- to:

  Canonical target shape. One of `"Sphere"`, `"Cylinder"`,
  `"ProlateSpheroid"`, or `"OblateSpheroid"`.

- method:

  Canonicalization rule. `"auto"` selects a shape-specific default:
  `"volume"` for spheres, `"length_volume"` for cylinders, and
  `"profile_l2"` for spheroids. `"volume"` preserves enclosed volume
  only. `"length_volume"` preserves body length and enclosed volume.
  `"profile_l2"` fits the target canonical profile directly to the
  source equivalent-radius profile by least squares.

- n_segments:

  Optional number of segments for the returned canonical shape. Defaults
  to the source shape segment count.

- diagnostics:

  Logical. If `FALSE` (default), return only the canonical `Shape`. If
  `TRUE`, return a list containing the canonical shape and fit
  diagnostics.

## Value

If `diagnostics = FALSE`, a canonical `Shape` object. If
`diagnostics = TRUE`, a list with elements:

- `shape`: the fitted canonical `Shape`

- `diagnostics`: a named list of source, target, and fit metrics

## Details

`canonicalize_shape()` is intentionally explicit. It does **not** run
automatically inside model calls such as
`target_strength(..., model = "psms")`, because there is no single
defensible canonicalization rule for every biological or segmented
target.

When the input shape contains both an explicit width profile and
upper/lower height profiles, the canonicalization step first reduces
those local cross-sections to an equal-area circular radius before
fitting the canonical surrogate. That keeps the reduction axisymmetric
while still honoring the original local cross-sectional area more
closely than a height-only radius.

The returned diagnostics report how the source and target compare in
length, enclosed volume, maximum equivalent radius, and
equivalent-radius profile root-mean-square error on the source x-grid.
Those diagnostics are the package's recommended way to judge whether a
canonical surrogate is defensible for the intended model comparison.

## See also

[`create_shape()`](https://brandynlucca.github.io/acousticTS/reference/create_shape.md),
[`arbitrary()`](https://brandynlucca.github.io/acousticTS/reference/arbitrary.md),
[`sphere()`](https://brandynlucca.github.io/acousticTS/reference/sphere.md),
[`cylinder()`](https://brandynlucca.github.io/acousticTS/reference/cylinder.md),
[`prolate_spheroid()`](https://brandynlucca.github.io/acousticTS/reference/prolate_spheroid.md),
[`oblate_spheroid()`](https://brandynlucca.github.io/acousticTS/reference/oblate_spheroid.md)

## Examples

``` r
bladder_like <- arbitrary(
  x_body = c(0, 0.01, 0.02, 0.03, 0.04),
  radius_body = c(0, 0.004, 0.006, 0.004, 0)
)

psms_shape <- canonicalize_shape(
  bladder_like,
  to = "ProlateSpheroid"
)

cyl_fit <- canonicalize_shape(
  bladder_like,
  to = "Cylinder",
  diagnostics = TRUE
)
cyl_fit$diagnostics$fit$radius_nrmse
#> [1] 0.4567582
```
