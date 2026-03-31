# Smooth a stored shape or scatterer profile

Smooth the stored geometry using a centered moving-average filter
applied to profile coordinates. This is useful for cleaning digitized
outlines before canonicalization or model runs. For canonical `Shape`
objects, the result is returned as an `Arbitrary` shape because the
edited profile is no longer guaranteed to preserve the canonical class
geometry.

## Usage

``` r
smooth_shape(
  object,
  span = 5,
  component = NULL,
  preserve_ends = TRUE,
  containment = c("warn", "error", "ignore")
)
```

## Arguments

- object:

  Shape or Scatterer object.

- span:

  Centered moving-average span. Even values are rounded up to the next
  odd integer.

- component:

  Optional component name for scatterers. Defaults to the primary
  geometry (`"body"` for most scatterers and `"shell"` for `ESS`).

- preserve_ends:

  Logical; whether to keep the first and last profile points fixed.

- containment:

  Containment policy used when a moved swimbladder or backbone is
  checked against its body: `"warn"`, `"error"`, or `"ignore"`.

## Value

The modified object. `Shape` inputs that are smoothed are returned as
`Arbitrary` shapes.

## See also

[`inflate_shape()`](https://brandynlucca.github.io/acousticTS/reference/inflate_shape.md),
[`resample_shape()`](https://brandynlucca.github.io/acousticTS/reference/resample_shape.md),
[`brake()`](https://brandynlucca.github.io/acousticTS/reference/brake.md),
[`reforge()`](https://brandynlucca.github.io/acousticTS/reference/reforge.md)

## Examples

``` r
shape_obj <- arbitrary(
  x_body = c(0, 0.01, 0.02, 0.03, 0.04),
  zU_body = c(0, 0.003, 0.006, 0.0035, 0),
  zL_body = c(0, -0.0025, -0.0055, -0.003, 0)
)
smoothed_shape <- smooth_shape(shape_obj, span = 3)
extract(smoothed_shape, c("shape_parameters", "n_segments"))
#> [1] 4
```
