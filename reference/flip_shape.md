# Flip a stored shape or scatterer profile across the x or z axis

Reverse axial orientation (`axis = "x"`) or mirror the geometry
dorsoventrally (`axis = "z"`) without manually editing the position
matrix.

## Usage

``` r
flip_shape(
  object,
  axis = c("x", "z"),
  component = NULL,
  containment = c("warn", "error", "ignore")
)
```

## Arguments

- object:

  Shape or Scatterer object.

- axis:

  Flip axis. Use `"x"` to reverse nose/tail orientation while keeping
  the existing x grid, or `"z"` to mirror the profile vertically.

- component:

  Optional component name for scatterers. Defaults to the primary
  geometry (`"body"` for most scatterers and `"shell"` for `ESS`).

- containment:

  Containment policy used when a moved swimbladder or backbone is
  checked against its body: `"warn"`, `"error"`, or `"ignore"`.

## Value

The modified object, returned as the same broad object type.

## See also

[`translate_shape()`](https://brandynlucca.github.io/acousticTS/reference/translate_shape.md),
[`reanchor_shape()`](https://brandynlucca.github.io/acousticTS/reference/reanchor_shape.md),
[`inflate_shape()`](https://brandynlucca.github.io/acousticTS/reference/inflate_shape.md),
[`smooth_shape()`](https://brandynlucca.github.io/acousticTS/reference/smooth_shape.md)

## Examples

``` r
shape_obj <- arbitrary(
  x_body = c(0, 0.01, 0.02, 0.03),
  radius_body = c(0, 0.004, 0.002, 0)
)
flipped_shape <- flip_shape(shape_obj, axis = "x")
head(extract(flipped_shape, c("position_matrix", "zU")))
#> [1] 0.000 0.002 0.004 0.000
```
