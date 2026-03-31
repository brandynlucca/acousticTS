# Translate a stored shape or scatterer component

Shift a `Shape` object or one geometry-bearing component of a
`Scatterer` without rebuilding it from scratch. This is most useful for
re-centering profiles, aligning a stored body to a preferred axial
origin, or nudging an internal component before model comparisons.

## Usage

``` r
translate_shape(
  object,
  x_offset = 0,
  y_offset = 0,
  z_offset = 0,
  component = NULL,
  containment = c("warn", "error", "ignore")
)
```

## Arguments

- object:

  Shape or Scatterer object.

- x_offset:

  Along-axis translation (m).

- y_offset:

  Lateral translation (m). This is ignored for profile-style geometries
  that do not store an explicit lateral centerline.

- z_offset:

  Vertical translation (m).

- component:

  Optional component name for scatterers. Defaults to the primary
  geometry (`"body"` for most scatterers and `"shell"` for `ESS`).

- containment:

  Containment policy used when a moved swimbladder or backbone is
  checked against its body: `"warn"`, `"error"`, or `"ignore"`.

## Value

The modified object, returned as the same broad object type.

## See also

[`reanchor_shape()`](https://brandynlucca.github.io/acousticTS/reference/reanchor_shape.md),
[`offset_component()`](https://brandynlucca.github.io/acousticTS/reference/offset_component.md),
[`brake()`](https://brandynlucca.github.io/acousticTS/reference/brake.md),
[`reforge()`](https://brandynlucca.github.io/acousticTS/reference/reforge.md),
[`extract()`](https://brandynlucca.github.io/acousticTS/reference/extract.md)

## Examples

``` r
shape_obj <- cylinder(length_body = 0.05, radius_body = 0.003, n_segments = 40)
moved_shape <- translate_shape(shape_obj, x_offset = 0.01)
range(extract(moved_shape, c("position_matrix", "x")))
#> [1] 0.01 0.06
```
