# Re-anchor a stored shape or scatterer component along the x axis

Translate a shape or scatterer component so that its nose, center, or
tail lies at a specified x position.

## Usage

``` r
reanchor_shape(
  object,
  anchor = c("nose", "center", "tail", "max_x", "min_x"),
  at = 0,
  component = NULL,
  containment = c("warn", "error", "ignore")
)
```

## Arguments

- object:

  Shape or Scatterer object.

- anchor:

  Anchor location used to define the translation target. `nose` is
  treated as the maximum stored x position and `tail` as the minimum.

- at:

  Target x location (m).

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
[`brake()`](https://brandynlucca.github.io/acousticTS/reference/brake.md),
[`reforge()`](https://brandynlucca.github.io/acousticTS/reference/reforge.md),
[`extract()`](https://brandynlucca.github.io/acousticTS/reference/extract.md)

## Examples

``` r
shape_obj <- prolate_spheroid(length_body = 0.04, radius_body = 0.004)
centered_shape <- reanchor_shape(shape_obj, anchor = "center", at = 0)
range(extract(centered_shape, c("position_matrix", "x")))
#> [1] -0.02  0.02
```
