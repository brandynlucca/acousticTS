# Resample a stored shape or scatterer profile to a new segment count

Re-discretize a shape or scatterer component along its x axis without
rebuilding it from scratch.

## Usage

``` r
resample_shape(
  object,
  n_segments,
  component = NULL,
  containment = c("warn", "error", "ignore")
)
```

## Arguments

- object:

  Shape or Scatterer object.

- n_segments:

  New number of intervals used to discretize the profile.

- component:

  Optional component name for scatterers. Defaults to the primary
  geometry (`"body"` for most scatterers and `"shell"` for `ESS`).

- containment:

  Containment policy used when a moved swimbladder or backbone is
  checked against its body: `"warn"`, `"error"`, or `"ignore"`.

## Value

The modified object, returned as the same broad object type.

## See also

[`smooth_shape()`](https://brandynlucca.github.io/acousticTS/reference/smooth_shape.md),
[`inflate_shape()`](https://brandynlucca.github.io/acousticTS/reference/inflate_shape.md),
[`reforge()`](https://brandynlucca.github.io/acousticTS/reference/reforge.md),
[`extract()`](https://brandynlucca.github.io/acousticTS/reference/extract.md)

## Examples

``` r
shape_obj <- sphere(radius_body = 0.01, n_segments = 20)
shape_fine <- resample_shape(shape_obj, n_segments = 80)
extract(shape_fine, c("shape_parameters", "n_segments"))
#> [1] 80
```
