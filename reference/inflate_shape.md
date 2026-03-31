# Locally widen, pinch, or taper a stored shape profile

Apply a localized scaling window to a stored profile. Values of
`scale > 1` inflate the selected region, while `scale < 1` pinch or
taper it. For canonical `Shape` objects, the result is returned as an
`Arbitrary` shape because the edited profile is no longer guaranteed to
preserve the canonical class geometry.

## Usage

``` r
inflate_shape(
  object,
  x_range = NULL,
  scale = 1,
  component = NULL,
  axis = c("radius", "width", "height"),
  profile = c("cosine", "linear", "box"),
  containment = c("warn", "error", "ignore")
)
```

## Arguments

- object:

  Shape or Scatterer object.

- x_range:

  Optional axial interval over which the local manipulation is applied.

- scale:

  Positive local scale factor.

- component:

  Optional component name for scatterers. Defaults to the primary
  geometry (`"body"` for most scatterers and `"shell"` for `ESS`).

- axis:

  Which profile dimension to scale.

- profile:

  Local window profile. `"cosine"` gives a smooth bump centered inside
  `x_range`, `"linear"` gives a triangular bump, and `"box"` applies a
  uniform factor inside the interval.

- containment:

  Containment policy used when a moved swimbladder or backbone is
  checked against its body: `"warn"`, `"error"`, or `"ignore"`.

## Value

The modified object. `Shape` inputs that are locally reshaped are
returned as `Arbitrary` shapes.

## See also

[`smooth_shape()`](https://brandynlucca.github.io/acousticTS/reference/smooth_shape.md),
[`resample_shape()`](https://brandynlucca.github.io/acousticTS/reference/resample_shape.md),
[`flip_shape()`](https://brandynlucca.github.io/acousticTS/reference/flip_shape.md),
[`brake()`](https://brandynlucca.github.io/acousticTS/reference/brake.md),
[`reforge()`](https://brandynlucca.github.io/acousticTS/reference/reforge.md)

## Examples

``` r
shape_obj <- cylinder(length_body = 0.05, radius_body = 0.003, n_segments = 80)
pinched_shape <- inflate_shape(
  shape_obj,
  x_range = c(0.015, 0.035),
  scale = 0.7
)
max(extract(pinched_shape, c("shape_parameters", "max_radius")))
#> [1] 0.003
```
