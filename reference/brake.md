# Bend a scatterer body or body component

Apply a smooth curvature transformation to an existing scatterer body or
to a list-like body component containing an `rpos` matrix. This is
useful when a target should keep the same broad identity while adopting
a curved centerline for model comparisons or sensitivity studies.

## Usage

``` r
brake(object, radius_curvature, mode = "ratio")
```

## Arguments

- object:

  Dataframe or scatterer-class object

- radius_curvature:

  Radius of curvature that can be parameterized either as a ratio
  relative to body length or actual measurement

- mode:

  Either "ratio" or "measurement"

## Value

A bent version of `object`, returned as the same broad object type with
updated geometry and curvature metadata.

## See also

[`extract()`](https://brandynlucca.github.io/acousticTS/reference/extract.md),
[`reforge()`](https://brandynlucca.github.io/acousticTS/reference/reforge.md),
[`translate_shape()`](https://brandynlucca.github.io/acousticTS/reference/translate_shape.md),
[`reanchor_shape()`](https://brandynlucca.github.io/acousticTS/reference/reanchor_shape.md),
[`inflate_shape()`](https://brandynlucca.github.io/acousticTS/reference/inflate_shape.md),
[`smooth_shape()`](https://brandynlucca.github.io/acousticTS/reference/smooth_shape.md),
[`resample_shape()`](https://brandynlucca.github.io/acousticTS/reference/resample_shape.md),
[`flip_shape()`](https://brandynlucca.github.io/acousticTS/reference/flip_shape.md),
[`offset_component()`](https://brandynlucca.github.io/acousticTS/reference/offset_component.md)

## Examples

``` r
shape_obj <- cylinder(
  length_body = 0.05,
  radius_body = 0.003,
  n_segments = 80
)
obj <- fls_generate(
  shape = shape_obj,
  density_body = 1045,
  sound_speed_body = 1520
)

bent_obj <- brake(obj, radius_curvature = 5)
head(extract(bent_obj, c("body", "rpos", "z")))
#> [1] -0.0012489587 -0.0011873402 -0.0011272768 -0.0010687689 -0.0010118167
#> [6] -0.0009564208
extract(bent_obj, c("shape_parameters", "radius_curvature_ratio"))
#> [1] 5

bent_body <- brake(extract(obj, "body"), radius_curvature = 0.35, mode = "measurement")
head(bent_body$rpos["z", ])
#> [1] -0.0008924776 -0.0008484293 -0.0008054944 -0.0007636730 -0.0007229653
#> [6] -0.0006833713
```
