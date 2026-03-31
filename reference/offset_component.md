# Offset an internal scatterer component without rebuilding the object

Translate an internal geometry-bearing component such as a swimbladder
or backbone while leaving the outer body unchanged.

## Usage

``` r
offset_component(
  object,
  component = c("bladder", "backbone"),
  x_offset = 0,
  z_offset = 0,
  containment = c("warn", "error", "ignore")
)
```

## Arguments

- object:

  Scatterer object containing an internal component.

- component:

  Internal component name. Currently most useful for `"bladder"` and
  `"backbone"`.

- x_offset:

  Along-axis translation (m).

- z_offset:

  Vertical translation (m).

- containment:

  Containment policy used when the shifted component is checked against
  the body: `"warn"`, `"error"`, or `"ignore"`.

## Value

The modified scatterer object.

## See also

[`translate_shape()`](https://brandynlucca.github.io/acousticTS/reference/translate_shape.md),
[`reanchor_shape()`](https://brandynlucca.github.io/acousticTS/reference/reanchor_shape.md),
[`reforge()`](https://brandynlucca.github.io/acousticTS/reference/reforge.md)

## Examples

``` r
fish <- sbf_generate(
  x_body = c(0, 0.1),
  w_body = c(0.006, 0.008),
  zU_body = c(0.001, 0.002),
  zL_body = c(-0.001, -0.002),
  x_bladder = c(0.02, 0.08),
  w_bladder = c(0, 0),
  zU_bladder = c(0.0012, 0.0012),
  zL_bladder = c(-0.0012, -0.0012),
  density_body = 1040,
  density_bladder = 1.2,
  sound_speed_body = 1500,
  sound_speed_bladder = 340
)
shifted_fish <- offset_component(fish, component = "bladder", x_offset = 0.003)
min(extract(shifted_fish, c("bladder", "rpos", "x_bladder")))
#> [1] 0.023
```
