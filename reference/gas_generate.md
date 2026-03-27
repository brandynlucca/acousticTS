# Generate a GAS-class object

Generate a GAS-class object

## Usage

``` r
gas_generate(shape = NULL, radius_body = NULL, h_fluid = 0.22,
  g_fluid = 0.0012, sound_speed_fluid = NULL, density_fluid = NULL,
  theta_body = pi/2, ID = NULL, radius_units = "m",
  theta_units = "radians", n_segments = 100, ...)
```

## Arguments

- shape:

  Pre-built `Shape` object describing the gas-filled geometry. If
  omitted, explicit profile coordinates such as `x_body`, `y_body`, and
  `z_body` are treated as the manual geometry pathway. Legacy character
  dispatch such as `"sphere"` is retained only for backward
  compatibility and is now deprecated.

- radius_body:

  Radius (m). For non-canonical shapes, this would be the maximum or
  mean radius at the scatterer midsection.

- h_fluid:

  Sound speed contrast of fluid relative to surrounding medium (h).

- g_fluid:

  Density contrast of fluid relative to surrounding density (g).

- sound_speed_fluid:

  Optional fluid sound speed (m/s).

- density_fluid:

  Optional fluid density (m/s).

- theta_body:

  Orientation of the target relative to the incident wave (radians).

- ID:

  Optional metadata identifier.

- radius_units:

  Compatibility argument. `gas_generate()` now assumes meters and
  ignores non-SI alternatives.

- theta_units:

  Compatibility argument. Scatterer constructors now assume radians and
  ignore non-SI alternatives.

- n_segments:

  Number of body segments.

- ...:

  Additional manual profile arguments or legacy canonical shape
  arguments used by the compatibility geometry pathway, such as
  `x_body`, `y_body`, `z_body`, `radius_body`, `length_body`, or
  `taper`.

## Value

GAS-class object

## Details

The preferred workflow is to supply a pre-built `Shape` object or
explicit profile coordinates and then describe the internal gas by
either contrasts (`g_fluid`, `h_fluid`) or absolute density/sound-speed
values (`density_fluid`, `sound_speed_fluid`). Character-based shape
dispatch remains available only as a compatibility pathway and is now
deprecated.

Scatterer constructors store geometry in meters and orientations in
radians. `radius_units` and `theta_units` are retained as compatibility
arguments, but non-SI values are normalized to the package-standard
representation.

## See also

[`GAS`](https://brandynlucca.github.io/acousticTS/reference/GAS-class.md)

## Examples

``` r
shape_gas <- sphere(radius_body = 0.01, n_segments = 60)
gas_generate(shape = shape_gas, g_fluid = 0.0012, h_fluid = 0.22)
#> GAS-object
#>  Gas- and fluid-filled scatterer 
#>  ID:UID
#> Body dimensions:
#>  Diameter:0.02 m
#>  Radius:0.01 m
#> Material properties:
#>  g: 0.0012
#>  h: 0.22
```
