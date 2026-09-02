# Generate a ELA-class object.

Generate a ELA-class object.

## Usage

``` r
ela_generate(
  shape,
  density_body = NULL,
  sound_speed_longitudinal_body = NULL,
  sound_speed_transversal_body = NULL,
  theta_body = pi/2,
  ID = NULL,
  length_units = "m",
  theta_units = "radians"
)
```

## Arguments

- shape:

  Pre-built `Shape` object describing the elastic target geometry.

- density_body:

  Optional body density (kg/m^3).

- sound_speed_longitudinal_body:

  Optional longitudinal wave speed (m/s).

- sound_speed_transversal_body:

  Optional transversal wave speed (m/s).

- theta_body:

  Body orientation relative to the incident wave (radians).

- ID:

  Optional metadata identifier.

- length_units:

  Compatibility argument. Scatterer constructors now assume meters and
  ignore non-SI alternatives.

- theta_units:

  Compatibility argument. Scatterer constructors now assume radians and
  ignore non-SI alternatives.

## Value

ELA-class object.

## Details

`ela_generate()` builds a generic solid-elastic scatterer under the
shared `ELA` class. Unlike `CAL`, it is not restricted to spheres, so it
can carry prolate, oblate, cylindrical, or arbitrary elastic shapes
together with the body density and longitudinal/transversal wave speeds
used by solid-elastic models such as `TMM`.

## See also

[`ELA`](https://brandynlucca.github.io/acousticTS/reference/ELA-class.md),
[`CAL`](https://brandynlucca.github.io/acousticTS/reference/CAL-class.md),
[`ESS`](https://brandynlucca.github.io/acousticTS/reference/ESS-class.md)

## Examples

``` r
elastic_sphere <- sphere(radius_body = 0.01, n_segments = 40)
ela_generate(
  shape = elastic_sphere,
  density_body = 7800,
  sound_speed_longitudinal_body = 5900,
  sound_speed_transversal_body = 3200
)
```
