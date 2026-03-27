# Create GAS object

Create GAS object

## Usage

``` r
gas_generate(
  shape = "sphere",
  radius_body = NULL,
  h_fluid = 0.22,
  g_fluid = 0.0012,
  sound_speed_fluid = NULL,
  density_fluid = NULL,
  theta_body = pi/2,
  ID = NULL,
  radius_units = "m",
  theta_units = "radians",
  n_segments = 100
)
```

## Arguments

- shape:

  Optional pre-made shape input. Default is a sphere.

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

  Orientation of the target relative to the transmit source
  (\\\theta\\). Broadside incidence is considered 90 degrees, or pi/2.
  Default value is pi/2; input should be in radians.

- ID:

  Optional metadata entry.

- radius_units:

  Diameter units. Defaults to "m".

- theta_units:

  Units used for orientation. Defaults to "radians".

- n_segments:

  Number of body segments.

## Value

GAS-class object

## See also

[`GAS`](https://brandynlucca.github.io/acousticTS/reference/GAS-class.md)
