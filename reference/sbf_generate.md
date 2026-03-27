# Generate a SBF-class object.

Generate a SBF-class object.

## Usage

``` r
sbf_generate(
  x_body = NULL,
  w_body = NULL,
  zU_body = NULL,
  zL_body = NULL,
  x_bladder = NULL,
  w_bladder = NULL,
  zU_bladder = NULL,
  zL_bladder = NULL,
  sound_speed_body = NULL,
  sound_speed_bladder = NULL,
  g_body = NULL,
  h_body = NULL,
  g_bladder = NULL,
  h_bladder = NULL,
  density_body = NULL,
  density_bladder = NULL,
  theta_body = pi/2,
  theta_bladder = pi/2,
  theta_units = "radians",
  length_units = "m",
  ID = NULL,
  body_shape = NULL,
  bladder_shape = NULL
)
```

## Arguments

- x_body:

  Vector containing along-body axis (m).

- w_body:

  Vector containing across-body axis (m).

- zU_body:

  Vector containing dorsal-body axis (m).

- zL_body:

  Vector containing ventral-body axis (m).

- x_bladder:

  Vector containing along-bladder axis (m).

- w_bladder:

  Vector containing across-bladder axis (m).

- zU_bladder:

  Vector containing dorsal-bladder axis (m).

- zL_bladder:

  Vector containing ventral-bladder axis (m).

- sound_speed_body:

  Flesh sound speed (c;_(body), m s⁻¹).

- sound_speed_bladder:

  Bladder sound speed (c, m \\s^-1\\).

- g_body:

  Body density contrast.

- h_body:

  Body sound speed contrast.

- g_bladder:

  Bladder density contrast.

- h_bladder:

  Bladder sound speed contrast.

- density_body:

  Flesh density (ρ_(body), kg m³).

- density_bladder:

  Bladder density (\\\rho\\, kg m³).

- theta_body:

  Angle of body relative to wavefront (\\\theta_body\\, radians).

- theta_bladder:

  Angle of bladder relative to wavefront (\\\theta_bladder\\, radians).

- theta_units:

  Compatibility argument. Scatterer constructors now assume radians and
  ignore non-SI alternatives.

- length_units:

  Compatibility argument. Scatterer constructors now assume meters and
  ignore non-SI alternatives.

- ID:

  Optional metadata identifier.

- body_shape:

  Optional pre-built Shape for the body; overrides
  x_body/w_body/zU_body/zL_body if supplied.

- bladder_shape:

  Optional pre-built Shape for the bladder; overrides
  x_bladder/w_bladder/zU_bladder/zL_bladder if supplied.

## Value

Generates a SBF-class object.

## Details

The recommended interface is to supply pre-built `Shape` objects through
`body_shape` and `bladder_shape`, then specify material properties
through either contrasts (`g_*`, `h_*`) or absolute density/sound-speed
values. Legacy coordinate-vector inputs are retained for backward
compatibility and are converted internally to `Shape` objects before the
scatterer is built. Character-based shape dispatch is deprecated; build
the component geometry first with a `Shape` constructor and then pass
the resulting `Shape` object.

Body/bladder width vectors (`w_body`, `w_bladder`) default to zeros when
missing. Resulting `rpos` matrices always carry the row names x_body,
w_body, zU_body, zL_body for the body and x_bladder, w_bladder,
zU_bladder, zL_bladder for the bladder.

## See also

[`SBF`](https://brandynlucca.github.io/acousticTS/reference/SBF-class.md)

## Examples

``` r
# Manual body/bladder coordinates
sbf_generate(
  x_body = seq(0, 0.1, length.out = 5),
  w_body = rep(0, 5),
  zU_body = seq(0.001, 0.002, length.out = 5),
  zL_body = -seq(0.001, 0.002, length.out = 5),
  x_bladder = seq(0.02, 0.09, length.out = 4),
  w_bladder = rep(0, 4),
  zU_bladder = rep(0.0015, 4),
  zL_bladder = rep(-0.0015, 4),
  sound_speed_body = 1500,
  sound_speed_bladder = 340,
  density_body = 1040,
  density_bladder = 1.2
)
#> SBF-object
#>  Swimbladdered fish (SBF) 
#>  ID:UID
#> Body dimensions:
#>  Length:0.1 m(n = 4 cylinders)
#>  Mean radius:0.0015 m
#>  Max radius:0.002 m
#> Bladder dimensions:
#>  Length:0.07 m(n = 3 cylinders)
#>  Mean radius:0.0015 m
#>  Max radius:0.0015 m
#> Body material properties:
#>  Density: 1040 kg m^-3 | Sound speed: 1500 m s^-1
#> Bladder fluid material properties:
#>  Density: 1.2 kg m^-3 | Sound speed: 340 m s^-1
#> Body orientation (relative to transducer face/axis):1.571 radians

# Using pre-built shapes
body_shape <- arbitrary(
  x_body = c(0, 0.1), zU_body = c(0.001, 0.002),
  zL_body = c(-0.001, -0.002)
)
bladder_shape <- arbitrary(
  x_bladder = c(0.02, 0.09), w_bladder = c(0, 0),
  zU_bladder = c(0.0015, 0.0015),
  zL_bladder = c(-0.0015, -0.0015)
)
sbf_generate(
  body_shape = body_shape, bladder_shape = bladder_shape,
  sound_speed_body = 1500, sound_speed_bladder = 340,
  density_body = 1040, density_bladder = 1.2
)
#> SBF-object
#>  Swimbladdered fish (SBF) 
#>  ID:UID
#> Body dimensions:
#>  Length:0.1 m(n = 1 cylinders)
#>  Mean radius:0.0015 m
#>  Max radius:0.002 m
#> Bladder dimensions:
#>  Length:0.07 m(n = 1 cylinders)
#>  Mean radius:0.0015 m
#>  Max radius:0.0015 m
#> Body material properties:
#>  Density: 1040 kg m^-3 | Sound speed: 1500 m s^-1
#> Bladder fluid material properties:
#>  Density: 1.2 kg m^-3 | Sound speed: 340 m s^-1
#> Body orientation (relative to transducer face/axis):1.571 radians
```
