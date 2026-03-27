# Initialize FLS-class object for TS modeling.

Initialize FLS-class object for TS modeling.

## Usage

``` r
dwba_curved_initialize(
  object,
  frequency,
  sound_speed_sw = 1500,
  density_sw = 1026,
  radius_curvature_ratio = NULL,
  theta = pi/2
)
```

## Arguments

- object:

  FLS-class object.

- frequency:

  Frequency (Hz).

- sound_speed_sw:

  Seawater sound speed.

- density_sw:

  Seawater density.

- radius_curvature_ratio:

  Radius of curvature ratio (length-to-curvature ratio).

- theta:

  Angle of incident soundwave (pi / 2 is broadside).
