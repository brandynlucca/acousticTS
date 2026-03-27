# Initialize GAS-object for modal series solution.

Initialize GAS-object for modal series solution.

## Usage

``` r
mss_anderson_initialize(
  object,
  frequency,
  radius = NULL,
  g_body = NULL,
  h_body = NULL,
  sound_speed_sw = 1500,
  density_sw = 1026,
  ka_limit = NULL
)
```

## Arguments

- object:

  GAS-class object.

- frequency:

  Frequency (Hz).

- radius:

  Radius of sphere (m).

- g_body:

  Density contrast for gas.

- h_body:

  Sound speed contrast for gas.

- sound_speed_sw:

  Seawater sound speed.

- density_sw:

  Seawater density.

- ka_limit:

  Modal series limit (i.e. max "m"). The default is the maximum ka + 10.
