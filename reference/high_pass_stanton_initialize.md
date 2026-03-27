# Initialize object for Stanton high-pass approximation

Initialize object for Stanton high-pass approximation

## Usage

``` r
high_pass_stanton_initialize(
  object,
  frequency,
  radius_shell,
  g_shell,
  h_shell,
  sound_speed_sw = 1500,
  density_sw = 1026
)
```

## Arguments

- object:

  ESS-class object.

- frequency:

  Frequency (Hz).

- radius_shell:

  Radius of shell.

- g_shell:

  Optional shell density contrast.

- h_shell:

  Optional shell sound speed contrast.

- sound_speed_sw:

  Seawater sound speed.

- density_sw:

  Seawater density.
