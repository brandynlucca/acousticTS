# Initialize SBF-class object for KRM calculations.

Initialize SBF-class object for KRM calculations.

## Usage

``` r
krm_initialize(
  object,
  frequency,
  sound_speed_sw = 1500,
  density_sw = 1026,
  density_body = NULL,
  density_swimbladder = NULL,
  sound_speed_body = NULL,
  sound_speed_swimbladder = NULL,
  theta_body = NULL,
  theta_swimbladder = NULL
)
```

## Arguments

- object:

  SBF-class object

- frequency:

  Frequency (Hz).

- sound_speed_sw:

  Seawater sound speed.

- density_sw:

  Seawater density.

- density_body:

  Optional flesh density input.

- density_swimbladder:

  Optional gas density input.

- sound_speed_body:

  Optional flesh sound speed input.

- sound_speed_swimbladder:

  Optional gas sound speed input.

- theta_body:

  Optional orientation input (relative to incident sound wave).

- theta_swimbladder:

  Optional orientation input (relative to incident sound wave).
