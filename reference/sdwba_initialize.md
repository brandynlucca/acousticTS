# Initialize FLS-class object for SDWBA modeling

Initialize FLS-class object for SDWBA modeling

## Usage

``` r
sdwba_initialize(
  object,
  frequency,
  sound_speed_sw = 1500,
  density_sw = 1026,
  n_iterations = 100,
  n_segments_init = 14,
  phase_sd_init = sqrt(2)/2,
  length_init = 0.03835,
  frequency_init = 120000
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

- n_iterations:

  Number of iterations to repeat SDWBA

- n_segments_init:

  Reference number of body segments

- phase_sd_init:

  Reference phase deviation

- length_init:

  Reference body length

- frequency_init:

  Reference frequency
