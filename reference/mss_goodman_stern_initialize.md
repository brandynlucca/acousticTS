# Initialize ESS-object for modal series solution.

Initialize ESS-object for modal series solution.

## Usage

``` r
mss_goodman_stern_initialize(
  object,
  frequency,
  sound_speed_sw = 1500,
  density_sw = 1026,
  m_limit = NULL
)
```

## Arguments

- object:

  ESS-class object.

- frequency:

  Frequency vector (Hz).

- sound_speed_sw:

  Seawater sound speed (m/s). Default: 1500.

- density_sw:

  Seawater density (kg/m³). Default: 1026.

- m_limit:

  Modal series limit (i.e. max "m"). The default is the maximum ka + 10.
