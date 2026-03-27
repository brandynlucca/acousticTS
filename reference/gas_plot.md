# Plotting for GAS-class objects

Plotting for GAS-class objects

## Usage

``` r
gas_plot(
  object,
  type = "shape",
  nudge_y = 1.01,
  nudge_x = 1.01,
  x_units = "frequency",
  y_units = "TS",
  ...
)
```

## Arguments

- object:

  GAS-class object.

- type:

  Toggle between body shape ("shape") or modeling results ("model")

- nudge_y:

  y-axis nudge.

- nudge_x:

  x-axis nudge.

- x_units:

  If "model" is selected, then toggle between frequency ("frequency",
  kHz) or ka ("ka").

- y_units:

  y-axis data selection (e.g. TS, sigma_bs – defaults to TS).

- ...:

  Additional plot inputs
