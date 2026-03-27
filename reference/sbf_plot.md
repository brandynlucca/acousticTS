# Plotting for SBF-class objects

Plotting for SBF-class objects

## Usage

``` r
sbf_plot(
  object,
  type = "shape",
  nudge_y = 1.05,
  nudge_x = 1.01,
  x_units = "frequency"
)
```

## Arguments

- object:

  SBF-class object.

- type:

  Toggle between body shape ("shape") or modeling results ("model")

- nudge_y:

  y-axis nudge.

- nudge_x:

  x-axis nudge.

- x_units:

  If "model" is selected, then toggle between frequency ("frequency",
  kHz) or ka ("ka").
