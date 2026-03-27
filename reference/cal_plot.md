# Plotting for CAL-class objects

Plotting for CAL-class objects

## Usage

``` r
cal_plot(
  object,
  type = "shape",
  nudge_y = 1.01,
  nudge_x = 1.01,
  x_units = "frequency",
  ...
)
```

## Arguments

- object:

  CAL-class object.

- type:

  Toggle between body shape ("shape") or modeling results ("model")

- nudge_y:

  y-axis nudge.

- nudge_x:

  x-axis nudge.

- x_units:

  If "model" is selected, then toggle between frequency ("frequency",
  kHz) or ka ("ka").

- ...:

  Additional plot inputs
