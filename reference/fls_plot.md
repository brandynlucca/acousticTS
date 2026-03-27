# Plotting for FLS-class objects

Plotting for FLS-class objects

## Usage

``` r
fls_plot(
  object,
  type = "shape",
  nudge_y = 1.05,
  nudge_x = 1.01,
  aspect_ratio = "manual",
  x_units = "frequency",
  y_units = "TS",
  ...
)
```

## Arguments

- object:

  FLS-class object.

- type:

  Toggle between body shape ("shape") or modeling results ("model")

- nudge_y:

  y-axis nudge.

- nudge_x:

  x-axis nudge.

- aspect_ratio:

  Aspect ratio setting ( defaults to "manual" for nudge_y and nudge_x to
  apply; otherwise, input "auto").

- x_units:

  If "model" is selected, then toggle between frequency ("frequency",
  kHz) or ka ("ka").

- y_units:

  y-axis data selection (e.g. TS, sigma_bs – defaults to TS)

- ...:

  Additional plot inputs
