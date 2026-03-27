# Method for what is printed for objects.

Method for what is printed for objects.

## Usage

``` r
plot(x, y, ...)

# S4 method for class 'Scatterer,missing'
plot(
  x,
  y,
  type = "shape",
  nudge_y = 1.1,
  nudge_x = 1.05,
  aspect_ratio = "manual",
  x_units = "frequency",
  y_units = "TS",
  ...
)
```

## Arguments

- x:

  Scatterer-class object.

- y:

  Ignored (required for plot method signature).

- ...:

  Additional plot inputs

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
