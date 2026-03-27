# Plot scatterer geometry, stored model results, or stored TMM scattering views

S3 plotting method for `Scatterer` objects. Depending on `type`, the
method draws the target geometry, one or more stored model curves, or a
stored single-target `TMM` scattering slice / map when the object was
computed with retained T-matrix state.

## Usage

``` r
# S3 method for class 'Scatterer'
plot(
  x,
  y = NULL,
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

  Ignored (required for the base
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) method
  signature).

- type:

  Plot mode. Use `"shape"` for the stored geometry, `"model"` for the
  currently stored model output, or `"scattering"` for stored `TMM`
  angular products.

- nudge_y:

  Multiplicative vertical padding applied to automatically derived
  y-axis limits.

- nudge_x:

  Multiplicative horizontal padding applied to automatically derived
  x-axis limits.

- aspect_ratio:

  Aspect-ratio control for supported shape plots. Use `"manual"` to
  honor `nudge_x` / `nudge_y`, or `"auto"` to let the plotting helper
  derive the aspect ratio directly.

- x_units:

  Horizontal-axis convention for `type = "model"`. Supported values
  depend on the stored model family and typically include `"frequency"`
  and geometry-scaled wavenumber variants such as `"ka"`.

- y_units:

  Stored response quantity for `type = "model"` or
  `type = "scattering"`. Typical values include `"TS"` and `"sigma_bs"`.

- ...:

  Additional arguments passed through to the class-specific plotting
  helpers. For stored `TMM` scattering plots, this includes controls
  such as `frequency`, `vary`, `polar`, and `heatmap`.

## Value

Called for its side effect of drawing a plot; returns the input
invisibly.

## Details

This method dispatches to the relevant scatterer-class plotting helper:
`cal_plot()`, `ess_plot()`, `sbf_plot()`, `bbf_plot()`, `fls_plot()`, or
`gas_plot()`. The supported `type` values therefore depend on what has
been stored on the object. For example, `type = "model"` requires the
object to already contain model output, and `type = "scattering"`
currently applies only to stored `TMM` results.

## See also

[`extract()`](https://brandynlucca.github.io/acousticTS/reference/extract.md),
[`target_strength()`](https://brandynlucca.github.io/acousticTS/reference/target_strength.md),
[`tmm_scattering()`](https://brandynlucca.github.io/acousticTS/reference/tmm_scattering.md),
[`tmm_scattering_grid()`](https://brandynlucca.github.io/acousticTS/reference/tmm_scattering_grid.md)
