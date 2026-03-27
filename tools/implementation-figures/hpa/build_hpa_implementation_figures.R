source("tools/implementation-figures/helpers/common.R")
impl_load_all()

fluid_sphere <- fls_generate(
  shape = sphere(radius_body = 4e-3, n_segments = 80),
  density_body = 1045,
  sound_speed_body = 1520,
  theta_body = pi / 2
)

frequency <- seq(10e3, 400e3, by = 2e3)
johnson_object <- target_strength(
  object = fluid_sphere,
  frequency = frequency,
  model = "hpa",
  method = "johnson"
)
stanton_object <- target_strength(
  object = fluid_sphere,
  frequency = frequency,
  model = "hpa",
  method = "stanton"
)

impl_with_png(
  impl_output_path("hpa", "hpa-model-plot.png"),
  {
    plot(stanton_object, type = "model")
  },
  width = 1100,
  height = 1100
)

impl_with_png(
  impl_output_path("hpa", "hpa-sphere-formulation-comparison.png"),
  {
    johnson_df <- extract(johnson_object, "model")$HPA
    stanton_df <- extract(stanton_object, "model")$HPA
    graphics::plot(
      johnson_df$frequency * 1e-3,
      johnson_df$TS,
      type = "l",
      lwd = 2.5,
      col = "black",
      xlab = "Frequency (kHz)",
      ylab = impl_ts_label()
    )
    graphics::lines(
      stanton_df$frequency * 1e-3,
      stanton_df$TS,
      col = "firebrick3",
      lwd = 2
    )
    graphics::legend(
      "bottomright",
      legend = c("Johnson", "Stanton"),
      col = c("black", "firebrick3"),
      lwd = c(2.5, 2),
      bty = "n"
    )
  }
)
