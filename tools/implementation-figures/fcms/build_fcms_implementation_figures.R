source("tools/implementation-figures/helpers/common.R")
impl_load_all()

make_fcms_object <- function(boundary = "liquid_filled") {
  obj <- fls_generate(
    shape = cylinder(
      length_body = 50e-3,
      radius_body = 5e-3,
      n_segments = 80
    ),
    density_body = 1045,
    sound_speed_body = 1520,
    theta_body = pi / 2
  )

  target_strength(
    object = obj,
    frequency = seq(10e3, 400e3, by = 1e3),
    model = "fcms",
    boundary = boundary
  )
}

liquid_object <- make_fcms_object("liquid_filled")
rigid_object <- make_fcms_object("fixed_rigid")

impl_with_png(
  impl_output_path("fcms", "fcms-shape-plot.png"),
  {
    plot(liquid_object, type = "shape")
  },
  width = 1100,
  height = 1100
)

impl_with_png(
  impl_output_path("fcms", "fcms-model-plot.png"),
  {
    plot(liquid_object, type = "model")
  },
  width = 1100,
  height = 1100
)

impl_with_png(
  impl_output_path("fcms", "fcms-boundary-comparison.png"),
  {
    liquid_df <- extract(liquid_object, "model")$FCMS
    rigid_df <- extract(rigid_object, "model")$FCMS
    graphics::plot(
      liquid_df$frequency * 1e-3,
      liquid_df$TS,
      type = "l",
      lwd = 2.5,
      col = "black",
      xlab = "Frequency (kHz)",
      ylab = impl_ts_label(),
      ylim = c(-105, -20)
    )
    graphics::lines(
      rigid_df$frequency * 1e-3,
      rigid_df$TS,
      col = "firebrick3",
      lwd = 2
    )
    graphics::legend(
      "bottomleft",
      legend = c("Liquid-filled", "Fixed rigid"),
      col = c("black", "firebrick3"),
      lwd = c(2.5, 2),
      bty = "n"
    )
  }
)
