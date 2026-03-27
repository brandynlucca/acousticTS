source("tools/implementation-figures/helpers/common.R")
impl_load_all()

make_sphms_object <- function(boundary = "gas_filled") {
  obj <- gas_generate(
    shape = sphere(radius_body = 1e-3),
    g_fluid = 0.0012,
    h_fluid = 0.22,
    theta_body = pi / 2
  )

  target_strength(
    object = obj,
    frequency = seq(1e3, 100e3, by = 1e3),
    model = "sphms",
    boundary = boundary
  )
}

gas_object <- make_sphms_object("gas_filled")
rigid_object <- make_sphms_object("fixed_rigid")
pr_object <- make_sphms_object("pressure_release")

impl_with_png(
  impl_output_path("sphms", "sphms-shape-plot.png"),
  {
    plot(gas_object, type = "shape")
  },
  width = 1100,
  height = 1100
)

impl_with_png(
  impl_output_path("sphms", "sphms-model-plot.png"),
  {
    plot(gas_object, type = "model")
  },
  width = 1100,
  height = 1100
)

impl_with_png(
  impl_output_path("sphms", "sphms-boundary-comparison.png"),
  {
    models <- list(
      `Fixed rigid` = rigid_object,
      `Pressure release` = pr_object,
      `Gas filled` = gas_object
    )
    colors <- c("black", "firebrick3", "dodgerblue3")
    graphics::plot(
      extract(models[[1]], "model")$SPHMS$frequency * 1e-3,
      extract(models[[1]], "model")$SPHMS$TS,
      type = "l",
      lwd = 2.5,
      col = colors[1],
      xlab = "Frequency (kHz)",
      ylab = impl_ts_label(),
      ylim = c(-140, -30)
    )
    for (i in 2:length(models)) {
      graphics::lines(
        extract(models[[i]], "model")$SPHMS$frequency * 1e-3,
        extract(models[[i]], "model")$SPHMS$TS,
        col = colors[i],
        lwd = 2
      )
    }
    graphics::legend(
      "bottomright",
      legend = names(models),
      col = colors,
      lwd = c(2.5, 2, 2),
      bty = "n"
    )
  }
)
