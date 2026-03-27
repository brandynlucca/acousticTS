source("tools/implementation-figures/helpers/common.R")
impl_load_all()

make_trcm_object <- function(radius_curvature_ratio = NULL) {
  obj <- fls_generate(
    shape = cylinder(
      length_body = 20e-3,
      radius_body = 1.5e-3,
      radius_curvature_ratio = radius_curvature_ratio,
      n_segments = 60
    ),
    density_body = 1045,
    sound_speed_body = 1520,
    theta_body = pi / 2
  )

  target_strength(
    object = obj,
    frequency = seq(10e3, 400e3, by = 2e3),
    model = "trcm"
  )
}

straight_object <- make_trcm_object()
bent_object <- make_trcm_object(radius_curvature_ratio = 1.5)

impl_with_png(
  impl_output_path("trcm", "trcm-shape-plot.png"),
  {
    plot(straight_object, type = "shape")
  },
  width = 1100,
  height = 1100
)

impl_with_png(
  impl_output_path("trcm", "trcm-model-plot.png"),
  {
    plot(straight_object, type = "model")
  },
  width = 1100,
  height = 1100
)

impl_with_png(
  impl_output_path("trcm", "trcm-straight-vs-bent.png"),
  {
    straight_df <- extract(straight_object, "model")$TRCM
    bent_df <- extract(bent_object, "model")$TRCM
    graphics::plot(
      straight_df$frequency * 1e-3,
      straight_df$TS,
      type = "l",
      lwd = 2.5,
      col = "black",
      xlab = "Frequency (kHz)",
      ylab = impl_ts_label()
    )
    graphics::lines(
      bent_df$frequency * 1e-3,
      bent_df$TS,
      col = "firebrick3",
      lwd = 2
    )
    graphics::legend(
      "topleft",
      legend = c("Straight", "Bent"),
      col = c("black", "firebrick3"),
      lwd = c(2.5, 2),
      bty = "n"
    )
  }
)
