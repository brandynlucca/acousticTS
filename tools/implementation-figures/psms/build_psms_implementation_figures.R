source("tools/implementation-figures/helpers/common.R")
impl_load_all()

make_psms_object <- function(boundary = "liquid_filled",
                             precision = "double",
                             simplify_Amn = TRUE) {
  obj <- fls_generate(
    shape = prolate_spheroid(
      length_body = 40e-3,
      radius_body = 4e-3,
      n_segments = 100
    ),
    density_body = 1045,
    sound_speed_body = 1520,
    theta_body = pi / 2
  )

  target_strength(
    object = obj,
    frequency = seq(10e3, 120e3, by = 2e3),
    model = "psms",
    boundary = boundary,
    phi_body = pi,
    n_integration = 96,
    simplify_Amn = simplify_Amn,
    precision = precision
  )
}

example_object <- make_psms_object()
rigid_object <- make_psms_object(boundary = "fixed_rigid")

impl_with_png(
  impl_output_path("psms", "psms-shape-plot.png"),
  {
    plot(example_object, type = "shape")
  },
  width = 1100,
  height = 1100
)

impl_with_png(
  impl_output_path("psms", "psms-boundary-comparison.png"),
  {
    liquid_df <- extract(example_object, "model")$PSMS
    rigid_df <- extract(rigid_object, "model")$PSMS
    graphics::plot(
      liquid_df$frequency * 1e-3,
      liquid_df$TS,
      type = "l",
      lwd = 2.5,
      col = "black",
      xlab = "Frequency (kHz)",
      ylab = impl_ts_label(),
      ylim = c(-100, -35)
    )
    graphics::lines(
      rigid_df$frequency * 1e-3,
      rigid_df$TS,
      col = "firebrick3",
      lwd = 2
    )
    graphics::legend(
      "bottomright",
      legend = c("Liquid-filled", "Fixed rigid"),
      col = c("black", "firebrick3"),
      lwd = c(2.5, 2),
      bty = "n"
    )
  }
)

# precision_freq <- c(12e3, 18e3, 38e3, 70e3, 100e3)
precision_freq <- seq(10e3, 100e3, 2e3)
precision_object <- fls_generate(
  shape = prolate_spheroid(
    length_body = 0.14,
    radius_body = 0.01,
    n_segments = 80
  ),
  theta_body = pi / 2,
  density_body = 1028.9,
  sound_speed_body = 1480.3
)

precision_double <- target_strength(
  object = precision_object,
  frequency = precision_freq,
  model = "psms",
  boundary = "liquid_filled",
  phi_body = pi,
  adaptive = FALSE,
  precision = "double",
  simplify_Amn = FALSE,
  n_integration = 96,
  density_sw = 1026.8,
  sound_speed_sw = 1477.3
)

precision_quad <- target_strength(
  object = precision_object,
  frequency = precision_freq,
  model = "psms",
  boundary = "liquid_filled",
  phi_body = pi,
  adaptive = FALSE,
  precision = "quad",
  simplify_Amn = FALSE,
  n_integration = 96,
  density_sw = 1026.8,
  sound_speed_sw = 1477.3
)

precision_df <- data.frame(
  frequency = precision_freq,
  k0b = 2 * pi * precision_freq * 0.01 / 1477.3,
  double_ts = precision_double@model$PSMS$TS,
  quad_ts = precision_quad@model$PSMS$TS
)
precision_df$delta_ts <- precision_df$double_ts - precision_df$quad_ts

impl_with_png(
  impl_output_path("psms", "psms-precision-drift.png"),
  {
    old_par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old_par), add = TRUE)
    graphics::par(mfrow = c(2, 1), mar = c(4, 4.5, 2, 1))
    graphics::plot(
      precision_df$k0b,
      precision_df$double_ts,
      type = "b",
      lwd = 2.5,
      pch = 16,
      col = "black",
      xlab = expression(k[0] * b),
      ylab = impl_ts_label()
    )
    graphics::lines(
      precision_df$k0b,
      precision_df$quad_ts,
      type = "b",
      pch = 1,
      col = "firebrick3",
      lwd = 2,
      lty = 2
    )
    graphics::legend(
      "bottomleft",
      legend = c("Double", "Quad"),
      col = c("black", "firebrick3"),
      lwd = c(2.5, 2),
      pch = c(16, 1),
      lty = c(1, 2),
      bty = "n"
    )
    graphics::plot(
      precision_df$k0b,
      precision_df$delta_ts,
      type = "b",
      lwd = 2,
      pch = 16,
      col = "firebrick3",
      xlab = expression(k[0] * b),
      ylab = impl_delta_label()
    )
    graphics::abline(h = 0, lty = 3, col = "grey50")
  },
  width = 1400,
  height = 1400
)
