source("tools/implementation-figures/helpers/common.R")
impl_load_all()

density_sw <- 1026.8
sound_speed_sw <- 1477.3

sphere_object <- fls_generate(
  shape = sphere(radius_body = 0.01),
  density_body = 1028.9,
  sound_speed_body = 1480.3,
  theta_body = pi / 2
)

sphere_object <- target_strength(
  object = sphere_object,
  frequency = seq(12e3, 120e3, by = 12e3),
  model = "tmm",
  boundary = "liquid_filled",
  density_sw = density_sw,
  sound_speed_sw = sound_speed_sw
)

impl_with_png(
  impl_output_path("tmm", "tmm-shape-plot.png"),
  {
    plot(sphere_object, type = "shape")
  },
  width = 1100,
  height = 1100
)

impl_with_png(
  impl_output_path("tmm", "tmm-model-plot.png"),
  {
    plot(sphere_object, type = "model")
  },
  width = 1100,
  height = 1100
)

prolate_store <- target_strength(
  object = fls_generate(
    shape = prolate_spheroid(
      length_body = 0.14,
      radius_body = 0.01,
      n_segments = 80
    ),
    density_body = 1028.9,
    sound_speed_body = 1480.3,
    theta_body = pi / 2
  ),
  frequency = seq(12e3, 100e3, 2e3),
  model = "tmm",
  boundary = "liquid_filled",
  density_sw = density_sw,
  sound_speed_sw = sound_speed_sw,
  store_t_matrix = TRUE
)

orientation_dist <- tmm_orientation_distribution(
  distribution = "uniform",
  lower = 0.45 * pi,
  upper = pi,
  n_theta = 7
)

orientation_avg <- tmm_average_orientation(
  object = prolate_store,
  distribution = orientation_dist
)

stored_df <- extract(prolate_store, "model")$TMM

impl_with_png(
  impl_output_path("tmm", "tmm-orientation-average-spectrum.png"),
  {
    graphics::plot(
      stored_df$frequency * 1e-3,
      stored_df$TS,
      type = "b",
      lwd = 2.5,
      pch = 16,
      col = "black",
      xlab = "Frequency (kHz)",
      ylab = impl_ts_label(),
      ylim = c(-105, -73)
    )
    graphics::lines(
      orientation_avg$frequency * 1e-3,
      orientation_avg$TS,
      type = "b",
      pch = 1,
      col = "firebrick3",
      lwd = 2,
      lty = 2
    )
    graphics::legend(
      "topright",
      legend = c("Stored monostatic", "Orientation-averaged"),
      col = c("black", "firebrick3"),
      lwd = c(2.5, 2),
      pch = c(16, 1),
      lty = c(1, 2),
      bty = "n"
    )
  }
)

impl_with_png(
  impl_output_path("tmm", "tmm-scattering-slice.png"),
  {
    plot(
      prolate_store,
      type = "scattering",
      frequency = 70e3,
      vary = "theta_scatter"
    )
  },
  width = 1400,
  height = 1100
)

sphere_frequency <- seq(12e3, 120e3, by = 4e3)
prolate_frequency <- seq(12e3, 120e3, by = 4e3)

sphere_plot_object <- fls_generate(
  shape = sphere(radius_body = 0.01, n_segments = 80),
  density_body = 1028.9,
  sound_speed_body = 1480.3,
  theta_body = pi / 2
)

sphere_liq <- list(
  tmm = target_strength(
    object = sphere_plot_object,
    frequency = sphere_frequency,
    model = "tmm",
    boundary = "liquid_filled",
    density_sw = density_sw,
    sound_speed_sw = sound_speed_sw
  ),
  reference = target_strength(
    object = sphere_plot_object,
    frequency = sphere_frequency,
    model = "sphms",
    boundary = "liquid_filled",
    density_sw = density_sw,
    sound_speed_sw = sound_speed_sw
  )
)

prolate_plot_object <- fls_generate(
  shape = prolate_spheroid(
    length_body = 0.14,
    radius_body = 0.01,
    n_segments = 80
  ),
  density_body = 1028.9,
  sound_speed_body = 1480.3,
  theta_body = pi / 2
)

prolate_liq <- list(
  tmm = target_strength(
    object = prolate_plot_object,
    frequency = prolate_frequency,
    model = "tmm",
    boundary = "liquid_filled",
    density_sw = density_sw,
    sound_speed_sw = sound_speed_sw
  ),
  reference = target_strength(
    object = prolate_plot_object,
    frequency = prolate_frequency,
    model = "psms",
    boundary = "liquid_filled",
    density_sw = density_sw,
    sound_speed_sw = sound_speed_sw,
    precision = "quad",
    simplify_Amn = FALSE
  )
)

impl_with_png(
  impl_output_path("tmm", "tmm-representative-spectra.png"),
  {
    old_par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old_par), add = TRUE)
    graphics::par(mfrow = c(2, 2), mar = c(4, 4.5, 2, 1))
    graphics::plot(
      sphere_frequency * 1e-3,
      sphere_liq$tmm@model$TMM$TS,
      type = "l",
      lwd = 2.5,
      xlab = "Frequency (kHz)",
      ylab = impl_ts_label(),
      main = "Sphere"
    )
    graphics::lines(
      sphere_frequency * 1e-3,
      sphere_liq$reference@model$SPHMS$TS,
      col = "firebrick3",
      lwd = 2,
      lty = 2
    )
    graphics::plot(
      sphere_frequency * 1e-3,
      sphere_liq$tmm@model$TMM$TS - sphere_liq$reference@model$SPHMS$TS,
      type = "l",
      lwd = 2,
      col = "firebrick3",
      xlab = "Frequency (kHz)",
      ylab = impl_delta_label(),
      main = "Sphere delta"
    )
    graphics::abline(h = 0, lty = 3, col = "grey50")
    graphics::plot(
      prolate_frequency * 1e-3,
      prolate_liq$tmm@model$TMM$TS,
      type = "b",
      lwd = 2.5,
      pch = 16,
      xlab = "Frequency (kHz)",
      ylab = impl_ts_label(),
      main = "Prolate"
    )
    graphics::lines(
      prolate_frequency * 1e-3,
      prolate_liq$reference@model$PSMS$TS,
      col = "dodgerblue3",
      lwd = 2,
      lty = 2,
      type = "b",
      pch = 1
    )
    graphics::plot(
      prolate_frequency * 1e-3,
      prolate_liq$tmm@model$TMM$TS -
        prolate_liq$reference@model$PSMS$TS,
      type = "b",
      lwd = 2,
      pch = 16,
      col = "firebrick3",
      xlab = "Frequency (kHz)",
      ylab = impl_delta_label(),
      main = "Prolate delta"
    )
    graphics::abline(h = 0, lty = 3, col = "grey50")
  },
  width = 1800,
  height = 1600
)
