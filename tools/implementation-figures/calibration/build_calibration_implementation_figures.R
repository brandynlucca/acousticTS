source("tools/implementation-figures/helpers/common.R")
impl_load_all()

cal_sphere <- cal_generate()
frequency <- seq(1e3, 600e3, 1e3)
cal_sphere <- target_strength(
  object = cal_sphere,
  frequency = frequency,
  model = "calibration"
)

for (x_units in c("frequency", "k_sw", "k_l", "k_t")) {
  file_name <- switch(x_units,
    frequency = "calibration-spectrum-frequency.png",
    k_sw = "calibration-spectrum-k-sw.png",
    k_l = "calibration-spectrum-k-l.png",
    k_t = "calibration-spectrum-k-t.png"
  )
  impl_with_png(
    impl_output_path("calibration", file_name),
    {
      plot(cal_sphere, type = "model", x_units = x_units)
    },
    width = 1100,
    height = 1100
  )
}

diameters <- c(0.023, 0.0381, 0.06)
diameter_models <- lapply(
  diameters,
  function(diameter) {
    obj <- cal_generate(diameter = diameter)
    target_strength(
      object = obj,
      frequency = seq(1e3, 360e3, 1e3),
      model = "calibration"
    )
  }
)

impl_with_png(
  impl_output_path("calibration", "calibration-diameter-comparison.png"),
  {
    old_par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old_par), add = TRUE)
    graphics::par(mar = c(4, 4.5, 2, 1))
    graphics::plot(
      extract(diameter_models[[1]], "model")$calibration$frequency * 1e-3,
      extract(diameter_models[[1]], "model")$calibration$TS,
      type = "l",
      lwd = 2.5,
      col = "black",
      xlab = "Frequency (kHz)",
      ylab = impl_ts_label(),
      ylim = c(-95, -30)
    )
    graphics::lines(
      extract(diameter_models[[2]], "model")$calibration$frequency * 1e-3,
      extract(diameter_models[[2]], "model")$calibration$TS,
      col = "firebrick3",
      lwd = 2
    )
    graphics::lines(
      extract(diameter_models[[3]], "model")$calibration$frequency * 1e-3,
      extract(diameter_models[[3]], "model")$calibration$TS,
      col = "dodgerblue3",
      lwd = 2
    )
    graphics::legend(
      "bottomright",
      legend = c("23.0 mm", "38.1 mm", "60.0 mm"),
      col = c("black", "firebrick3", "dodgerblue3"),
      lwd = c(2.5, 2, 2),
      bty = "n"
    )
  }
)

materials <- c("WC", "steel", "Al")
material_models <- lapply(
  materials,
  function(material) {
    obj <- cal_generate(material = material, diameter = 38.1e-3)
    target_strength(
      object = obj,
      frequency = seq(1e3, 360e3, 1e3),
      model = "calibration"
    )
  }
)

impl_with_png(
  impl_output_path("calibration", "calibration-material-comparison.png"),
  {
    old_par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old_par), add = TRUE)
    graphics::par(mar = c(4, 4.5, 2, 1))
    graphics::plot(
      extract(material_models[[1]], "model")$calibration$frequency * 1e-3,
      extract(material_models[[1]], "model")$calibration$TS,
      type = "l",
      lwd = 2.5,
      col = "black",
      xlab = "Frequency (kHz)",
      ylab = impl_ts_label()
    )
    graphics::lines(
      extract(material_models[[2]], "model")$calibration$frequency * 1e-3,
      extract(material_models[[2]], "model")$calibration$TS,
      col = "firebrick3",
      lwd = 2
    )
    graphics::lines(
      extract(material_models[[3]], "model")$calibration$frequency * 1e-3,
      extract(material_models[[3]], "model")$calibration$TS,
      col = "dodgerblue3",
      lwd = 2
    )
    graphics::legend(
      "bottomright",
      legend = c("WC", "Steel", "Al"),
      col = c("black", "firebrick3", "dodgerblue3"),
      lwd = c(2.5, 2, 2),
      bty = "n"
    )
  }
)

compare_df <- impl_sort_compare_df(
  utils::read.csv(impl_data_path("calibration_wc381_package_compare.csv"))
)

impl_with_png(
  impl_output_path("calibration", "calibration-external-comparison.png"),
  {
    old_par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old_par), add = TRUE)
    graphics::par(mfrow = c(2, 1), mar = c(4, 4.5, 2, 1))
    graphics::plot(
      compare_df$frequency * 1e-3,
      compare_df$acousticTS,
      type = "l",
      lwd = 2.5,
      col = "black",
      xlab = "",
      ylab = impl_ts_label()
    )
    graphics::lines(
      compare_df$frequency * 1e-3,
      compare_df$echoSMs,
      col = "firebrick3",
      lwd = 2
    )
    graphics::lines(
      compare_df$frequency * 1e-3,
      compare_df$sphereTS,
      col = "dodgerblue3",
      lwd = 2,
      lty = 2
    )
    graphics::lines(
      compare_df$frequency * 1e-3,
      compare_df$noaa_applet,
      col = "darkgreen",
      lwd = 2,
      lty = 3
    )
    graphics::legend(
      "bottomright",
      legend = c("acousticTS", "echoSMs", "sphereTS", "NOAA"),
      col = c("black", "firebrick3", "dodgerblue3", "darkgreen"),
      lwd = c(2.5, 2, 2, 2),
      lty = c(1, 1, 2, 3),
      bty = "n"
    )

    graphics::plot(
      compare_df$frequency * 1e-3,
      compare_df$acousticTS - compare_df$echoSMs,
      type = "l",
      lwd = 2,
      col = "firebrick3",
      xlab = "Frequency (kHz)",
      ylab = impl_delta_label()
    )
    graphics::lines(
      compare_df$frequency * 1e-3,
      compare_df$acousticTS - compare_df$sphereTS,
      col = "dodgerblue3",
      lwd = 2,
      lty = 2
    )
    graphics::lines(
      compare_df$frequency * 1e-3,
      compare_df$acousticTS - compare_df$noaa_applet,
      col = "darkgreen",
      lwd = 2,
      lty = 3
    )
    graphics::abline(h = 0, lty = 3, col = "grey50")
  }
)
