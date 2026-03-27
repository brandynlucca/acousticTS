source("tools/implementation-figures/helpers/common.R")
impl_load_all()

density_sw <- 1026.8
sound_speed_sw <- 1477.3

straight_shape <- cylinder(
  length_body = 10.5e-3,
  radius_body = 1e-3,
  n_segments = 401
)

straight_object <- fls_generate(
  shape = straight_shape,
  density_body = density_sw * 1.0357,
  sound_speed_body = sound_speed_sw * 1.0279,
  theta_body = pi / 2
)

bent_object <- brake(
  straight_object,
  radius_curvature = 1.5
)

frequency <- seq(12e3, 400e3, by = 2e3)

straight_object <- target_strength(
  object = straight_object,
  frequency = frequency,
  model = "bcms",
  density_sw = density_sw,
  sound_speed_sw = sound_speed_sw
)

bent_object <- target_strength(
  object = bent_object,
  frequency = frequency,
  model = "bcms",
  density_sw = density_sw,
  sound_speed_sw = sound_speed_sw
)

impl_with_png(
  impl_output_path("bcms", "bcms-shape-plot.png"),
  {
    old_par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old_par), add = TRUE)
    graphics::par(mfrow = c(1, 2))
    plot(straight_object, type = "shape")
    graphics::title("Straight cylinder")
    plot(bent_object, type = "shape")
    graphics::title("Bent cylinder")
  },
  width = 1800,
  height = 900
)

impl_with_png(
  impl_output_path("bcms", "bcms-model-plot.png"),
  {
    straight_ts <- extract(straight_object, "model")$BCMS
    bent_ts <- extract(bent_object, "model")$BCMS
    graphics::plot(
      straight_ts$frequency * 1e-3,
      straight_ts$TS,
      type = "l",
      lwd = 2.5,
      col = "black",
      xlab = "Frequency (kHz)",
      ylab = impl_ts_label()
    )
    graphics::lines(
      bent_ts$frequency * 1e-3,
      bent_ts$TS,
      col = "firebrick3",
      lwd = 2
    )
    graphics::legend(
      "bottomright",
      legend = c("Straight", "Bent"),
      col = c("black", "firebrick3"),
      lwd = c(2.5, 2),
      bty = "n"
    )
  }
)

compare_df <- impl_sort_compare_df(
  utils::read.csv(impl_data_path("bcms_reference_compare.csv"))
)

impl_with_png(
  impl_output_path("bcms", "bcms-reference-comparison.png"),
  {
    old_par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old_par), add = TRUE)
    graphics::par(mfrow = c(2, 2), mar = c(4, 4.5, 2, 1))
    for (case_name in c("straight", "bent")) {
      case_df <- compare_df[compare_df$case == case_name, ]
      case_df <- case_df[order(case_df$frequency_hz), , drop = FALSE]
      graphics::plot(
        case_df$frequency_khz,
        case_df$TS_reference,
        type = "l",
        lwd = 2.5,
        col = "black",
        xlab = "",
        ylab = impl_ts_label(),
        main = tools::toTitleCase(case_name)
      )
      graphics::lines(
        case_df$frequency_khz,
        case_df$TS_acousticts,
        col = "firebrick3",
        lwd = 2
      )
      graphics::plot(
        case_df$frequency_khz,
        case_df$TS_acousticts - case_df$TS_reference,
        type = "l",
        lwd = 2,
        col = "firebrick3",
        xlab = "Frequency (kHz)",
        ylab = impl_delta_label()
      )
      graphics::abline(h = 0, lty = 3, col = "grey50")
    }
  },
  width = 1800,
  height = 1400
)
