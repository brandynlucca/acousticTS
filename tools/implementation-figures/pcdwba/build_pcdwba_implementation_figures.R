source("tools/implementation-figures/helpers/common.R")
impl_load_all()

compare_df <- utils::read.csv(impl_data_path("pcdwba_reference_compare.csv"))

impl_with_png(
  impl_output_path("pcdwba", "pcdwba-reference-comparison.png"),
  {
    old_par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old_par), add = TRUE)
    graphics::par(mfrow = c(2, 1), mar = c(4, 4.5, 2, 1))
    graphics::plot(
      compare_df$frequency_khz,
      compare_df$TS_acousticts,
      type = "l",
      lwd = 2.5,
      col = "black",
      xlab = "",
      ylab = impl_ts_label()
    )
    graphics::lines(
      compare_df$frequency_khz,
      compare_df$TS_zooscatr,
      col = "firebrick3",
      lwd = 2
    )
    graphics::lines(
      compare_df$frequency_khz,
      compare_df$TS_echopop,
      col = "dodgerblue3",
      lwd = 2,
      lty = 2
    )
    graphics::legend(
      "bottomleft",
      legend = c("acousticTS", "ZooScatR", "echopop"),
      col = c("black", "firebrick3", "dodgerblue3"),
      lwd = c(2.5, 2, 2),
      lty = c(1, 1, 2),
      bty = "n"
    )

    graphics::plot(
      compare_df$frequency_khz,
      compare_df$TS_acousticts - compare_df$TS_zooscatr,
      type = "l",
      lwd = 2,
      col = "firebrick3",
      xlab = "Frequency (kHz)",
      ylab = impl_delta_label(),
      ylim = c(-0.07, 0.07)
    )
    graphics::lines(
      compare_df$frequency_khz,
      compare_df$TS_acousticts - compare_df$TS_echopop,
      col = "dodgerblue3",
      lwd = 2,
      lty = 2
    )
    graphics::abline(h = 0, lty = 3, col = "grey50")
  }
)
