source("tools/implementation-figures/helpers/common.R")
impl_load_all()

make_sdwba_object <- function(phase_sd_init = sqrt(2) / 2) {
  obj <- fls_generate(
    shape = cylinder(
      length_body = 15e-3,
      radius_body = 2e-3,
      n_segments = 50
    ),
    g_body = 1.058,
    h_body = 1.058,
    theta_body = pi / 2
  )

  target_strength(
    object = obj,
    frequency = seq(10e3, 200e3, by = 2e3),
    model = "sdwba",
    n_iterations = 100,
    n_segments_init = 14,
    phase_sd_init = phase_sd_init,
    length_init = 15e-3,
    frequency_init = 120e3
  )
}

baseline_object <- make_sdwba_object()

impl_with_png(
  impl_output_path("sdwba", "sdwba-model-plot.png"),
  {
    plot(baseline_object, type = "model")
  },
  width = 1100,
  height = 1100
)

phase_models <- list(
  `Lower phase SD` = make_sdwba_object(phase_sd_init = 0.25),
  `Higher phase SD` = make_sdwba_object(phase_sd_init = 0.75)
)

impl_with_png(
  impl_output_path("sdwba", "sdwba-phase-comparison.png"),
  {
    colors <- c("black", "firebrick3")
    graphics::plot(
      extract(phase_models[[1]], "model")$SDWBA$frequency * 1e-3,
      extract(phase_models[[1]], "model")$SDWBA$TS,
      type = "l",
      lwd = 2.5,
      col = colors[1],
      xlab = "Frequency (kHz)",
      ylab = impl_ts_label()
    )
    graphics::lines(
      extract(phase_models[[2]], "model")$SDWBA$frequency * 1e-3,
      extract(phase_models[[2]], "model")$SDWBA$TS,
      col = colors[2],
      lwd = 2
    )
    graphics::legend(
      "bottomright",
      legend = names(phase_models),
      col = colors,
      lwd = c(2.5, 2),
      bty = "n"
    )
  }
)
