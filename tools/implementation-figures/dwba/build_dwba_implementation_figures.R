source("tools/implementation-figures/helpers/common.R")
impl_load_all()

make_dwba_object <- function(theta_body = pi / 2) {
  fls_generate(
    shape = cylinder(
      length_body = 15e-3,
      radius_body = 2e-3,
      n_segments = 50
    ),
    g_body = 1.03,
    h_body = 1.03,
    theta_body = theta_body
  )
}

frequency <- seq(10e3, 400e3, by = 2e3)
example_object <- target_strength(
  object = make_dwba_object(),
  frequency = frequency,
  model = "dwba"
)

impl_with_png(
  impl_output_path("dwba", "dwba-shape-plot.png"),
  {
    plot(make_dwba_object(), type = "shape")
  },
  width = 1100,
  height = 1100
)

impl_with_png(
  impl_output_path("dwba", "dwba-model-plot.png"),
  {
    plot(example_object, type = "model")
  },
  width = 1100,
  height = 1100
)

orientations <- list(
  Broadside = pi / 2,
  Oblique = pi / 4,
  `End-on` = 0
)

orientation_models <- lapply(
  orientations,
  function(theta_body) {
    target_strength(
      object = make_dwba_object(theta_body = theta_body),
      frequency = frequency,
      model = "dwba"
    )
  }
)

impl_with_png(
  impl_output_path("dwba", "dwba-orientation-comparison.png"),
  {
    colors <- c("black", "firebrick3", "dodgerblue3")
    graphics::plot(
      extract(orientation_models[[1]], "model")$DWBA$frequency * 1e-3,
      extract(orientation_models[[1]], "model")$DWBA$TS,
      type = "l",
      lwd = 2.5,
      col = colors[1],
      xlab = "Frequency (kHz)",
      ylab = impl_ts_label(),
      ylim = c(-140, -65)
    )
    for (i in 2:length(orientation_models)) {
      graphics::lines(
        extract(orientation_models[[i]], "model")$DWBA$frequency * 1e-3,
        extract(orientation_models[[i]], "model")$DWBA$TS,
        col = colors[i],
        lwd = 2
      )
    }
    graphics::legend(
      "bottomleft",
      legend = names(orientation_models),
      col = colors,
      lwd = c(2.5, 2, 2),
      bty = "n"
    )
  }
)
