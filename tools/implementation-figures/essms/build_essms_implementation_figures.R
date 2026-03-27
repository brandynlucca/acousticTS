source("tools/implementation-figures/helpers/common.R")
impl_load_all()

make_essms_object <- function(shell_thickness = 0.8e-3) {
  ess_generate(
    shape = sphere(radius_body = 10e-3, n_segments = 80),
    radius_shell = 10e-3,
    shell_thickness = shell_thickness,
    density_shell = 2565,
    sound_speed_shell = 3750,
    density_fluid = 1077.3,
    sound_speed_fluid = 1575,
    E = 7.0e10,
    nu = 0.32
  )
}

baseline_object <- make_essms_object()
baseline_object <- target_strength(
  object = baseline_object,
  frequency = seq(10e3, 200e3, by = 1e3),
  model = "essms"
)

impl_with_png(
  impl_output_path("essms", "essms-shape-plot.png"),
  {
    plot(make_essms_object(), type = "shape")
  },
  width = 1100,
  height = 1100
)

thicknesses <- c(0.5e-3, 0.8e-3, 1.2e-3)
thickness_models <- lapply(
  thicknesses,
  function(shell_thickness) {
    target_strength(
      object = make_essms_object(shell_thickness = shell_thickness),
      frequency = seq(10e3, 200e3, by = 1e3),
      model = "essms"
    )
  }
)

impl_with_png(
  impl_output_path("essms", "essms-thickness-comparison.png"),
  {
    colors <- c("black", "firebrick3", "dodgerblue3")
    graphics::plot(
      extract(thickness_models[[1]], "model")$ESSMS$frequency * 1e-3,
      extract(thickness_models[[1]], "model")$ESSMS$TS,
      type = "l",
      lwd = 2.5,
      col = colors[1],
      xlab = "Frequency (kHz)",
      ylab = impl_ts_label(),
      ylim = c(-75, -25)
    )
    for (i in 2:length(thickness_models)) {
      graphics::lines(
        extract(thickness_models[[i]], "model")$ESSMS$frequency * 1e-3,
        extract(thickness_models[[i]], "model")$ESSMS$TS,
        col = colors[i],
        lwd = 2
      )
    }
    graphics::legend(
      "bottom",
      legend = c("0.5 mm", "0.8 mm", "1.2 mm"),
      col = colors,
      lwd = c(2.5, 2, 2),
      bty = "n"
    )
  }
)
