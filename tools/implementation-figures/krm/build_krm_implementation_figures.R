source("tools/implementation-figures/helpers/common.R")
impl_load_all()

body_shape <- arbitrary(
  x_body = c(0.00, 0.04, 0.09, 0.15, 0.22, 0.28),
  w_body = c(0.00, 0.018, 0.028, 0.030, 0.018, 0.00),
  zU_body = c(0.000, 0.010, 0.016, 0.017, 0.010, 0.000),
  zL_body = c(0.000, -0.010, -0.016, -0.017, -0.010, 0.000)
)

bladder_shape <- arbitrary(
  x_bladder = c(0.06, 0.10, 0.14, 0.18, 0.21),
  w_bladder = c(0.00, 0.008, 0.012, 0.008, 0.00),
  zU_bladder = c(0.000, 0.004, 0.006, 0.004, 0.000),
  zL_bladder = c(0.000, -0.004, -0.006, -0.004, 0.000)
)

fish_object <- sbf_generate(
  body_shape = body_shape,
  bladder_shape = bladder_shape,
  density_body = 1070,
  sound_speed_body = 1570,
  density_bladder = 1.2,
  sound_speed_bladder = 340,
  theta_body = pi / 2,
  theta_bladder = pi / 2
)

fish_object <- target_strength(
  object = fish_object,
  frequency = seq(10e3, 400e3, by = 2e3),
  model = "krm"
)

body_only_object <- fls_generate(
  shape = body_shape,
  density_body = 1070,
  sound_speed_body = 1570,
  theta_body = pi / 2
)

body_only_object <- target_strength(
  object = body_only_object,
  frequency = seq(10e3, 400e3, by = 2e3),
  model = "krm"
)

data(sardine, package = "acousticTS")
variant_frequency <- seq(0.5e3, 100e3, by = 0.5e3)

variant_models <- lapply(
  c("lowcontrast", "body_embedded", "mixed"),
  function(krm_variant) {
    target_strength(
      object = sardine,
      frequency = variant_frequency,
      model = "krm",
      krm_variant = krm_variant
    )
  }
)

names(variant_models) <- c("lowcontrast", "body_embedded", "mixed")

impl_with_png(
  impl_output_path("krm", "krm-shape-plot.png"),
  {
    plot(fish_object, type = "shape")
  },
  width = 1400,
  height = 1100
)

impl_with_png(
  impl_output_path("krm", "krm-model-plot.png"),
  {
    plot(fish_object, type = "model")
  },
  width = 1100,
  height = 1100
)

impl_with_png(
  impl_output_path("krm", "krm-body-vs-swimbladder.png"),
  {
    fish_df <- extract(fish_object, "model")$KRM
    body_df <- extract(body_only_object, "model")$KRM
    graphics::plot(
      fish_df$frequency * 1e-3,
      body_df$TS,
      type = "l",
      lwd = 2.5,
      col = "black",
      xlab = "Frequency (kHz)",
      ylab = impl_ts_label(),
      ylim = c(-70, -20)
    )
    graphics::lines(
      fish_df$frequency * 1e-3,
      fish_df$TS,
      col = "firebrick3",
      lwd = 2
    )
    graphics::legend(
      "bottomleft",
      legend = c("Body only", "Body + swimbladder"),
      col = c("black", "firebrick3"),
      lwd = c(2.5, 2),
      bty = "n"
    )
  }
)

impl_with_png(
  impl_output_path("krm", "krm-variant-comparison.png"),
  {
    colors <- c("black", "firebrick3", "dodgerblue3")
    graphics::plot(
      extract(variant_models[[1]], "model")$KRM$frequency * 1e-3,
      extract(variant_models[[1]], "model")$KRM$TS,
      type = "l",
      lwd = 2.5,
      col = colors[1],
      xlab = "Frequency (kHz)",
      ylab = impl_ts_label(),
      ylim = c(-60, -30)
    )
    for (i in 2:length(variant_models)) {
      graphics::lines(
        extract(variant_models[[i]], "model")$KRM$frequency * 1e-3,
        extract(variant_models[[i]], "model")$KRM$TS,
        col = colors[i],
        lwd = 2
      )
    }
    graphics::legend(
      "bottomleft",
      legend = names(variant_models),
      col = colors,
      lwd = c(2.5, 2, 2),
      bty = "n"
    )
  }
)
