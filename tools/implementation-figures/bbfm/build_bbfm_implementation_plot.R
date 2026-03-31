devtools::load_all(".")

density_sw <- 1026.8
sound_speed_sw <- 1477.3

body_shape <- cylinder(
  length_body = 0.08,
  radius_body = 0.003,
  n_segments = 40
)

backbone_shape <- cylinder(
  length_body = 0.05,
  radius_body = 0.0006,
  n_segments = 20
)

bbf_object <- bbf_generate(
  body_shape = body_shape,
  backbone_shape = backbone_shape,
  density_body = 1070,
  sound_speed_body = 1570,
  density_backbone = 1900,
  sound_speed_longitudinal_backbone = 3500,
  sound_speed_transversal_backbone = 1700,
  x_offset_backbone = 0.015
)

bbf_object <- target_strength(
  object = bbf_object,
  frequency = seq(10e3, 400e3, by = 1e3),
  model = "bbfm",
  density_sw = density_sw,
  sound_speed_sw = sound_speed_sw
)

bbfm_out <- extract(bbf_object, "model")$BBFM

body_object <- methods::new(
  "FLS",
  metadata = list(ID = "body"),
  model_parameters = list(),
  model = list(),
  body = extract(bbf_object, "body"),
  shape_parameters = extract(bbf_object, c("shape_parameters", "body"))
)

body_object <- target_strength(
  object = body_object,
  frequency = bbfm_out$frequency,
  model = "dwba",
  density_sw = density_sw,
  sound_speed_sw = sound_speed_sw
)

backbone_object <- methods::new(
  "FLS",
  metadata = list(ID = "backbone"),
  model_parameters = list(),
  model = list(),
  body = extract(bbf_object, "backbone"),
  shape_parameters = extract(bbf_object, c("shape_parameters", "backbone"))
)

backbone_object <- target_strength(
  object = backbone_object,
  frequency = bbfm_out$frequency,
  model = "ecms",
  density_sw = density_sw,
  sound_speed_sw = sound_speed_sw,
  density_body = extract(bbf_object, c("backbone", "density")),
  sound_speed_longitudinal_body = extract(
    bbf_object,
    c("backbone", "sound_speed_longitudinal")
  ),
  sound_speed_transversal_body = extract(
    bbf_object,
    c("backbone", "sound_speed_transversal")
  )
)

backbone_body <- extract(backbone_object, "body")
x_center <- mean(range(backbone_body$rpos["x", ], na.rm = TRUE))
z_center <- mean(range(backbone_body$rpos["z", ], na.rm = TRUE))
phase_shift <- exp(
  2i * acousticTS::wavenumber(bbfm_out$frequency, sound_speed_sw) *
    (x_center * cos(backbone_body$theta) + z_center * sin(backbone_body$theta))
)

reconstructed_fbs <-
  extract(body_object, "model")$DWBA$f_bs +
  extract(backbone_object, "model")$ECMS$f_bs * phase_shift

bbfm_check <- data.frame(
  frequency_kHz = bbfm_out$frequency * 1e-3,
  delta_f_bs = Mod(bbfm_out$f_bs - reconstructed_fbs),
  delta_TS_dB = bbfm_out$TS - 10 * log10(abs(reconstructed_fbs)^2)
)

png("vignettes/bbfm/bbfm-example-spectrum.png", width = 1600, height = 1200, res = 200)
old_par <- par(no.readonly = TRUE)
par(mfrow = c(2, 1), mar = c(4, 4.5, 2, 1))

plot(
  bbfm_out$frequency * 1e-3,
  bbfm_out$TS,
  type = "l",
  lwd = 2.5,
  col = "black",
  xlab = "",
  ylab = expression(Target ~ strength ~ (dB ~ re. ~ 1 ~ m^2))
)
lines(bbfm_out$frequency * 1e-3, bbfm_out$TS_body, col = "#0b6e4f", lwd = 2)
lines(
  bbfm_out$frequency * 1e-3,
  bbfm_out$TS_backbone,
  col = "#b2472f",
  lwd = 2,
  lty = 2
)
legend(
  "bottomright",
  legend = c("Composite", "Flesh body", "Backbone"),
  col = c("black", "#0b6e4f", "#b2472f"),
  lwd = c(2.5, 2, 2),
  lty = c(1, 1, 2),
  bty = "n"
)

plot(
  bbfm_check$frequency_kHz,
  bbfm_check$delta_TS_dB,
  type = "l",
  lwd = 2,
  col = "#7a3e9d",
  xlab = "Frequency (kHz)",
  ylab = expression(Delta * TS ~ (dB))
)
abline(h = 0, lty = 3, col = "grey40")

par(old_par)
dev.off()
