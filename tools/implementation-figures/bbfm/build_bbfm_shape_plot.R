devtools::load_all('.')

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

png('vignettes/bbfm/bbfm-shape-plot.png', width = 1200, height = 1200, res = 200)
plot(bbf_object, type = 'shape')
dev.off()
