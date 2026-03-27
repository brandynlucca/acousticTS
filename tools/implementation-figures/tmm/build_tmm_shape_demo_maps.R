devtools::load_all(getwd())

density_sw <- 1026.8
sound_speed_sw <- 1477.3

render_scattering_pair <- function(object,
                                   heatmap_file,
                                   polar_file,
                                   frequency,
                                   n_theta = 61,
                                   n_phi = 121,
                                   width = 1800,
                                   height = 1600,
                                   res = 200) {
  png(filename = heatmap_file, width = width, height = height, res = res)
  plot(
    object,
    type = "scattering",
    frequency = frequency,
    heatmap = TRUE,
    n_theta = n_theta,
    n_phi = n_phi
  )
  dev.off()

  png(filename = polar_file, width = width, height = height, res = res)
  plot(
    object,
    type = "scattering",
    frequency = frequency,
    polar = TRUE,
    n_theta = n_theta,
    n_phi = n_phi
  )
  dev.off()
}

sphere_object <- target_strength(
  object = fls_generate(
    shape = sphere(radius_body = 0.01, n_segments = 80),
    g_body = 1,
    h_body = 1,
    theta_body = pi / 2
  ),
  frequency = 70e3,
  model = "tmm",
  boundary = "pressure_release",
  density_sw = density_sw,
  sound_speed_sw = sound_speed_sw,
  store_t_matrix = TRUE
)

oblate_object <- target_strength(
  object = fls_generate(
    shape = oblate_spheroid(length_body = 0.012, radius_body = 0.01, n_segments = 80),
    density_body = 1028.9,
    sound_speed_body = 1480.3,
    theta_body = pi / 2
  ),
  frequency = 70e3,
  model = "tmm",
  boundary = "liquid_filled",
  density_sw = density_sw,
  sound_speed_sw = sound_speed_sw,
  store_t_matrix = TRUE
)

prolate_object <- target_strength(
  object = fls_generate(
    shape = prolate_spheroid(length_body = 0.14, radius_body = 0.01, n_segments = 80),
    density_body = 1028.9,
    sound_speed_body = 1480.3,
    theta_body = pi / 2
  ),
  frequency = 70e3,
  model = "tmm",
  boundary = "liquid_filled",
  density_sw = density_sw,
  sound_speed_sw = sound_speed_sw,
  store_t_matrix = TRUE
)

render_scattering_pair(
  sphere_object,
  heatmap_file = "vignettes/tmm/tmm-sphere-heatmap.png",
  polar_file = "vignettes/tmm/tmm-sphere-polar.png",
  frequency = 70e3
)

render_scattering_pair(
  oblate_object,
  heatmap_file = "vignettes/tmm/tmm-oblate-heatmap.png",
  polar_file = "vignettes/tmm/tmm-oblate-polar.png",
  frequency = 70e3
)

render_scattering_pair(
  prolate_object,
  heatmap_file = "vignettes/tmm/tmm-prolate-heatmap.png",
  polar_file = "vignettes/tmm/tmm-prolate-polar.png",
  frequency = 70e3
)
