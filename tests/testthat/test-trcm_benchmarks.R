test_that("TRCM remains close to straight and bent cylinder reference problems", {
  density_sw <- 1026.8
  sound_speed_sw <- 1477.3
  g_body <- 1.0357
  h_body <- 1.0279
  length_body <- 10.5e-3
  radius_body <- 1e-3
  radius_curvature_ratio <- 1.5
  ka <- seq(0.5, 5, length.out = 12)
  frequency <- ka * sound_speed_sw / (2 * pi * radius_body)

  straight <- fls_generate(
    shape = cylinder(
      length_body = length_body,
      radius_body = radius_body,
      n_segments = 401
    ),
    density_body = density_sw * g_body,
    sound_speed_body = sound_speed_sw * h_body,
    theta_body = pi / 2
  )

  bent <- fls_generate(
    shape = cylinder(
      length_body = length_body,
      radius_body = radius_body,
      radius_curvature_ratio = radius_curvature_ratio,
      n_segments = 401
    ),
    density_body = density_sw * g_body,
    sound_speed_body = sound_speed_sw * h_body,
    theta_body = pi / 2
  )

  straight_fcms <- target_strength(
    straight,
    frequency,
    "FCMS",
    boundary = "liquid_filled",
    density_sw = density_sw,
    sound_speed_sw = sound_speed_sw
  )

  straight_trcm <- target_strength(
    straight,
    frequency,
    "TRCM",
    density_sw = density_sw,
    sound_speed_sw = sound_speed_sw
  )

  bent_trcm_fresnel <- target_strength(
    bent,
    frequency,
    "TRCM",
    density_sw = density_sw,
    sound_speed_sw = sound_speed_sw,
    stationary_phase = FALSE
  )

  bent_trcm_stationary <- target_strength(
    bent,
    frequency,
    "TRCM",
    density_sw = density_sw,
    sound_speed_sw = sound_speed_sw,
    stationary_phase = TRUE
  )

  lebc <- .trcm_equivalent_length_fresnel(
    k1 = 2 * pi * frequency / sound_speed_sw,
    l = length_body,
    a = radius_body,
    rho_c = radius_curvature_ratio * length_body
  )

  bent_reference_ts <- 20 * log10(abs(lebc * straight_fcms@model$FCMS$f_bs / length_body))

  expect_lt(
    mean(abs(straight_trcm@model$TRCM$TS - straight_fcms@model$FCMS$TS)),
    1
  )
  expect_lt(
    mean(abs(bent_trcm_fresnel@model$TRCM$TS - bent_reference_ts)),
    1
  )
  expect_lt(
    mean(abs(bent_trcm_stationary@model$TRCM$TS - bent_reference_ts)),
    2
  )
})
