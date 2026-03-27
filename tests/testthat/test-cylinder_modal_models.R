library(acousticTS)

test_that("BCMS reproduces the exact straight-cylinder kernel and bent correction", {
  density_sw <- 1026.8
  sound_speed_sw <- 1477.3
  density_body <- density_sw * 1.0357
  sound_speed_body <- sound_speed_sw * 1.0279
  length_body <- 10.5e-3
  radius_body <- 1e-3
  radius_curvature_ratio <- 1.5
  frequency <- seq(80e3, 200e3, length.out = 5)

  straight <- fls_generate(
    shape = cylinder(
      length_body = length_body,
      radius_body = radius_body,
      n_segments = 401
    ),
    density_body = density_body,
    sound_speed_body = sound_speed_body,
    theta_body = pi / 2
  )

  bent <- fls_generate(
    shape = cylinder(
      length_body = length_body,
      radius_body = radius_body,
      radius_curvature_ratio = radius_curvature_ratio,
      n_segments = 401
    ),
    density_body = density_body,
    sound_speed_body = sound_speed_body,
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

  straight_bcms <- target_strength(
    straight,
    frequency,
    "BCMS",
    boundary = "liquid_filled",
    density_sw = density_sw,
    sound_speed_sw = sound_speed_sw
  )

  bent_bcms <- target_strength(
    bent,
    frequency,
    "BCMS",
    boundary = "liquid_filled",
    density_sw = density_sw,
    sound_speed_sw = sound_speed_sw,
    stationary_phase = FALSE
  )

  lebc <- .bcms_equivalent_length_fresnel(
    k1 = 2 * pi * frequency / sound_speed_sw,
    l = length_body,
    a = radius_body,
    rho_c = radius_curvature_ratio * length_body
  )

  bent_reference_ts <- 20 * log10(
    abs(lebc * straight_fcms@model$FCMS$f_bs / length_body)
  )

  expect_equal(
    straight_bcms@model$BCMS$TS,
    straight_fcms@model$FCMS$TS,
    tolerance = 1e-10
  )
  expect_equal(
    bent_bcms@model$BCMS$TS,
    bent_reference_ts,
    tolerance = 1e-10
  )
})

test_that("BCMS infers default boundaries and guards limited options", {
  density_sw <- 1026.8
  sound_speed_sw <- 1477.3
  frequency <- c(38e3, 70e3)

  liquid_object <- fls_generate(
    shape = cylinder(
      length_body = 0.03,
      radius_body = 0.003,
      n_segments = 101
    ),
    density_body = density_sw * 1.03,
    sound_speed_body = sound_speed_sw * 1.02,
    theta_body = pi / 2
  )
  gas_object <- gas_generate(
    shape = cylinder(
      length_body = 0.03,
      radius_body = 0.003,
      n_segments = 101
    ),
    density_fluid = 1.24,
    sound_speed_fluid = 345,
    theta_body = pi / 2
  )

  liquid_default <- target_strength(
    liquid_object,
    frequency = frequency,
    model = "bcms",
    density_sw = density_sw,
    sound_speed_sw = sound_speed_sw
  )
  liquid_explicit <- target_strength(
    liquid_object,
    frequency = frequency,
    model = "bcms",
    boundary = "liquid_filled",
    density_sw = density_sw,
    sound_speed_sw = sound_speed_sw
  )
  gas_default <- target_strength(
    gas_object,
    frequency = frequency,
    model = "bcms",
    density_sw = density_sw,
    sound_speed_sw = sound_speed_sw
  )
  gas_explicit <- target_strength(
    gas_object,
    frequency = frequency,
    model = "bcms",
    boundary = "gas_filled",
    density_sw = density_sw,
    sound_speed_sw = sound_speed_sw
  )

  expect_equal(liquid_default@model$BCMS$TS, liquid_explicit@model$BCMS$TS, tolerance = 1e-12)
  expect_equal(gas_default@model$BCMS$TS, gas_explicit@model$BCMS$TS, tolerance = 1e-12)

  expect_warning(
    target_strength(
      fls_generate(
        shape = cylinder(
          length_body = 0.03,
          radius_body = 0.003,
          radius_curvature_ratio = 1.8,
          n_segments = 101
        ),
        density_body = density_sw * 1.03,
        sound_speed_body = sound_speed_sw * 1.02,
        theta_body = 0.7 * pi
      ),
      frequency = 38e3,
      model = "bcms",
      boundary = "liquid_filled",
      density_sw = density_sw,
      sound_speed_sw = sound_speed_sw
    ),
    "broadside or near-broadside incidence"
  )

  expect_warning(
    target_strength(
      liquid_object,
      frequency = 38e3,
      model = "bcms",
      boundary = "liquid_filled",
      density_sw = density_sw,
      sound_speed_sw = sound_speed_sw,
      stationary_phase = TRUE
    ),
    "'stationary_phase' is ignored"
  )

  expect_error(
    target_strength(
      fls_generate(
        shape = sphere(radius_body = 0.01, n_segments = 80),
        density_body = density_sw * 1.03,
        sound_speed_body = sound_speed_sw * 1.02
      ),
      frequency = 38e3,
      model = "bcms",
      boundary = "liquid_filled",
      density_sw = density_sw,
      sound_speed_sw = sound_speed_sw
    ),
    "requires scatterer to be shape-type 'Cylinder'"
  )
})

test_that("brake() propagates cylinder curvature metadata into BCMS inputs", {
  density_sw <- 1026.8
  sound_speed_sw <- 1477.3

  straight <- fls_generate(
    shape = cylinder(
      length_body = 10.5e-3,
      radius_body = 1e-3,
      n_segments = 401
    ),
    density_body = density_sw * 1.0357,
    sound_speed_body = sound_speed_sw * 1.0279,
    theta_body = pi / 2
  )

  bent <- brake(straight, radius_curvature = 1.5)

  expect_equal(bent@shape_parameters$radius_curvature_ratio, 1.5)
  expect_equal(bent@body$radius_curvature_ratio, 1.5)

  bent_bcms <- target_strength(
    object = bent,
    frequency = 38e3,
    model = "bcms",
    density_sw = density_sw,
    sound_speed_sw = sound_speed_sw
  )

  expect_true(
    isTRUE(bent_bcms@model_parameters$BCMS$body$is_bent)
  )
  expect_equal(
    bent_bcms@model_parameters$BCMS$body$radius_curvature,
    1.5 * bent_bcms@shape_parameters$length
  )
})

test_that("ECMS reproduces straight and bent elastic-cylinder coherence bookkeeping", {
  length_body <- 0.04
  radius_body <- 0.005
  radius_curvature_ratio <- 1.8
  freqs <- c(38e3, 120e3, 200e3)

  straight <- fls_generate(
    shape = cylinder(
      length_body = length_body,
      radius_body = radius_body,
      n_segments = 201
    ),
    density_body = 2800,
    sound_speed_body = 1500,
    theta_body = pi / 2
  )

  bent <- fls_generate(
    shape = cylinder(
      length_body = length_body,
      radius_body = radius_body,
      radius_curvature_ratio = radius_curvature_ratio,
      n_segments = 201
    ),
    density_body = 2800,
    sound_speed_body = 1500,
    theta_body = pi / 2
  )

  straight_out <- target_strength(
    straight,
    model = "ECMS",
    frequency = freqs,
    density_sw = 1026.8,
    sound_speed_sw = 1477.3,
    sound_speed_longitudinal_body = 6398,
    sound_speed_transversal_body = 3122
  )

  bent_out <- target_strength(
    bent,
    model = "ECMS",
    frequency = freqs,
    density_sw = 1026.8,
    sound_speed_sw = 1477.3,
    sound_speed_longitudinal_body = 6398,
    sound_speed_transversal_body = 3122
  )

  lebc <- .ecms_equivalent_length_fresnel(
    k1 = 2 * pi * freqs / 1477.3,
    l = length_body,
    a = radius_body,
    rho_c = radius_curvature_ratio * length_body
  )

  bent_reference_ts <- 20 * log10(
    abs(lebc * straight_out@model$ECMS$f_bs / length_body)
  )

  expect_true("ECMS" %in% names(straight_out@model))
  expect_true(all(is.finite(straight_out@model$ECMS$TS)))
  expect_true(all(straight_out@model$ECMS$sigma_bs > 0))
  expect_equal(
    bent_out@model$ECMS$TS,
    bent_reference_ts,
    tolerance = 1e-10
  )
})

test_that("ECMS accepts ESS cylinders directly", {
  elastic_shell_cylinder <- ess_generate(
    shape = cylinder(
      length_body = 0.04,
      radius_body = 0.005,
      n_segments = 201
    ),
    density_shell = 2800,
    theta_shell = pi / 2
  )

  out <- target_strength(
    elastic_shell_cylinder,
    frequency = c(38e3, 120e3, 200e3),
    model = "ECMS",
    density_sw = 1026.8,
    sound_speed_sw = 1477.3,
    sound_speed_longitudinal_body = 6398,
    sound_speed_transversal_body = 3122
  )

  expect_true("ECMS" %in% names(out@model))
  expect_true(all(is.finite(out@model$ECMS$TS)))
  expect_true(all(out@model$ECMS$sigma_bs > 0))
})

test_that("ECMS enforces required elastic inputs and limited-angle scope", {
  density_sw <- 1026.8
  sound_speed_sw <- 1477.3
  frequency <- c(38e3, 120e3)

  density_missing <- fls_generate(
    shape = cylinder(
      length_body = 0.04,
      radius_body = 0.005,
      n_segments = 201
    ),
    g_body = 1.02,
    h_body = 1.01,
    theta_body = pi / 2
  )
  speed_missing <- fls_generate(
    shape = cylinder(
      length_body = 0.04,
      radius_body = 0.005,
      n_segments = 201
    ),
    density_body = 2800,
    sound_speed_body = 1500,
    theta_body = pi / 2
  )

  expect_error(
    target_strength(
      density_missing,
      frequency = frequency,
      model = "ecms",
      density_sw = density_sw,
      sound_speed_sw = sound_speed_sw,
      sound_speed_longitudinal_body = 6398,
      sound_speed_transversal_body = 3122
    ),
    "requires the elastic-cylinder density"
  )

  expect_error(
    target_strength(
      speed_missing,
      frequency = frequency,
      model = "ecms",
      density_sw = density_sw,
      sound_speed_sw = sound_speed_sw
    ),
    "requires both the longitudinal and transversal wave speeds"
  )

  expect_warning(
    target_strength(
      fls_generate(
        shape = cylinder(
          length_body = 0.04,
          radius_body = 0.005,
          n_segments = 201
        ),
        density_body = 2800,
        sound_speed_body = 1500,
        theta_body = 0.7 * pi
      ),
      frequency = frequency,
      model = "ecms",
      density_sw = density_sw,
      sound_speed_sw = sound_speed_sw,
      sound_speed_longitudinal_body = 6398,
      sound_speed_transversal_body = 3122
    ),
    "broadside or near-broadside incidence"
  )

  expect_error(
    suppressWarnings(target_strength(
      fls_generate(
        shape = cylinder(
          length_body = 0.04,
          radius_body = 0.005,
          n_segments = 201
        ),
        density_body = 2800,
        sound_speed_body = 1500,
        theta_body = 0
      ),
      frequency = 38e3,
      model = "ecms",
      density_sw = density_sw,
      sound_speed_sw = sound_speed_sw,
      sound_speed_longitudinal_body = 6398,
      sound_speed_transversal_body = 3122
    )),
    "not implemented for end-on incidence"
  )

  expect_error(
    target_strength(
      fls_generate(
        shape = sphere(radius_body = 0.01, n_segments = 80),
        density_body = 2800,
        sound_speed_body = 1500,
        theta_body = pi / 2
      ),
      frequency = 38e3,
      model = "ecms",
      density_sw = density_sw,
      sound_speed_sw = sound_speed_sw,
      density_body = 2800,
      sound_speed_longitudinal_body = 6398,
      sound_speed_transversal_body = 3122
    ),
    "requires scatterer to be shape-type 'Cylinder'"
  )
})
