library(acousticTS)

test_that("TMM matches SPHMS across multiple sphere sizes and boundary cases", {
  density_sw <- 1026.8
  sound_speed_sw <- 1477.3
  frequency <- c(12e3, 38e3, 120e3, 180e3)
  radii <- c(0.005, 0.01, 0.018)

  compare_sphere_case <- function(object, boundary, tolerance = 1e-9) {
    tmm_obj <- target_strength(
      object = object,
      frequency = frequency,
      model = "tmm",
      boundary = boundary,
      density_sw = density_sw,
      sound_speed_sw = sound_speed_sw
    )
    sphms_obj <- target_strength(
      object = object,
      frequency = frequency,
      model = "sphms",
      boundary = boundary,
      density_sw = density_sw,
      sound_speed_sw = sound_speed_sw
    )

    expect_equal(
      tmm_obj@model$TMM$f_bs,
      sphms_obj@model$SPHMS$f_bs,
      tolerance = tolerance
    )
    expect_equal(
      tmm_obj@model$TMM$TS,
      sphms_obj@model$SPHMS$TS,
      tolerance = tolerance
    )
  }

  for (radius in radii) {
    compare_sphere_case(
      fls_generate(
        shape = sphere(radius_body = radius, n_segments = 80),
        g_body = 1,
        h_body = 1
      ),
      boundary = "fixed_rigid"
    )

    compare_sphere_case(
      fls_generate(
        shape = sphere(radius_body = radius, n_segments = 80),
        g_body = 1,
        h_body = 1
      ),
      boundary = "pressure_release"
    )

    compare_sphere_case(
      fls_generate(
        shape = sphere(radius_body = radius, n_segments = 80),
        g_body = 1028.9 / density_sw,
        h_body = 1480.3 / sound_speed_sw
      ),
      boundary = "liquid_filled"
    )

    compare_sphere_case(
      gas_generate(
        shape = sphere(radius_body = radius, n_segments = 80),
        density_fluid = 1.24,
        sound_speed_fluid = 345
      ),
      boundary = "gas_filled"
    )
  }
})

test_that("TMM can retain the frequency-specific T-matrix blocks", {
  density_sw <- 1026.8
  sound_speed_sw <- 1477.3
  object <- fls_generate(
    shape = sphere(radius_body = 0.01, n_segments = 80),
    g_body = 1,
    h_body = 1
  )

  object <- target_strength(
    object = object,
    frequency = 38e3,
    model = "tmm",
    boundary = "pressure_release",
    density_sw = density_sw,
    sound_speed_sw = sound_speed_sw,
    store_t_matrix = TRUE
  )

  t_store <- object@model_parameters$TMM$parameters$t_matrix
  expect_length(t_store, 1)
  expect_true(is.list(t_store[[1]]))
  expect_length(t_store[[1]], object@model$TMM$n_max[1] + 1)
  expect_true(all(vapply(t_store[[1]], function(x) !is.null(x$T), logical(1))))
})

test_that("TMM supports oblate spheroids through the spherical-coordinate branch", {
  density_sw <- 1026.8
  sound_speed_sw <- 1477.3
  frequency <- c(12e3, 38e3, 70e3, 120e3)

  # Sphere-limit sanity check: when the polar and equatorial radii match, the
  # oblate shape should reduce to the exact sphere response.
  oblate_sphere <- fls_generate(
    shape = oblate_spheroid(length_body = 0.02, radius_body = 0.01, n_segments = 80),
    density_body = 1028.9,
    sound_speed_body = 1480.3,
    theta_body = pi / 2
  )
  sphere_object <- fls_generate(
    shape = sphere(radius_body = 0.01, n_segments = 80),
    density_body = 1028.9,
    sound_speed_body = 1480.3,
    theta_body = pi / 2
  )

  oblate_tmm <- target_strength(
    object = oblate_sphere,
    frequency = frequency,
    model = "tmm",
    boundary = "liquid_filled",
    density_sw = density_sw,
    sound_speed_sw = sound_speed_sw,
    store_t_matrix = TRUE
  )
  sphere_sphms <- target_strength(
    object = sphere_object,
    frequency = frequency,
    model = "sphms",
    boundary = "liquid_filled",
    density_sw = density_sw,
    sound_speed_sw = sound_speed_sw
  )

  expect_equal(oblate_tmm@model$TMM$TS, sphere_sphms@model$SPHMS$TS, tolerance = 1e-9)

  t_store <- oblate_tmm@model_parameters$TMM$parameters$t_matrix
  expect_length(t_store, length(frequency))
  expect_true(all(vapply(t_store, is.list, logical(1))))
  expect_true(all(vapply(t_store[[1]], function(x) !is.null(x$T), logical(1))))
})

test_that("Oblate TMM retained blocks support the same post-processing helpers", {
  density_sw <- 1026.8
  sound_speed_sw <- 1477.3
  object <- target_strength(
    object = fls_generate(
      shape = oblate_spheroid(length_body = 0.012, radius_body = 0.01, n_segments = 80),
      density_body = 1028.9,
      sound_speed_body = 1480.3,
      theta_body = pi / 2
    ),
    frequency = c(38e3, 70e3),
    model = "tmm",
    boundary = "liquid_filled",
    density_sw = density_sw,
    sound_speed_sw = sound_speed_sw,
    store_t_matrix = TRUE
  )

  mono <- tmm_scattering(object)
  avg <- tmm_average_orientation(
    object,
    distribution = tmm_orientation_distribution(
      distribution = "uniform",
      lower = 0.45 * pi,
      upper = pi,
      n_theta = 9
    )
  )
  grid <- tmm_scattering_grid(
    object,
    frequency = 70e3,
    theta_body = pi / 2,
    phi_body = pi,
    theta_scatter = pi / 2,
    phi_scatter = 2 * pi
  )
  summary <- tmm_bistatic_summary(
    object,
    frequency = 70e3,
    n_theta = 21,
    n_phi = 41,
    n_psi = 41
  )
  bundle <- tmm_products(
    object,
    frequency = 70e3,
    orientation = tmm_orientation_distribution(
      distribution = "uniform",
      lower = 0.45 * pi,
      upper = pi,
      n_theta = 9
    ),
    bistatic_summary = TRUE,
    n_theta = 21,
    n_phi = 41,
    n_psi = 41
  )

  expect_true(all(is.finite(mono$sigma_scat)))
  expect_true(all(is.finite(avg$sigma_bs)))
  expect_equal(grid$sigma_scat[1, 1], mono$sigma_scat[2], tolerance = 1e-8)
  expect_equal(summary$metrics$sigma_bs, mono$sigma_scat[2], tolerance = 1e-10)
  expect_equal(
    bundle$bistatic_summary$metrics$peak_sigma_scat,
    summary$metrics$peak_sigma_scat,
    tolerance = 1e-10
  )
})

test_that("TMM matches PSMS across multiple prolate spheroids and boundary cases", {
  density_sw <- 1026.8
  sound_speed_sw <- 1477.3
  shape_specs <- list(
    list(
      shape = prolate_spheroid(
        length_body = 0.06,
        radius_body = 0.008,
        n_segments = 80
      ),
      frequency = c(12e3, 18e3, 38e3, 70e3, 100e3, 150e3)
    ),
    list(
      shape = prolate_spheroid(
        length_body = 0.14,
        radius_body = 0.01,
        n_segments = 80
      ),
      frequency = c(12e3, 18e3, 38e3, 70e3, 100e3)
    )
  )

  compare_prolate_case <- function(object, boundary, frequency, max_delta) {
    tmm_obj <- target_strength(
      object = object,
      frequency = frequency,
      model = "tmm",
      boundary = boundary,
      density_sw = density_sw,
      sound_speed_sw = sound_speed_sw
    )
    psms_obj <- target_strength(
      object = object,
      frequency = frequency,
      model = "psms",
      boundary = boundary,
      density_sw = density_sw,
      sound_speed_sw = sound_speed_sw,
      precision = if (boundary %in% c("liquid_filled", "gas_filled")) "quad" else "double",
      simplify_Amn = FALSE
    )

    expect_true(all(is.finite(tmm_obj@model$TMM$TS)))
    expect_true(all(is.finite(tmm_obj@model$TMM$n_max)))
    expect_true(all(tmm_obj@model$TMM$n_max >= 1))
    expect_lt(max(abs(tmm_obj@model$TMM$TS - psms_obj@model$PSMS$TS)), max_delta)
  }

  for (spec in shape_specs) {
    compare_prolate_case(
      fls_generate(shape = spec$shape, g_body = 1, h_body = 1),
      boundary = "fixed_rigid",
      frequency = spec$frequency,
      max_delta = 1e-10
    )

    compare_prolate_case(
      fls_generate(shape = spec$shape, g_body = 1, h_body = 1),
      boundary = "pressure_release",
      frequency = spec$frequency,
      max_delta = 1e-10
    )

    compare_prolate_case(
      fls_generate(
        shape = spec$shape,
        g_body = 1028.9 / density_sw,
        h_body = 1480.3 / sound_speed_sw
      ),
      boundary = "liquid_filled",
      frequency = spec$frequency,
      max_delta = 1e-10
    )

    compare_prolate_case(
      gas_generate(
        shape = spec$shape,
        density_fluid = 1.24,
        sound_speed_fluid = 345
      ),
      boundary = "gas_filled",
      frequency = spec$frequency,
      max_delta = 1e-10
    )
  }
})

test_that("TMM warns when prolate spheroids request spherical-branch options", {
  density_sw <- 1026.8
  sound_speed_sw <- 1477.3
  object <- fls_generate(
    shape = prolate_spheroid(length_body = 0.14, radius_body = 0.01, n_segments = 80),
    g_body = 1,
    h_body = 1
  )

  expect_warning(
    target_strength(
      object,
      frequency = 38e3,
      model = "tmm",
      boundary = "fixed_rigid",
      density_sw = density_sw,
      sound_speed_sw = sound_speed_sw,
      n_max = 30,
      store_t_matrix = TRUE
    ),
    "n_max' is ignored"
  )
})

test_that("TMM can retain the frequency-specific T-matrix blocks for prolates", {
  density_sw <- 1026.8
  sound_speed_sw <- 1477.3
  object <- fls_generate(
    shape = prolate_spheroid(length_body = 0.14, radius_body = 0.01, n_segments = 80),
    density_body = 1028.9,
    sound_speed_body = 1480.3,
    theta_body = pi / 2
  )

  object <- target_strength(
    object = object,
    frequency = c(38e3, 100e3),
    model = "tmm",
    boundary = "liquid_filled",
    density_sw = density_sw,
    sound_speed_sw = sound_speed_sw,
    store_t_matrix = TRUE
  )

  t_store <- object@model_parameters$TMM$parameters$t_matrix
  expect_length(t_store, 2)
  expect_true(all(vapply(t_store, is.list, logical(1))))
  expect_true(all(vapply(t_store[[1]], function(x) !is.null(x$T), logical(1))))
})

test_that("Stored prolate TMM keeps the spheroidal modal cutoff", {
  density_sw <- 1026.8
  sound_speed_sw <- 1480
  object <- fls_generate(
    shape = prolate_spheroid(length_body = 0.07, radius_body = 0.01, n_segments = 80),
    density_body = 1024,
    sound_speed_body = 1480,
    theta_body = pi / 2
  )

  object <- target_strength(
    object = object,
    frequency = 38e3,
    model = "tmm",
    boundary = "pressure_release",
    density_sw = density_sw,
    sound_speed_sw = sound_speed_sw,
    store_t_matrix = TRUE
  )

  acoustics <- object@model_parameters$TMM$parameters$acoustics
  expected_n_max <- acoustics$m_max + ceiling(0.5 * acoustics$chi_sw)

  expect_equal(acoustics$n_max, expected_n_max)
})

test_that("tmm_scattering reproduces stored monostatic results", {
  density_sw <- 1026.8
  sound_speed_sw <- 1477.3

  sphere_object <- target_strength(
    object = fls_generate(
      shape = sphere(radius_body = 0.01, n_segments = 80),
      density_body = 1028.9,
      sound_speed_body = 1480.3,
      theta_body = pi / 2
    ),
    frequency = c(38e3, 70e3),
    model = "tmm",
    boundary = "liquid_filled",
    density_sw = density_sw,
    sound_speed_sw = sound_speed_sw,
    store_t_matrix = TRUE
  )
  sphere_scat <- tmm_scattering(sphere_object)
  expect_true(all(is.finite(sphere_scat$sigma_scat)))

  prolate_object <- target_strength(
    object = fls_generate(
      shape = prolate_spheroid(length_body = 0.14, radius_body = 0.01, n_segments = 80),
      density_body = 1028.9,
      sound_speed_body = 1480.3,
      theta_body = pi / 2
    ),
    frequency = c(38e3, 70e3),
    model = "tmm",
    boundary = "liquid_filled",
    density_sw = density_sw,
    sound_speed_sw = sound_speed_sw,
    store_t_matrix = TRUE
  )
  prolate_scat <- tmm_scattering(prolate_object)
  expect_equal(prolate_scat$f_scat, prolate_object@model$TMM$f_bs, tolerance = 1e-10)
})

test_that("Stored and unstored prolate TMM runs agree and stay tied to the PSMS benchmark", {
  density_sw <- 1026.8
  sound_speed_sw <- 1477.3
  frequency <- c(38e3, 70e3, 100e3)
  object <- fls_generate(
    shape = prolate_spheroid(length_body = 0.14, radius_body = 0.01, n_segments = 80),
    density_body = 1028.9,
    sound_speed_body = 1480.3,
    theta_body = pi / 2
  )

  tmm_direct <- target_strength(
    object = object,
    frequency = frequency,
    model = "tmm",
    boundary = "liquid_filled",
    density_sw = density_sw,
    sound_speed_sw = sound_speed_sw,
    store_t_matrix = FALSE
  )
  tmm_stored <- target_strength(
    object = object,
    frequency = frequency,
    model = "tmm",
    boundary = "liquid_filled",
    density_sw = density_sw,
    sound_speed_sw = sound_speed_sw,
    store_t_matrix = TRUE
  )
  psms_obj <- target_strength(
    object = object,
    frequency = frequency,
    model = "psms",
    boundary = "liquid_filled",
    density_sw = density_sw,
    sound_speed_sw = sound_speed_sw,
    precision = "quad",
    simplify_Amn = FALSE
  )
  expected_f_bs <- (-2i / wavenumber(frequency, sound_speed_sw)) * psms_obj@model$PSMS$f_bs

  expect_lt(max(Mod(tmm_direct@model$TMM$f_bs - expected_f_bs)), 1e-6)
  expect_lt(max(Mod(tmm_stored@model$TMM$f_bs - expected_f_bs)), 1e-6)
  expect_lt(max(Mod(tmm_stored@model$TMM$f_bs - tmm_direct@model$TMM$f_bs)), 1e-6)
  expect_equal(tmm_stored@model$TMM$sigma_bs, tmm_direct@model$TMM$sigma_bs, tolerance = 1e-6)
  expect_equal(tmm_stored@model$TMM$TS, tmm_direct@model$TMM$TS, tolerance = 5e-4)
})

test_that("Stored prolate TMM matches the exact general-angle pressure-release field", {
  density_sw <- 1026.8
  sound_speed_sw <- 1480
  theta_body <- pi / 2
  phi_body <- pi / 2
  scatter_angles <- data.frame(
    theta_scatter = c(pi / 2, pi / 2, pi / 2, pi / 3, 2 * pi / 3),
    phi_scatter = c(0, pi / 2, pi, pi / 2, 3 * pi / 2)
  )

  object <- target_strength(
    object = fls_generate(
      shape = prolate_spheroid(length_body = 0.07, radius_body = 0.01, n_segments = 80),
      g_body = 1,
      h_body = 1,
      theta_body = theta_body
    ),
    frequency = c(38e3, 70e3),
    model = "tmm",
    boundary = "pressure_release",
    density_sw = density_sw,
    sound_speed_sw = sound_speed_sw,
    store_t_matrix = TRUE
  )

  acoustics <- object@model_parameters$TMM$parameters$acoustics
  body <- object@model_parameters$TMM$body
  medium <- object@model_parameters$TMM$medium
  quad <- gauss_legendre(n = 96, a = -1, b = 1)

  for (i in seq_len(nrow(scatter_angles))) {
    body$theta_body <- theta_body
    body$phi_body <- phi_body
    body$theta_scatter <- scatter_angles$theta_scatter[i]
    body$phi_scatter <- scatter_angles$phi_scatter[i]

    exact_raw <- prolate_spheroid_fbs(
      acoustics = acoustics,
      body = body,
      medium = medium,
      integration_pts = quad,
      precision = "double",
      Amn_method = "Amn_pressure_release"
    )
    exact_f <- (-2i / acoustics$k_sw) * exact_raw
    stored_f <- tmm_scattering(
      object,
      theta_body = theta_body,
      phi_body = phi_body,
      theta_scatter = body$theta_scatter,
      phi_scatter = body$phi_scatter
    )$f_scat

    expect_equal(stored_f, exact_f, tolerance = 1e-10)
  }
})

test_that("Stored prolate TMM matches the exact general-angle rigid field", {
  density_sw <- 1026.8
  sound_speed_sw <- 1480
  theta_body <- pi / 2
  phi_body <- pi / 2
  scatter_angles <- data.frame(
    theta_scatter = c(pi / 2, pi / 2, pi / 2, pi / 3, 2 * pi / 3),
    phi_scatter = c(0, pi / 2, pi, pi / 2, 3 * pi / 2)
  )

  object <- target_strength(
    object = fls_generate(
      shape = prolate_spheroid(length_body = 0.07, radius_body = 0.01, n_segments = 80),
      g_body = 1,
      h_body = 1,
      theta_body = theta_body
    ),
    frequency = c(38e3, 70e3),
    model = "tmm",
    boundary = "fixed_rigid",
    density_sw = density_sw,
    sound_speed_sw = sound_speed_sw,
    store_t_matrix = TRUE
  )

  acoustics <- object@model_parameters$TMM$parameters$acoustics
  body <- object@model_parameters$TMM$body
  medium <- object@model_parameters$TMM$medium
  quad <- gauss_legendre(n = 96, a = -1, b = 1)

  for (i in seq_len(nrow(scatter_angles))) {
    body$theta_body <- theta_body
    body$phi_body <- phi_body
    body$theta_scatter <- scatter_angles$theta_scatter[i]
    body$phi_scatter <- scatter_angles$phi_scatter[i]

    exact_raw <- prolate_spheroid_fbs(
      acoustics = acoustics,
      body = body,
      medium = medium,
      integration_pts = quad,
      precision = "double",
      Amn_method = "Amn_fixed_rigid"
    )
    exact_f <- (-2i / acoustics$k_sw) * exact_raw
    stored_f <- tmm_scattering(
      object,
      theta_body = theta_body,
      phi_body = phi_body,
      theta_scatter = body$theta_scatter,
      phi_scatter = body$phi_scatter
    )$f_scat

    expect_equal(stored_f, exact_f, tolerance = 1e-10)
  }
})

test_that("tmm_average_orientation reduces to the supplied single angle", {
  density_sw <- 1026.8
  sound_speed_sw <- 1477.3
  object <- target_strength(
    object = fls_generate(
      shape = prolate_spheroid(length_body = 0.14, radius_body = 0.01, n_segments = 80),
      density_body = 1028.9,
      sound_speed_body = 1480.3,
      theta_body = pi / 2
    ),
    frequency = c(38e3, 70e3),
    model = "tmm",
    boundary = "liquid_filled",
    density_sw = density_sw,
    sound_speed_sw = sound_speed_sw,
    store_t_matrix = TRUE
  )

  one_angle <- tmm_scattering(object, theta_body = pi / 2)
  avg_angle <- tmm_average_orientation(object, theta_body = pi / 2)
  expect_equal(avg_angle$sigma_bs, one_angle$sigma_scat, tolerance = 1e-10)
  expect_equal(avg_angle$TS, db(avg_angle$sigma_bs), tolerance = 1e-12)
})

test_that("TMM orientation distributions normalize and agree with explicit averaging", {
  density_sw <- 1026.8
  sound_speed_sw <- 1477.3
  object <- target_strength(
    object = fls_generate(
      shape = prolate_spheroid(length_body = 0.14, radius_body = 0.01, n_segments = 80),
      density_body = 1028.9,
      sound_speed_body = 1480.3,
      theta_body = pi / 2
    ),
    frequency = c(38e3, 70e3),
    model = "tmm",
    boundary = "liquid_filled",
    density_sw = density_sw,
    sound_speed_sw = sound_speed_sw,
    store_t_matrix = TRUE
  )

  quadrature_dist <- tmm_orientation_distribution(
    distribution = "quadrature",
    theta_body = c(0.6 * pi, 0.75 * pi, pi),
    weights = c(1, 2, 1),
    phi_body = pi
  )
  expect_equal(sum(quadrature_dist$weights), 1, tolerance = 1e-12)

  avg_from_dist <- tmm_average_orientation(object, distribution = quadrature_dist)
  avg_explicit <- tmm_average_orientation(
    object,
    theta_body = quadrature_dist$theta_body,
    weights = quadrature_dist$weights,
    phi_body = quadrature_dist$phi_body
  )
  expect_equal(avg_from_dist$sigma_bs, avg_explicit$sigma_bs, tolerance = 1e-10)
  expect_equal(avg_from_dist$TS, avg_explicit$TS, tolerance = 1e-10)

  pdf_dist <- tmm_orientation_distribution(
    distribution = "pdf",
    theta_body = seq(0.5 * pi, pi, length.out = 9),
    pdf = function(theta) sin(theta)^2
  )
  expect_equal(sum(pdf_dist$weights), 1, tolerance = 1e-12)
  expect_true(all(pdf_dist$weights >= 0))

  trunc_dist <- tmm_orientation_distribution(
    distribution = "truncated_normal",
    mean_theta = 0.75 * pi,
    sd_theta = 0.08 * pi,
    lower = 0.5 * pi,
    upper = pi,
    n_theta = 11
  )
  expect_equal(sum(trunc_dist$weights), 1, tolerance = 1e-12)
  expect_true(all(trunc_dist$theta_body >= 0.5 * pi))
  expect_true(all(trunc_dist$theta_body <= pi))
})

test_that("TMM scattering plots work from stored blocks", {
  suppressPlot <- function(expr) {
    pdf(NULL)
    on.exit(dev.off())
    invisible(force(expr))
  }

  density_sw <- 1026.8
  sound_speed_sw <- 1477.3

  sphere_object <- target_strength(
    object = fls_generate(
      shape = sphere(radius_body = 0.01, n_segments = 80),
      density_body = 1028.9,
      sound_speed_body = 1480.3,
      theta_body = pi / 2
    ),
    frequency = c(38e3, 70e3),
    model = "tmm",
    boundary = "liquid_filled",
    density_sw = density_sw,
    sound_speed_sw = sound_speed_sw,
    store_t_matrix = TRUE
  )
  expect_error(
    suppressPlot(
      plot(
        sphere_object,
        type = "scattering",
        frequency = 38e3,
        vary = "theta_scatter"
      )
    ),
    NA
  )
  expect_error(
    suppressPlot(
      plot(
        sphere_object,
        type = "scattering",
        frequency = 38e3,
        heatmap = TRUE,
        n_theta = 31,
        n_phi = 61
      )
    ),
    NA
  )

  prolate_object <- target_strength(
    object = fls_generate(
      shape = prolate_spheroid(length_body = 0.14, radius_body = 0.01, n_segments = 80),
      density_body = 1028.9,
      sound_speed_body = 1480.3,
      theta_body = pi / 2
    ),
    frequency = c(38e3, 70e3),
    model = "tmm",
    boundary = "liquid_filled",
    density_sw = density_sw,
    sound_speed_sw = sound_speed_sw,
    store_t_matrix = TRUE
  )
  expect_error(
    suppressPlot(
      plot(
        prolate_object,
        type = "scattering",
        frequency = 70e3,
        vary = "phi_scatter",
        quantity = "sigma_scat"
      )
    ),
    NA
  )
  expect_error(
    suppressPlot(
      plot(
        prolate_object,
        type = "scattering",
        frequency = 70e3,
        polar = TRUE,
        n_theta = 31,
        n_phi = 61
      )
    ),
    NA
  )
})

test_that("TMM scattering plots require stored blocks", {
  density_sw <- 1026.8
  sound_speed_sw <- 1477.3

  object <- target_strength(
    object = fls_generate(
      shape = sphere(radius_body = 0.01, n_segments = 80),
      density_body = 1028.9,
      sound_speed_body = 1480.3,
      theta_body = pi / 2
    ),
    frequency = 38e3,
    model = "tmm",
    boundary = "liquid_filled",
    density_sw = density_sw,
    sound_speed_sw = sound_speed_sw
  )

  expect_error(
    plot(object, type = "scattering", frequency = 38e3),
    "Stored T-matrix blocks are required"
  )
})

test_that("tmm_scattering_grid reproduces stored monostatic scattering", {
  density_sw <- 1026.8
  sound_speed_sw <- 1477.3

  object <- target_strength(
    object = fls_generate(
      shape = prolate_spheroid(length_body = 0.14, radius_body = 0.01, n_segments = 80),
      density_body = 1028.9,
      sound_speed_body = 1480.3,
      theta_body = pi / 2
    ),
    frequency = c(38e3, 70e3),
    model = "tmm",
    boundary = "liquid_filled",
    density_sw = density_sw,
    sound_speed_sw = sound_speed_sw,
    store_t_matrix = TRUE
  )

  mono <- tmm_scattering(
    object,
    theta_body = pi / 2,
    phi_body = pi,
    theta_scatter = pi / 2,
    phi_scatter = 2 * pi
  )
  grid <- tmm_scattering_grid(
    object,
    frequency = 70e3,
    theta_body = pi / 2,
    phi_body = pi,
    theta_scatter = pi / 2,
    phi_scatter = 2 * pi
  )

  expect_equal(grid$frequency, 70e3)
  expect_equal(grid$sigma_scat[1, 1], mono$sigma_scat[2], tolerance = 1e-10)
  expect_equal(grid$sigma_scat_dB[1, 1], mono$sigma_scat_dB[2], tolerance = 1e-10)
})

test_that("stored spherical TMM point scattering stays on the monostatic amplitude scale", {
  density_sw <- 1026.8
  sound_speed_sw <- 1477.3

  object <- target_strength(
    object = fls_generate(
      shape = sphere(radius_body = 0.01, n_segments = 80),
      g_body = 1,
      h_body = 1,
      theta_body = pi / 2
    ),
    frequency = c(38e3, 70e3),
    model = "tmm",
    boundary = "pressure_release",
    density_sw = density_sw,
    sound_speed_sw = sound_speed_sw,
    store_t_matrix = TRUE
  )

  mono <- tmm_scattering(
    object,
    theta_body = pi / 2,
    phi_body = pi / 2,
    theta_scatter = pi / 2,
    phi_scatter = 3 * pi / 2
  )
  grid <- tmm_scattering_grid(
    object,
    frequency = 70e3,
    theta_body = pi / 2,
    phi_body = pi / 2,
    theta_scatter = pi / 2,
    phi_scatter = 3 * pi / 2
  )

  expect_equal(
    mono$f_scat[2],
    object@model$TMM$f_bs[2],
    tolerance = 1e-10
  )
  expect_equal(
    mono$sigma_scat[2],
    object@model$TMM$sigma_bs[2],
    tolerance = 1e-10
  )
  expect_equal(
    grid$f_scat[1, 1],
    mono$f_scat[2],
    tolerance = 1e-10
  )
})

test_that("TMM bistatic summaries and product bundles are internally consistent", {
  density_sw <- 1026.8
  sound_speed_sw <- 1477.3

  object <- target_strength(
    object = fls_generate(
      shape = prolate_spheroid(length_body = 0.14, radius_body = 0.01, n_segments = 80),
      density_body = 1028.9,
      sound_speed_body = 1480.3,
      theta_body = pi / 2
    ),
    frequency = c(38e3, 70e3, 100e3),
    model = "tmm",
    boundary = "liquid_filled",
    density_sw = density_sw,
    sound_speed_sw = sound_speed_sw,
    store_t_matrix = TRUE
  )

  summary_70 <- tmm_bistatic_summary(
    object,
    frequency = 70e3,
    n_theta = 21,
    n_phi = 41,
    n_psi = 41,
    include_grid = TRUE
  )
  monostatic_70 <- tmm_scattering(
    object,
    theta_body = pi / 2,
    phi_body = pi,
    theta_scatter = pi / 2,
    phi_scatter = 2 * pi
  )

  expect_equal(summary_70$metrics$sigma_bs, monostatic_70$sigma_scat[2], tolerance = 1e-10)
  expect_equal(summary_70$metrics$TS, db(summary_70$metrics$sigma_bs), tolerance = 1e-12)
  expect_equal(
    summary_70$metrics$peak_sigma_scat,
    max(summary_70$grid$sigma_scat),
    tolerance = 1e-10
  )
  expect_true(is.finite(summary_70$metrics$backscatter_lobe_width))
  expect_true(all(summary_70$sector_integrals$integrated_sigma_scat >= 0))

  orientation <- tmm_orientation_distribution(
    distribution = "uniform",
    lower = 0.5 * pi,
    upper = pi,
    n_theta = 9
  )
  products <- tmm_products(
    object,
    frequency = 70e3,
    orientation = orientation,
    bistatic_summary = TRUE,
    n_theta = 21,
    n_phi = 41,
    n_psi = 41
  )

  expect_equal(products$monostatic$f_scat, tmm_scattering(object)$f_scat, tolerance = 1e-10)
  expect_equal(
    products$orientation_average$sigma_bs,
    tmm_average_orientation(object, distribution = orientation)$sigma_bs,
    tolerance = 1e-10
  )
  expect_equal(
    products$bistatic_summary$metrics$peak_sigma_scat,
    summary_70$metrics$peak_sigma_scat,
    tolerance = 1e-10
  )
})

test_that("Default cylinder TMM matches FCMS through the cylindrical backend", {
  old_opt <- options(acousticTS.warn_tmm_cylinder = FALSE)
  on.exit(options(old_opt), add = TRUE)
  density_sw <- 1026.8
  sound_speed_sw <- 1477.3
  frequency <- c(12e3, 18e3, 38e3, 70e3, 100e3, 150e3, 200e3)

  compare_cylinder_case <- function(object, boundary, tolerance = 1e-12) {
    tmm_obj <- target_strength(
      object = object,
      frequency = frequency,
      model = "tmm",
      boundary = boundary,
      density_sw = density_sw,
      sound_speed_sw = sound_speed_sw
    )
    fcms_obj <- target_strength(
      object = object,
      frequency = frequency,
      model = "fcms",
      boundary = boundary,
      density_sw = density_sw,
      sound_speed_sw = sound_speed_sw
    )

    expect_true(all(is.finite(tmm_obj@model$TMM$TS)))
    expect_true(all(tmm_obj@model_parameters$TMM$parameters$coordinate_system == "cylindrical"))
    expect_equal(tmm_obj@model$TMM$f_bs, fcms_obj@model$FCMS$f_bs, tolerance = tolerance)
    expect_equal(tmm_obj@model$TMM$TS, fcms_obj@model$FCMS$TS, tolerance = tolerance)
  }

  cylinder_fls <- fls_generate(
    shape = cylinder(length_body = 0.07, radius_body = 0.01, n_segments = 80),
    density_body = 1028.9,
    sound_speed_body = 1480.3,
    theta_body = pi / 2
  )
  cylinder_gas <- gas_generate(
    shape = cylinder(length_body = 0.07, radius_body = 0.01, n_segments = 80),
    density_fluid = 1.24,
    sound_speed_fluid = 345,
    theta_body = pi / 2
  )

  compare_cylinder_case(cylinder_fls, "fixed_rigid")
  compare_cylinder_case(cylinder_fls, "pressure_release")
  compare_cylinder_case(cylinder_fls, "liquid_filled")
  compare_cylinder_case(cylinder_gas, "gas_filled")
})

test_that("Cylinder TMM warns that the branch is still experimental", {
  density_sw <- 1026.8
  sound_speed_sw <- 1477.3

  expect_warning(
    target_strength(
      object = fls_generate(
        shape = cylinder(length_body = 0.07, radius_body = 0.01, n_segments = 80),
        density_body = 1028.9,
        sound_speed_body = 1480.3,
        theta_body = pi / 2
      ),
      frequency = 38e3,
      model = "tmm",
      boundary = "pressure_release",
      density_sw = density_sw,
      sound_speed_sw = sound_speed_sw
    ),
    "Cylinder 'TMM' support remains experimental"
  )
})

test_that("Cylinder TMM stores cylindrical-family retained state when requested", {
  old_opt <- options(acousticTS.warn_tmm_cylinder = FALSE)
  on.exit(options(old_opt), add = TRUE)
  density_sw <- 1026.8
  sound_speed_sw <- 1477.3
  object <- target_strength(
    object = fls_generate(
      shape = cylinder(length_body = 0.07, radius_body = 0.01, n_segments = 80),
      density_body = 1028.9,
      sound_speed_body = 1480.3,
      theta_body = pi / 2
    ),
    frequency = c(38e3, 70e3),
    model = "tmm",
    boundary = "fixed_rigid",
    density_sw = density_sw,
      sound_speed_sw = sound_speed_sw,
      store_t_matrix = TRUE
  )

  t_store <- object@model_parameters$TMM$parameters$t_matrix

  expect_identical(object@model_parameters$TMM$parameters$coordinate_system, "cylindrical")
  expect_true(is.list(t_store))
  expect_length(t_store, 2)
  expect_true(all(vapply(t_store, is.list, logical(1))))
  expect_true(all(vapply(t_store, function(x) identical(x$family, "cylindrical_mono"), logical(1))))
})

test_that("Default cylinder TMM matches FCMS across incident-angle sweeps", {
  old_opt <- options(acousticTS.warn_tmm_cylinder = FALSE)
  on.exit(options(old_opt), add = TRUE)
  density_sw <- 1026.8
  sound_speed_sw <- 1477.3
  frequency <- c(38e3, 70e3, 120e3, 200e3)
  theta_vals <- c(0.55 * pi, 0.7 * pi, 0.85 * pi, pi - 0.05)

  compare_angle_case <- function(boundary, gas = FALSE) {
    for (theta_i in theta_vals) {
      object <- if (gas) {
        gas_generate(
          shape = cylinder(length_body = 0.07, radius_body = 0.01, n_segments = 80),
          density_fluid = 1.24,
          sound_speed_fluid = 345,
          theta_body = theta_i
        )
      } else {
        fls_generate(
          shape = cylinder(length_body = 0.07, radius_body = 0.01, n_segments = 80),
          density_body = 1028.9,
          sound_speed_body = 1480.3,
          theta_body = theta_i
        )
      }

      tmm_obj <- target_strength(
        object = object,
        frequency = frequency,
        model = "tmm",
        boundary = boundary,
        density_sw = density_sw,
        sound_speed_sw = sound_speed_sw
      )
      fcms_obj <- target_strength(
        object = object,
        frequency = frequency,
        model = "fcms",
        boundary = boundary,
        density_sw = density_sw,
        sound_speed_sw = sound_speed_sw
      )

      expect_equal(tmm_obj@model$TMM$f_bs, fcms_obj@model$FCMS$f_bs, tolerance = 1e-12)
      expect_equal(tmm_obj@model$TMM$TS, fcms_obj@model$FCMS$TS, tolerance = 1e-12)
    }
  }

  compare_angle_case("fixed_rigid")
  compare_angle_case("pressure_release")
  compare_angle_case("liquid_filled")
  compare_angle_case("gas_filled", gas = TRUE)
})

test_that("Cylinder TMM retained state supports exact monostatic reuse only", {
  old_opt <- options(acousticTS.warn_tmm_cylinder = FALSE)
  on.exit(options(old_opt), add = TRUE)
  density_sw <- 1026.8
  sound_speed_sw <- 1477.3
  object <- target_strength(
    object = fls_generate(
      shape = cylinder(length_body = 0.07, radius_body = 0.01, n_segments = 80),
      density_body = 1028.9,
      sound_speed_body = 1480.3,
      theta_body = pi / 2
    ),
    frequency = c(38e3, 70e3),
    model = "tmm",
    boundary = "fixed_rigid",
    density_sw = density_sw,
      sound_speed_sw = sound_speed_sw,
      store_t_matrix = TRUE
  )

  mono <- tmm_scattering(object)
  avg <- tmm_average_orientation(
    object,
    distribution = tmm_orientation_distribution(
      distribution = "uniform",
      lower = 0.45 * pi,
      upper = pi,
      n_theta = 9
    )
  )
  bundle <- tmm_products(
    object,
    frequency = 70e3,
    orientation = tmm_orientation_distribution(
      distribution = "uniform",
      lower = 0.45 * pi,
      upper = pi,
      n_theta = 9
    ),
    bistatic_summary = FALSE
  )

  expect_true(all(is.finite(mono$sigma_scat)))
  expect_true(all(is.finite(avg$sigma_bs)))
  expect_equal(mono$f_scat, object@model$TMM$f_bs, tolerance = 1e-12)
  expect_equal(bundle$monostatic$f_scat, mono$f_scat, tolerance = 1e-12)
  expect_equal(
    bundle$orientation_average$sigma_bs,
    avg$sigma_bs,
    tolerance = 1e-12
  )
  expect_error(
    tmm_scattering_grid(
      object,
      frequency = 70e3,
      theta_body = pi / 2,
      phi_body = pi,
      theta_scatter = pi / 2,
      phi_scatter = 2 * pi
    ),
    "Stored cylindrical TMM grid evaluations are not available yet"
  )
  expect_error(
    tmm_bistatic_summary(
      object,
      frequency = 70e3,
      n_theta = 21,
      n_phi = 41,
      n_psi = 41
    ),
    "Stored cylindrical TMM bistatic summaries are not available yet"
  )
  expect_error(
    tmm_products(
      object,
      frequency = 70e3,
      orientation = tmm_orientation_distribution(
        distribution = "uniform",
        lower = 0.45 * pi,
        upper = pi,
        n_theta = 9
      ),
      bistatic_summary = TRUE,
      n_theta = 21,
      n_phi = 41,
      n_psi = 41
    ),
    "Stored cylindrical TMM bistatic summaries are not available yet"
  )
})

test_that("Cylinder TMM rejects unsupported general-angle stored scattering requests", {
  old_opt <- options(acousticTS.warn_tmm_cylinder = FALSE)
  on.exit(options(old_opt), add = TRUE)
  density_sw <- 1026.8
  sound_speed_sw <- 1477.3

  object <- target_strength(
    object = fls_generate(
      shape = cylinder(length_body = 0.07, radius_body = 0.01, n_segments = 80),
      density_body = 1028.9,
      sound_speed_body = 1480.3,
      theta_body = pi / 2
    ),
    frequency = c(38e3, 70e3),
    model = "tmm",
    boundary = "fixed_rigid",
    density_sw = density_sw,
    sound_speed_sw = sound_speed_sw,
    store_t_matrix = TRUE
  )

  expect_error(
    tmm_scattering(
      object,
      theta_body = pi / 2,
      phi_body = pi,
      theta_scatter = pi / 2,
      phi_scatter = pi / 2
    ),
    "currently available only for the exact monostatic direction"
  )

  expect_error(
    suppressWarnings(plot(
      object,
      type = "scattering",
      frequency = 70e3,
      polar = TRUE,
      n_theta = 21,
      n_phi = 41
    )),
    "Stored cylindrical TMM grid evaluations are not available yet"
  )
})

test_that("tmm_diagnostics provides physics-based checks for stored TMM objects", {
  old_opt <- options(acousticTS.warn_tmm_cylinder = FALSE)
  on.exit(options(old_opt), add = TRUE)
  density_sw <- 1026.8
  sound_speed_sw <- 1477.3

  sphere_obj <- target_strength(
    object = fls_generate(
      shape = sphere(radius_body = 0.01, n_segments = 80),
      g_body = 1,
      h_body = 1,
      theta_body = pi / 2
    ),
    frequency = c(38e3, 70e3),
    model = "tmm",
    boundary = "fixed_rigid",
    density_sw = density_sw,
    sound_speed_sw = sound_speed_sw,
    store_t_matrix = TRUE
  )

  prolate_obj <- target_strength(
    object = fls_generate(
      shape = prolate_spheroid(
        length_body = 0.14,
        radius_body = 0.01,
        n_segments = 80
      ),
      density_body = 1028.9,
      sound_speed_body = 1480.3,
      theta_body = pi / 2
    ),
    frequency = c(38e3, 70e3),
    model = "tmm",
    boundary = "liquid_filled",
    density_sw = density_sw,
    sound_speed_sw = sound_speed_sw,
    store_t_matrix = TRUE
  )

  cylinder_obj <- target_strength(
    object = fls_generate(
      shape = cylinder(length_body = 0.07, radius_body = 0.01, n_segments = 80),
      density_body = 1028.9,
      sound_speed_body = 1480.3,
      theta_body = pi / 2
    ),
    frequency = c(38e3, 70e3),
    model = "tmm",
    boundary = "fixed_rigid",
    density_sw = density_sw,
      sound_speed_sw = sound_speed_sw,
      store_t_matrix = TRUE
  )

  sphere_diag <- tmm_diagnostics(sphere_obj, n_theta = 31, n_phi = 61)
  prolate_diag <- tmm_diagnostics(prolate_obj, n_theta = 31, n_phi = 61)
  cylinder_diag <- tmm_diagnostics(cylinder_obj, n_theta = 31, n_phi = 61)

  expect_true(all(c("shape", "coordinate_system", "boundary", "frequency") %in% names(sphere_diag$summary)))
  expect_equal(length(sphere_diag$block_metrics), 2L)
  expect_true(all(is.finite(sphere_diag$summary$monostatic_rel_residual)))
  expect_true(all(sphere_diag$summary$reciprocity_rel_residual < 1e-10))
  expect_true(all(is.finite(sphere_diag$summary$optical_theorem_rel_residual)))
  expect_true(all(is.na(sphere_diag$summary$continuation_max_abs_second_diff_TS)))

  expect_true(all(prolate_diag$summary$monostatic_rel_residual < 1e-10))
  expect_true(all(prolate_diag$summary$reciprocity_rel_residual < 1e-5))
  expect_true(all(is.na(prolate_diag$summary$optical_theorem_rel_residual)))
  expect_true(all(is.finite(prolate_diag$summary$continuation_max_abs_step_TS)))
  expect_true(all(is.finite(prolate_diag$summary$continuation_max_abs_second_diff_TS)))
  expect_false(any(prolate_diag$summary$continuation_any_nonfinite))
  expect_true(is.data.frame(prolate_diag$continuation))
  expect_equal(min(prolate_diag$continuation$aspect_ratio), 1, tolerance = 1e-12)
  expect_equal(
    max(prolate_diag$continuation$aspect_ratio),
    prolate_diag$summary$continuation_target_aspect_ratio[1],
    tolerance = 1e-12
  )

  expect_true(all(is.finite(cylinder_diag$summary$monostatic_rel_residual)))
  expect_true(all(is.na(cylinder_diag$summary$reciprocity_rel_residual)))
  expect_true(all(is.na(cylinder_diag$summary$optical_theorem_rel_residual)))
  expect_true(all(is.na(cylinder_diag$summary$min_block_rcond)))
  expect_true(all(is.na(cylinder_diag$summary$max_block_transpose_residual)))
  expect_true(all(vapply(cylinder_diag$block_metrics, nrow, integer(1)) == 0L))
  expect_true(all(is.na(cylinder_diag$summary$continuation_max_abs_second_diff_TS)))
})

test_that("TMM sphere-to-spheroid continuation starts from the exact equal-volume sphere", {
  density_sw <- 1026.8
  sound_speed_sw <- 1480

  prolate_obj <- target_strength(
    object = fls_generate(
      shape = prolate_spheroid(length_body = 0.07, radius_body = 0.01, n_segments = 80),
      g_body = 1,
      h_body = 1,
      theta_body = pi / 2
    ),
    frequency = c(38e3, 70e3),
    model = "tmm",
    boundary = "pressure_release",
    density_sw = density_sw,
    sound_speed_sw = sound_speed_sw,
    store_t_matrix = TRUE
  )
  diag <- tmm_diagnostics(prolate_obj, continuation_steps = 5L, n_theta = 21, n_phi = 41)
  first_step <- diag$continuation[diag$continuation$step == 1, ]

  r_eq <- ((0.07 / 2) * 0.01^2)^(1 / 3)
  sphere_obj <- target_strength(
    object = fls_generate(
      shape = sphere(radius_body = r_eq, n_segments = 80),
      g_body = 1,
      h_body = 1,
      theta_body = pi / 2
    ),
    frequency = c(38e3, 70e3),
    model = "tmm",
    boundary = "pressure_release",
    density_sw = density_sw,
    sound_speed_sw = sound_speed_sw
  )

  expect_true(all(first_step$shape == "Sphere"))
  expect_equal(first_step$TS, sphere_obj@model$TMM$TS, tolerance = 1e-10)
})

test_that("TMM rejects unsupported shapes", {
  object <- fls_generate(
    shape = arbitrary(
      x_body = seq(0, 0.04, length.out = 50),
      zU_body = 0.005 * sin(seq(0, pi, length.out = 50)),
      zL_body = -0.005 * sin(seq(0, pi, length.out = 50)),
      n_segments = 50
    ),
    density_body = 1028.9,
    sound_speed_body = 1480.3
  )

  expect_error(
    target_strength(object, frequency = 38e3, model = "tmm", boundary = "liquid_filled"),
    "supports the following shape-types"
  )
})
