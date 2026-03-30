library(acousticTS)

test_that("TMM scattering helper vectors, weights, and orientation utilities cover edge branches", {
  expect_error(
    acousticTS:::.tmm_angle_vector(
      value = NULL,
      default = NULL,
      name = "theta_body"
    ),
    "'theta_body' could not be resolved."
  )
  expect_equal(
    acousticTS:::.tmm_angle_vector(
      value = NULL,
      default = c(0, pi),
      n_default = 3,
      lower = 0,
      upper = pi,
      name = "theta_body"
    ),
    c(0, pi / 2, pi),
    tolerance = 1e-12
  )
  expect_equal(
    acousticTS:::.tmm_angle_vector(
      value = NULL,
      default = c(0, pi),
      name = "theta_body"
    ),
    c(0, pi)
  )
  expect_error(
    acousticTS:::.tmm_angle_vector(c(-0.1, 0.2), lower = 0, name = "theta_body"),
    "'theta_body' must be >= 0 radians."
  )
  expect_error(
    acousticTS:::.tmm_angle_vector(c(0.1, 3.5), upper = pi, name = "theta_body"),
    "'theta_body' must be <= "
  )

  expect_equal(acousticTS:::.tmm_interval_weights(0.5), 1)
  expect_error(
    acousticTS:::.tmm_interval_weights("bad"),
    "'theta_body' must be a non-empty numeric vector of finite values."
  )
  expect_equal(
    acousticTS:::.tmm_interval_weights(c(0, 1, 2), lower = 0, upper = 2),
    c(0.5, 1, 0.5),
    tolerance = 1e-12
  )
  expect_error(
    acousticTS:::.tmm_interval_weights(c(1, 1), name = "theta_body"),
    "'theta_body' must be strictly increasing."
  )
  expect_equal(acousticTS:::.tmm_grid_edges(0.5, 0, 1), c(0, 1), tolerance = 1e-12)

  expect_error(
    acousticTS:::.tmm_validate_orientation_n_theta(0),
    "'n_theta' must be a single positive integer."
  )
  expect_equal(acousticTS:::.tmm_validate_orientation_n_theta(3), 3L)

  expect_error(
    acousticTS:::.tmm_validate_orientation_interval(-0.1, pi / 2, "uniform"),
    "'uniform' requires a finite interval"
  )
  expect_equal(
    acousticTS:::.tmm_validate_orientation_interval(0, pi / 2, "uniform")$upper,
    pi / 2,
    tolerance = 1e-12
  )

  expect_error(
    acousticTS:::.tmm_validate_orientation_normal(pi / 2, 0),
    "'mean_theta' and 'sd_theta' must be finite numeric scalars"
  )
  expect_equal(
    acousticTS:::.tmm_validate_orientation_normal(pi / 2, 0.2)$sd_theta,
    0.2,
    tolerance = 1e-12
  )

  quadrature_default <- acousticTS:::.tmm_orientation_quadrature(
    theta_body = c(0, pi / 2, pi),
    weights = NULL
  )
  expect_equal(quadrature_default$weights, rep(1 / 3, 3), tolerance = 1e-12)
  expect_error(
    acousticTS:::.tmm_orientation_quadrature(
      theta_body = c(0, pi / 2),
      weights = c(1, -1)
    ),
    "'weights' must be a non-negative numeric vector"
  )

  pdf_numeric <- acousticTS:::.tmm_orientation_pdf(
    theta_body = c(0, pi / 2, pi),
    pdf = c(1, 2, 1)
  )
  pdf_function <- acousticTS:::.tmm_orientation_pdf(
    theta_body = c(0, pi / 2, pi),
    pdf = function(theta) rep(1, length(theta))
  )
  expect_equal(sum(pdf_numeric$weights), 1, tolerance = 1e-12)
  expect_equal(sum(pdf_function$weights), 1, tolerance = 1e-12)

  uniform <- acousticTS:::.tmm_orientation_uniform(0, pi / 2, 5)
  normal_density <- acousticTS:::.tmm_orientation_normal_density(pi / 2, 0.2, 5)
  truncated <- acousticTS:::.tmm_orientation_truncated_normal(
    mean_theta = pi / 2,
    sd_theta = 0.2,
    lower = pi / 4,
    upper = 3 * pi / 4,
    n_theta = 5
  )
  expect_equal(sum(uniform$weights), 1, tolerance = 1e-12)
  expect_equal(sum(normal_density$weights), 1, tolerance = 1e-12)
  expect_equal(sum(truncated$weights), 1, tolerance = 1e-12)

  expect_error(
    acousticTS:::.tmm_resolve_orientation_phi("bad", 2),
    "'phi_body' must be a finite numeric scalar or vector."
  )
  expect_equal(
    acousticTS:::.tmm_resolve_orientation_phi(pi, 3),
    rep(pi, 3),
    tolerance = 1e-12
  )
})

test_that("TMM scattering helper frequency and cylindrical branches cover fallbacks", {
  expect_equal(acousticTS:::.tmm_plot_frequency_index(NULL, 38e3), 1L)
  expect_error(
    acousticTS:::.tmm_plot_frequency_index(NULL, c(38e3, 70e3)),
    "requires a scalar 'frequency' input"
  )
  expect_error(
    acousticTS:::.tmm_plot_frequency_index("bad", c(38e3, 70e3)),
    "'frequency' must be a single finite value in Hz."
  )
  expect_warning(
    idx <- acousticTS:::.tmm_plot_frequency_index(39e3, c(38e3, 70e3)),
    "Using the nearest stored frequency"
  )
  expect_equal(idx, 1L)

  local({
    testthat::local_mocked_bindings(
      .fcms_bm_fluid = function(...) matrix(c(1 + 0i, 2 + 0i), ncol = 1),
      .fcms_bm_fixed_rigid = function(...) c(1 + 0i, 2 + 0i),
      .fcms_bm_pressure_release = function(...) c(1 + 0i, 2 + 0i),
      .package = "acousticTS"
    )

    fluid_f <- acousticTS:::.tmm_cylindrical_monostatic_f_bs(
      acoustics_row = list(k_sw = 2, n_max = 1),
      body_defaults = list(g_body = 1.1, h_body = 1.2),
      shape_parameters = list(radius = c(0.01, 0.01), length = 0.04),
      boundary = "liquid_filled",
      theta_body = pi / 2
    )
    rigid_f <- acousticTS:::.tmm_cylindrical_monostatic_f_bs(
      acoustics_row = list(k_sw = 2, n_max = 1),
      body_defaults = list(g_body = 1.1, h_body = 1.2),
      shape_parameters = list(radius = c(0.01, 0.01), length = 0.04),
      boundary = "fixed_rigid",
      theta_body = pi / 2
    )
    prelease_f <- acousticTS:::.tmm_cylindrical_monostatic_f_bs(
      acoustics_row = list(k_sw = 2, n_max = 1),
      body_defaults = list(g_body = 1.1, h_body = 1.2),
      shape_parameters = list(radius = c(0.01, 0.01), length = 0.04),
      boundary = "pressure_release",
      theta_body = pi / 2
    )
    rigid_zero_length <- acousticTS:::.tmm_cylindrical_monostatic_f_bs(
      acoustics_row = list(k_sw = 2, n_max = 1),
      body_defaults = list(g_body = 1.1, h_body = 1.2),
      shape_parameters = list(radius = c(0.01, 0.01), length = 0),
      boundary = "fixed_rigid",
      theta_body = pi / 2
    )

    expect_true(is.complex(fluid_f))
    expect_true(is.complex(rigid_f))
    expect_true(is.complex(prelease_f))
    expect_true(is.complex(rigid_zero_length))
    expect_error(
      acousticTS:::.tmm_cylindrical_monostatic_f_bs(
        acoustics_row = list(k_sw = 2, n_max = 1),
        body_defaults = list(g_body = 1.1, h_body = 1.2),
        shape_parameters = list(radius = c(0.01, 0.01), length = 0.04),
        boundary = "unsupported",
        theta_body = pi / 2
      ),
      "Unsupported boundary for cylindrical TMM branch."
    )
  })

  expect_error(
    acousticTS:::.tmm_scattering_points(
      model_params = list(parameters = list(
        coordinate_system = "unsupported",
        acoustics = data.frame(frequency = 38e3)
      )),
      frequency_idx = 1L,
      shape_parameters = list(),
      theta_body = 0,
      phi_body = 0,
      theta_scatter = pi,
      phi_scatter = pi
    ),
    "Unsupported stored TMM coordinate system."
  )
})
