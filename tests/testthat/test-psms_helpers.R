library(acousticTS)

test_that("psms helper validators and quadrature selectors enforce the documented contracts", {
  expect_error(
    acousticTS:::.psms_validate_shape(list(shape = "Sphere")),
    "requires scatterer to be shape-type 'ProlateSpheroid'"
  )
  expect_identical(acousticTS:::.psms_validate_boundary("gas_filled"), "gas_filled")
  expect_error(acousticTS:::.psms_validate_boundary("invalid"), "Only the following values")
  expect_identical(acousticTS:::.psms_validate_precision("double"), "double")
  expect_error(acousticTS:::.psms_validate_precision("single"), "must be either 'double' or 'quad'")
  expect_identical(acousticTS:::.psms_validate_adaptive(TRUE), TRUE)
  expect_error(acousticTS:::.psms_validate_adaptive(NA), "'adaptive' must be either TRUE or FALSE")

  expect_identical(
    acousticTS:::.psms_Amn_method("liquid_filled", TRUE),
    "Amn_fluid_simplify"
  )
  expect_identical(
    acousticTS:::.psms_Amn_method("fixed_rigid", FALSE),
    "Amn_fixed_rigid"
  )

  obj <- fls_generate(
    shape = prolate_spheroid(
      length_body = 0.08,
      radius_body = 0.01,
      n_segments = 40
    ),
    density_body = 1070,
    sound_speed_body = 1570,
    theta_body = 1.2
  )
  shape <- extract(obj, "shape_parameters")
  state <- acousticTS:::.psms_body_state(obj, sound_speed_sw = 1500, density_sw = 1026)
  model_params <- acousticTS:::.psms_model_parameters(
    frequency = c(38000, 70000),
    sound_speed_sw = 1500,
    body_h = state$body$h,
    precision = "double",
    n_integration = 96L,
    adaptive = FALSE
  )
  body_params <- acousticTS:::.psms_body_parameters(
    scatterer_shape = shape,
    body = state$body,
    phi_body = pi,
    sound_speed_sw = 1500,
    density_sw = 1026
  )

  expect_equal(model_params$acoustics$frequency, c(38000, 70000))
  expect_identical(model_params$precision, "double")
  expect_false(model_params$adaptive)
  expect_true(body_params$xi > 1)
  expect_equal(body_params$phi_scatter, 2 * pi)
  expect_equal(body_params$density, 1070, tolerance = 1e-10)
  expect_equal(body_params$sound_speed, 1570, tolerance = 1e-10)
  expect_true(body_params$q > 0)

  expect_identical(
    acousticTS:::.psms_n_integration(NULL, adaptive = FALSE, Amn_method = "Amn_fixed_rigid"),
    96L
  )
  expect_warning(
    expect_true(is.na(acousticTS:::.psms_n_integration(
      48L,
      adaptive = TRUE,
      Amn_method = "Amn_fluid"
    ))),
    "is ignored when 'adaptive = TRUE'"
  )
  expect_error(
    acousticTS:::.psms_n_integration(1.5, adaptive = FALSE, Amn_method = "Amn_fluid"),
    "'n_integration' must be a single positive integer"
  )

  adaptive_double <- acousticTS:::.psms_adaptive_n_integration(
    chi_sw = c(1, 400),
    chi_body = c(4, 900),
    m_max = c(3, 18),
    n_max = c(8, 40),
    precision = "double"
  )
  adaptive_quad <- acousticTS:::.psms_adaptive_n_integration(
    chi_sw = c(1, 400),
    chi_body = c(4, 900),
    m_max = c(3, 18),
    n_max = c(8, 40),
    precision = "quad"
  )

  expect_type(adaptive_double, "integer")
  expect_true(all(adaptive_double %% 8L == 0L))
  expect_true(all(adaptive_double >= 24L & adaptive_double <= 96L))
  expect_true(all(adaptive_quad %% 8L == 0L))
  expect_true(all(adaptive_quad >= 32L & adaptive_quad <= 96L))
  expect_true(all(adaptive_quad >= adaptive_double))
})

test_that("psms adaptive kernel helpers group repeated quadrature orders and dispatch cleanly", {
  calls <- list()
  acoustics <- data.frame(
    frequency = c(38000, 70000, 120000),
    chi_sw = c(1, 2, 3),
    chi_body = c(1, 2, 3),
    m_max = c(2, 2, 2),
    n_max = c(4, 4, 4)
  )

  local({
    testthat::local_mocked_bindings(
      .psms_adaptive_n_integration = function(...) c(24L, 32L, 24L),
      .prolate_spheroidal_kernels_fixed = function(acoustics,
                                                   body,
                                                   medium,
                                                   boundary_method,
                                                   n_integration,
                                                   precision,
                                                   adaptive) {
        calls[[length(calls) + 1L]] <<- list(
          frequencies = acoustics$frequency,
          n_integration = n_integration,
          precision = precision,
          adaptive = adaptive,
          boundary_method = boundary_method
        )
        rep(as.complex(n_integration), nrow(acoustics))
      },
      .package = "acousticTS"
    )

    grouped <- acousticTS:::.prolate_spheroidal_kernels_adaptive(
      acoustics = acoustics,
      body = list(),
      medium = data.frame(),
      boundary_method = "Amn_fluid",
      precision = "double",
      adaptive = TRUE
    )

    expect_equal(as.vector(grouped), c(24 + 0i, 32 + 0i, 24 + 0i))
    expect_equal(attr(grouped, "n_integration"), c(24L, 32L, 24L))
    expect_length(calls, 2)
    expect_equal(sort(vapply(calls, `[[`, integer(1), "n_integration")), c(24L, 32L))
  })

  local({
    testthat::local_mocked_bindings(
      .prolate_spheroidal_kernels_adaptive = function(...) rep(1 + 2i, 2),
      .prolate_spheroidal_kernels_fixed = function(acoustics,
                                                   body,
                                                   medium,
                                                   boundary_method,
                                                   n_integration,
                                                   precision,
                                                   adaptive) {
        structure(
          rep(2 + 0i, nrow(acoustics)),
          n_integration = n_integration,
          adaptive = adaptive,
          precision = precision,
          boundary_method = boundary_method
        )
      },
      .package = "acousticTS"
    )

    adaptive_out <- acousticTS:::prolate_spheroidal_kernels(
      acoustics = acoustics[1:2, , drop = FALSE],
      body = list(),
      medium = data.frame(),
      boundary_method = "Amn_fluid",
      precision = "double",
      adaptive = TRUE
    )
    fixed_out <- acousticTS:::prolate_spheroidal_kernels(
      acoustics = acoustics[1:2, , drop = FALSE],
      body = list(),
      medium = data.frame(),
      boundary_method = "Amn_fixed_rigid",
      n_integration = NA_integer_,
      precision = "quad",
      adaptive = FALSE
    )

    expect_equal(adaptive_out, rep(1 + 2i, 2))
    expect_equal(as.vector(fixed_out), rep(2 + 0i, 2))
  })
})
