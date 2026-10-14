library(acousticTS)

test_that("full gas spheroids agree with the independent native reference", {
  object <- gas_generate(
    shape = prolate_spheroid(length_body = 0.2, radius_body = 0.01),
    density_fluid = 1.24, sound_speed_fluid = 345
  )
  # Prol_Spheroid vectorized 995929a with precision_fbs = 1e-8 =================
  reference <- c(-26.549891089759697, -24.53103140836401)
  precisions <- c("double", if (acousticTS:::.quad_precision_available()) "quad")
  for (precision in precisions) {
    # The higher-contrast retained system requires quad at the second case =====
    selected <- if (precision == "double") 1L else 1:2
    result <- target_strength(
      object, c(18e3, 44e3)[selected], "psms", boundary = "gas_filled",
      precision = precision, n_integration = 128,
      density_sw = 1026.8, sound_speed_sw = 1477.3
    )
    # The reference terminates its azimuthal sum before the retained limit =====
    expect_lt(max(abs(result@model$PSMS$TS - reference[selected])), 1e-4)
    expect_lt(abs(result@model$PSMS$TS[1] - reference[1]), 2e-8)
  }
})

test_that("full gas coefficient overlaps agree with converged quad quadrature", {
  skip_if_not(acousticTS:::.quad_precision_available())
  object <- gas_generate(
    shape = prolate_spheroid(length_body = 0.14, radius_body = 0.01),
    density_fluid = 1.24, sound_speed_fluid = 345
  )
  # Native quad Gauss-Legendre rules at 384 and 512 points agree here ==========
  result <- target_strength(
    object, 200e3, "psms", boundary = "gas_filled", precision = "quad",
    density_sw = 1026.8, sound_speed_sw = 1477.3
  )
  expect_equal(result@model$PSMS$TS, -29.0909564966414, tolerance = 1e-10)

  # Full gas overlaps do not depend on the liquid quadrature setting ===========
  coarse <- target_strength(
    object, 38e3, "psms", boundary = "gas_filled", precision = "quad",
    n_integration = 8, density_sw = 1026.8, sound_speed_sw = 1477.3
  )
  fine <- target_strength(
    object, 38e3, "psms", boundary = "gas_filled", precision = "quad",
    n_integration = 192, density_sw = 1026.8, sound_speed_sw = 1477.3
  )
  expect_identical(coarse@model$PSMS$TS, fine@model$PSMS$TS)
})

test_that("full fluid solves retain the final single-degree azimuthal block", {
  skip_if_not(acousticTS:::.quad_precision_available())
  expected <- c(liquid_filled = -77.20002046221784,
    gas_filled = -28.634043251185428)
  for (boundary in names(expected)) {
    shape <- prolate_spheroid(length_body = 0.14, radius_body = 0.01)
    object <- if (boundary == "gas_filled") {
      gas_generate(shape, density_fluid = 1.24, sound_speed_fluid = 345)
    } else {
      fls_generate(shape, density_body = 1028.9, sound_speed_body = 1480.3)
    }
    initialized <- acousticTS:::psms_initialize(
      object, 38000, boundary = boundary, precision = "quad",
      adaptive = FALSE, simplify_Amn = FALSE, n_integration = 96,
      density_sw = 1026.8, sound_speed_sw = 1477.3
    )
    model <- initialized@model_parameters$PSMS
    acoustics <- model$parameters$acoustics
    acoustics$n_max <- acoustics$m_max <- 48L
    field <- acousticTS:::prolate_spheroidal_kernels(
      acoustics, model$body, model$medium, model$parameters$Amn_method,
      n_integration = 96, precision = "quad", adaptive = FALSE
    )
    actual <- 20 * log10(Mod(-2i / acoustics$k_sw * field))
    expect_equal(as.numeric(actual), unname(expected[boundary]), tolerance = 2e-9)
  }
})

test_that("unavailable quad builds reject quad requests", {
  local_mocked_bindings(
    .quad_precision_available = function() FALSE,
    .package = "acousticTS"
  )
  expect_error(
    Smn(0, 0, 1, 0, precision = "quad"),
    "Native quad precision is unavailable"
  )
  expect_error(
    Rmn(0, 0, 1, 1.5, precision = "quad"),
    "Native quad precision is unavailable"
  )
  expect_equal(
    acousticTS:::.validate_quad_precision_available("double"),
    "double"
  )
})

test_that("spheroidal fluid strategies agree at modal convergence", {
  acoustics <- data.frame(chi_sw = 0.8, chi_body = 0.7, m_max = 8L, n_max = 12L)
  body <- data.frame(
    xi = 1.5,
    theta_body = pi / 3,
    theta_scatter = 2 * pi / 3,
    phi_body = 0,
    phi_scatter = pi,
    density = 1100
  )
  medium <- data.frame(density = 1000)
  quadrature <- gauss_legendre(48, -1, 1)
  precisions <- c(
    "double",
    if (acousticTS:::.quad_precision_available()) "quad"
  )
  for (precision in precisions) {
    for (method in c(
      "Amn_fluid",
      "Amn_fluid_gas",
      "Amn_fluid_simplify",
      "Amn_fixed_rigid",
      "Amn_pressure_release"
    )) {
      reference <- acousticTS:::prolate_spheroid_fbs(
        acoustics,
        body,
        medium,
        quadrature,
        precision,
        method,
        adaptive = FALSE,
        vectorized = FALSE
      )
      expect_true(all(is.finite(reference)))
      for (vectorized in c(FALSE, TRUE)) {
        actual <- acousticTS:::prolate_spheroid_fbs(
          acoustics,
          body,
          medium,
          quadrature,
          precision,
          method,
          adaptive = TRUE,
          vectorized = vectorized
        )
        expect_equal(actual, reference, tolerance = 1e-07)
      }
      vectorized <- acousticTS:::prolate_spheroid_fbs(
        acoustics,
        body,
        medium,
        quadrature,
        precision,
        method,
        adaptive = FALSE,
        vectorized = TRUE
      )
      expect_equal(vectorized, reference, tolerance = 1e-07)
      # A general receive angle bypasses the monostatic shortcuts. =============
      general_body <- body
      general_body$theta_scatter <- 0.7
      general_body$phi_scatter <- 0.4
      direct <- acousticTS:::prolate_spheroid_fbs(
        acoustics,
        general_body,
        medium,
        quadrature,
        precision,
        method,
        adaptive = FALSE,
        vectorized = FALSE
      )
      retained <- acousticTS:::prolate_spheroid_tmatrix_cpp(
        acoustics,
        general_body,
        medium,
        quadrature,
        precision,
        method
      )
      expect_equal(retained$f_scat, direct, tolerance = 1e-07)
    }
  }
})

test_that("quad monostatic kernels agree with retained bistatic amplitudes", {
  skip_if_not(acousticTS:::.quad_precision_available())
  medium <- data.frame(density = 1000)
  quadrature <- gauss_legendre(48, -1, 1)
  for (angle in c(0, pi / 6, pi / 2)) {
    acoustics <- data.frame(
      chi_sw = 5, chi_body = 4.375,
      m_max = if (angle == 0) 0L else 8L, n_max = 16L
    )
    body <- data.frame(
      xi = 1.5, theta_body = angle, theta_scatter = pi - angle,
      phi_body = pi, phi_scatter = 2 * pi, density = 1100
    )
    for (method in c(
      "Amn_fixed_rigid", "Amn_pressure_release", "Amn_fluid",
      "Amn_fluid_gas", "Amn_fluid_simplify"
    )) {
      retained <- acousticTS:::prolate_spheroid_tmatrix_cpp(
        acoustics, body, medium, quadrature, "quad", method
      )$f_scat
      for (adaptive in c(FALSE, TRUE)) {
        actual <- acousticTS:::prolate_spheroid_fbs(
          acoustics, body, medium, quadrature, "quad", method,
          adaptive = adaptive
        )
        expect_equal(actual, retained, tolerance = 1e-7,
          info = paste(method, angle, adaptive))
      }
    }
  }
})

test_that("near-backscatter azimuths retain their bistatic dependence", {
  # cos(pi + delta) can round within the former backscatter tolerance even =====
  # though delta is well resolved in the double-precision input angles. ========
  acoustics <- data.frame(chi_sw = 20, chi_body = 17.5, m_max = 20L, n_max = 35L)
  body <- data.frame(
    xi = 1.5, theta_body = pi / 2, theta_scatter = pi / 2,
    phi_body = 0, phi_scatter = pi + 1.5e-7, density = 1100
  )
  medium <- data.frame(density = 1000)
  quadrature <- gauss_legendre(48, -1, 1)
  actual <- acousticTS:::prolate_spheroid_fbs(
    acoustics, body, medium, quadrature, "double", "Amn_pressure_release"
  )
  retained <- acousticTS:::prolate_spheroid_tmatrix_cpp(
    acoustics, body, medium, quadrature, "double", "Amn_pressure_release"
  )$f_scat
  expect_equal(actual, retained, tolerance = 1e-13)
})

test_that("axial spheroidal scattering agrees under reversal", {
  # Axial incidence excites only m = 0 =========================================
  acoustics <- data.frame(chi_sw = 0.8, chi_body = 0.7, m_max = 0L, n_max = 12L)
  body <- data.frame(
    xi = 1.5,
    theta_body = 0,
    theta_scatter = pi,
    phi_body = 0,
    phi_scatter = pi,
    density = 1100
  )
  reversed <- body
  reversed$theta_body <- pi
  # Low component of pi represented as the sum of two doubles ==================
  reversed$theta_scatter <- 1.2246467991473532e-16
  medium <- data.frame(density = 1000)
  quadrature <- gauss_legendre(48, -1, 1)
  precisions <- c(
    "double",
    if (acousticTS:::.quad_precision_available()) "quad"
  )
  for (precision in precisions) {
    for (method in c(
      "Amn_fluid",
      "Amn_fluid_gas",
      "Amn_fluid_simplify",
      "Amn_fixed_rigid",
      "Amn_pressure_release"
    )) {
      reference <- acousticTS:::prolate_spheroid_fbs(
        acoustics,
        body,
        medium,
        quadrature,
        precision,
        method,
        adaptive = FALSE
      )
      expect_true(all(is.finite(reference)))
      for (adaptive in c(FALSE, TRUE)) {
        actual <- acousticTS:::prolate_spheroid_fbs(
          acoustics,
          reversed,
          medium,
          quadrature,
          precision,
          method,
          adaptive = adaptive
        )
        expect_equal(actual, reference, tolerance = 1e-07)
      }
      # Rigid, pressure-release and diagonal fluid formulations can retain =====
      # the unexcited azimuthal orders at exact axial incidence. ===============
      if (method != "Amn_fluid") {
        extra_modes <- acoustics
        extra_modes$m_max <- 4L
        actual <- acousticTS:::prolate_spheroid_fbs(
          extra_modes, reversed, medium, quadrature, precision, method,
          adaptive = FALSE
        )
        expect_equal(actual, reference, tolerance = 1e-7)
      }
    }
  }
})

test_that(
  paste0(
    "psms helper validators and quadrature selectors enforce the documented ",
    "contracts"
  ),
  {
    expect_error(
      acousticTS:::.psms_validate_shape(list(shape = "Sphere")),
      "requires scatterer to be shape-type 'ProlateSpheroid'"
    )
    expect_identical(
      acousticTS:::.psms_validate_boundary("gas_filled"),
      "gas_filled"
    )
    expect_error(
      acousticTS:::.psms_validate_boundary("invalid"),
      "Only the following values"
    )
    expect_identical(acousticTS:::.psms_validate_precision("double"), "double")
    expect_error(
      acousticTS:::.psms_validate_precision("single"),
      "must be either 'double' or 'quad'"
    )
    if (acousticTS:::.quad_precision_available()) {
      expect_identical(acousticTS:::.psms_validate_precision("quad"), "quad")
    } else {
      expect_error(
        acousticTS:::.psms_validate_precision("quad"),
        "Native quad precision is unavailable"
      )
    }
    expect_identical(acousticTS:::.psms_validate_adaptive(TRUE), TRUE)
    expect_error(
      acousticTS:::.psms_validate_adaptive(NA),
      "'adaptive' must be either TRUE or FALSE"
    )

    expect_identical(
      acousticTS:::.psms_Amn_method("liquid_filled", TRUE),
      "Amn_fluid_simplify"
    )
    expect_identical(
      acousticTS:::.psms_Amn_method("gas_filled", FALSE),
      "Amn_fluid_gas"
    )
    expect_identical(
      acousticTS:::.psms_Amn_method("fixed_rigid", FALSE),
      "Amn_fixed_rigid"
    )
    expect_true(acousticTS:::.psms_is_full_fluid_method("Amn_fluid"))
    expect_true(acousticTS:::.psms_is_full_fluid_method("Amn_fluid_gas"))
    expect_false(acousticTS:::.psms_is_full_fluid_method("Amn_fluid_simplify"))
    expect_identical(
      acousticTS:::.psms_n_integration_label("Amn_fluid"),
      "full liquid-filled PSMS solve"
    )
    expect_identical(
      acousticTS:::.psms_n_integration_label("Amn_fluid_gas"),
      "full gas-filled PSMS solve"
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
    state <- acousticTS:::.psms_body_state(obj,
      sound_speed_sw = 1500,
      density_sw = 1026
    )
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
      acousticTS:::.psms_n_integration(NULL,
        adaptive = FALSE, Amn_method =
          "Amn_fixed_rigid"
      ),
      96L
    )
    expect_warning(
      expect_true(is.na(acousticTS:::.psms_n_integration(
        48L,
        adaptive = TRUE,
        Amn_method = "Amn_fluid"
      ))),
      "full liquid-filled PSMS solve"
    )
    expect_identical(
      acousticTS:::.psms_n_integration(
        48L,
        adaptive = TRUE,
        Amn_method = "Amn_fluid_gas"
      ),
      48L
    )
    expect_error(
      acousticTS:::.psms_n_integration(1.5,
        adaptive = FALSE, Amn_method =
          "Amn_fluid"
      ),
      "'n_integration' must be a single positive integer"
    )

    gas_model_params <- within(model_params, {
      precision <- "quad"
      adaptive <- TRUE
      Amn_method <- "Amn_fluid_gas"
    })
    gas_promoted <- acousticTS:::.psms_promote_gas_modal_limits(
      model_params = acousticTS:::.psms_complete_acoustics(
        model_params = gas_model_params,
        body_params = body_params,
        scatterer_shape = shape
      ),
      body_params = body_params,
      boundary = "gas_filled"
    )
    liquid_unchanged <- acousticTS:::.psms_promote_gas_modal_limits(
      model_params = acousticTS:::.psms_complete_acoustics(
        model_params = gas_model_params,
        body_params = body_params,
        scatterer_shape = shape
      ),
      body_params = body_params,
      boundary = "liquid_filled"
    )
    expect_true(all(
      gas_promoted$acoustics$m_max >= liquid_unchanged$acoustics$m_max
    ))
    expect_true(all(
      gas_promoted$acoustics$n_max >= liquid_unchanged$acoustics$n_max
    ))

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

    skip_if_not(
      acousticTS:::.quad_precision_available(),
      "native quad precision is unavailable"
    )

    gas_obj <- gas_generate(
      shape = prolate_spheroid(
        length_body = 0.14,
        radius_body = 0.01,
        n_segments = 40
      ),
      density_fluid = 1.24,
      sound_speed_fluid = 345,
      theta_body = pi / 2
    )
    expect_no_warning(
      gas_init <- acousticTS:::psms_initialize(
        gas_obj,
        frequency = c(12000, 38000),
        boundary = "gas_filled",
        precision = "quad",
        simplify_Amn = FALSE,
        adaptive = TRUE,
        density_sw = 1026.8,
        sound_speed_sw = 1477.3
      )
    )
    expect_identical(
      gas_init@model_parameters$PSMS$parameters$Amn_method,
      "Amn_fluid_gas"
    )

    expect_no_warning(
      gas_fit_full <- target_strength(
        gas_obj,
        frequency = c(12000, 38000),
        model = "psms",
        boundary = "gas_filled",
        precision = "quad",
        simplify_Amn = FALSE,
        adaptive = TRUE,
        density_sw = 1026.8,
        sound_speed_sw = 1477.3,
        phi_body = pi
      )
    )
    gas_fit_simplified <- target_strength(
      gas_obj,
      frequency = c(12000, 38000),
      model = "psms",
      boundary = "gas_filled",
      precision = "quad",
      simplify_Amn = TRUE,
      adaptive = TRUE,
      density_sw = 1026.8,
      sound_speed_sw = 1477.3,
      phi_body = pi
    )
    expect_true(all(is.finite(gas_fit_full@model$PSMS$TS)))
    expect_gt(max(abs(gas_fit_full@model$PSMS$TS -
      gas_fit_simplified@model$PSMS$TS)), 1e-4)
    expect_equal(
      gas_fit_simplified@model$PSMS$TS,
      c(-30.162331, -28.635226),
      tolerance = 1e-5
    )
  }
)

test_that(
  paste0(
    "psms adaptive kernel helpers group repeated quadrature orders and ",
    "dispatch cleanly"
  ),
  {
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
      expect_equal(
        sort(vapply(calls, `[[`, integer(1), "n_integration")),
        c(24L, 32L)
      )
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
      adaptive_out_gas <- acousticTS:::prolate_spheroidal_kernels(
        acoustics = acoustics[1:2, , drop = FALSE],
        body = list(),
        medium = data.frame(),
        boundary_method = "Amn_fluid_gas",
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
      expect_equal(as.vector(adaptive_out_gas), rep(2 + 0i, 2))
      expect_identical(attr(adaptive_out_gas, "n_integration"), 96)
      expect_equal(as.vector(fixed_out), rep(2 + 0i, 2))
    })
  }
)
