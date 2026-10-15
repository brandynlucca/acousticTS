library(acousticTS)

test_that("prolate coordinates preserve exact semiaxes and legacy metadata", {
  geometry <- acousticTS:::.prolate_spheroidal_geometry
  legacy <- list(length = 0.14, radius = 0.01)
  expected_xi <- 1 / sqrt(1 - (legacy$radius / (legacy$length / 2))^2)
  expected_q <- legacy$length / 2 / expected_xi
  fallback <- geometry(legacy)
  expect_identical(fallback$xi, expected_xi)
  expect_identical(fallback$q, expected_q)
  exact <- c(legacy, list(semimajor_length = 0.07, semiminor_length = 0.01))
  expect_identical(geometry(exact), fallback)
  for (axis in c("semimajor_length", "semiminor_length")) {
    partial <- exact
    partial[[axis]] <- NULL
    expect_identical(geometry(partial), fallback)
  }

  # Invalid supplied values must not trigger the missing-metadata fallback =====
  invalid_values <- list(NULL, numeric(), NA_real_, NaN, Inf, -Inf, 0,
    -0.01, c(0.01, 0.02), "0.01", 0.01 + 0i)
  for (axis in c("semimajor_length", "semiminor_length")) {
    for (value in invalid_values) {
      invalid <- exact
      invalid[axis] <- list(value)
      expect_error(geometry(invalid), "finite scalar semiaxes")
    }
  }
  expect_error(geometry(list(length = 0.14, radius = 0.07)), "0 < b < a")
  expect_error(geometry(list(length = 0.14, radius = 0.08)), "0 < b < a")
  expect_error(geometry(list(length = 0.14, radius = 0)), "0 < b < a")

  # Coordinate representation does not certify subsequent radial accuracy ======
  expect_error(geometry(list(semimajor_length = 1, semiminor_length = 2^-27)),
    "R input arithmetic")
  near_sphere <- geometry(list(semimajor_length = 1,
    semiminor_length = 1 - .Machine$double.eps / 2))
  expect_identical(near_sphere$xi, 2^26)
  expect_identical(near_sphere$q, 2^-26)
})

test_that("PSMS geometry and native fields do not depend on profile sampling", {
  segments <- c(99L, 100L, 101L)
  objects <- lapply(segments, function(count) {
    fls_generate(prolate_spheroid(semimajor_length = 0.07,
      semiminor_length = 0.01, n_segments = count),
      g_body = 1.03, h_body = 1.02, theta_body = 37 * pi / 180)
  })
  sampled_radii <- vapply(objects, function(object) {
    extract(object, "shape_parameters")$radius
  }, numeric(1))
  expect_lt(sampled_radii[1], sampled_radii[2])
  expect_lt(sampled_radii[3], sampled_radii[2])
  precisions <- c("double", if (acousticTS:::.quad_precision_available()) "quad")
  for (precision in precisions) {
    initialized <- lapply(objects, function(object) {
      acousticTS:::psms_initialize(object, 12000, boundary = "pressure_release",
        precision = precision, density_sw = 1026.8, sound_speed_sw = 1477.3)
    })
    reference <- initialized[[2]]@model_parameters$PSMS
    expected_xi <- 1 / sqrt(1 - (0.01 / 0.07)^2)
    expect_identical(reference$body$xi, expected_xi)
    expect_identical(reference$body$q, 0.07 / expected_xi)
    expect_identical(reference$parameters$acoustics$m_max,
      ceiling(2 * reference$parameters$acoustics$k_sw * 0.01))
    for (index in c(1L, 3L)) {
      actual <- initialized[[index]]@model_parameters$PSMS
      expect_identical(actual$body, reference$body)
      expect_identical(actual$parameters$acoustics,
        reference$parameters$acoustics)
    }
    fields <- lapply(initialized, function(object) acousticTS:::PSMS(object))
    expect_true(all(is.finite(fields[[2]]@model$PSMS$f_bs)))
    expect_gt(Mod(fields[[2]]@model$PSMS$f_bs), 0)
    for (index in c(1L, 3L)) {
      expect_identical(fields[[index]]@model$PSMS$f_bs,
        fields[[2]]@model$PSMS$f_bs)
      expect_identical(extract(fields[[index]], "shape_parameters"),
        extract(objects[[index]], "shape_parameters"))
    }
  }
})

test_that("FLS reforge preserves prolate semiaxes through scaling and sampling", {
  object <- fls_generate(prolate_spheroid(semimajor_length = 0.07,
    semiminor_length = 0.01, n_segments = 99L),
    g_body = 1.03, h_body = 1.02, theta_body = 37 * pi / 180)
  resampled <- reforge(object, n_segments_body = 101L)
  resampled_shape <- extract(resampled, "shape_parameters")
  expect_identical(resampled_shape$semimajor_length, 0.07)
  expect_identical(resampled_shape$semiminor_length, 0.01)
  scales <- c(length = 1.5, radius = 0.75)
  resized <- reforge(resampled, body_scale = scales, isometric_body = FALSE)
  resized_shape <- extract(resized, "shape_parameters")
  expect_identical(resized_shape$semimajor_length, 0.07 * scales[["length"]])
  expect_identical(resized_shape$semiminor_length, 0.01 * scales[["radius"]])

  # The metadata correction leaves profile scaling unchanged for other models ==
  control <- resampled
  control@shape_parameters$shape <- "Arbitrary"
  control <- reforge(control, body_scale = scales, isometric_body = FALSE)
  expect_identical(resized@body$rpos, control@body$rpos)
  expect_identical(resized@body$radius, control@body$radius)
  expect_identical(resized_shape$length, control@shape_parameters$length)
  expect_identical(resized_shape$radius, control@shape_parameters$radius)

  # Compare the rescaled native field with a newly constructed exact spheroid ==
  rebuilt <- fls_generate(prolate_spheroid(
    semimajor_length = 0.07 * scales[["length"]],
    semiminor_length = 0.01 * scales[["radius"]], n_segments = 100L),
    g_body = 1.03, h_body = 1.02, theta_body = 37 * pi / 180)
  initialize <- function(value) {
    acousticTS:::psms_initialize(value, 12000, boundary = "pressure_release",
      precision = "double", density_sw = 1026.8, sound_speed_sw = 1477.3)
  }
  resized <- initialize(resized)
  rebuilt <- initialize(rebuilt)
  expect_identical(resized@model_parameters$PSMS$body,
    rebuilt@model_parameters$PSMS$body)
  expect_identical(resized@model_parameters$PSMS$parameters$acoustics,
    rebuilt@model_parameters$PSMS$parameters$acoustics)
  actual <- acousticTS:::PSMS(resized)
  reference <- acousticTS:::PSMS(rebuilt)
  expect_true(all(is.finite(reference@model$PSMS$f_bs)))
  expect_gt(Mod(reference@model$PSMS$f_bs), 0)
  expect_identical(actual@model$PSMS$f_bs, reference@model$PSMS$f_bs)

  # Older objects keep the explicit fallback when semiaxis metadata is absent ==
  legacy <- object
  legacy@shape_parameters$semimajor_length <- NULL
  legacy@shape_parameters$semiminor_length <- NULL
  legacy <- reforge(legacy, body_scale = scales, isometric_body = FALSE)
  expect_false("semimajor_length" %in% names(legacy@shape_parameters))
  expect_false("semiminor_length" %in% names(legacy@shape_parameters))
  geometry <- acousticTS:::.prolate_spheroidal_geometry(legacy@shape_parameters)
  expect_identical(geometry$semimajor_length, legacy@shape_parameters$length / 2)
  expect_identical(geometry$semiminor_length, legacy@shape_parameters$radius)
})

test_that("PSMS rejects degenerate prolate semiaxes before native evaluation", {
  object <- fls_generate(prolate_spheroid(semimajor_length = 0.07,
    semiminor_length = 0.01), g_body = 1.03, h_body = 1.02)
  for (axes in list(c(0.07, 0.07), c(0.07, 0.08), c(0, 0.01),
    c(0.07, 0), c(-0.07, 0.01), c(0.07, Inf), c(NA_real_, 0.01))) {
    invalid <- object
    invalid@shape_parameters$semimajor_length <- axes[1]
    invalid@shape_parameters$semiminor_length <- axes[2]
    expect_error(acousticTS:::psms_initialize(invalid, 12000),
      "finite scalar semiaxes")
  }
})

# Check degree and azimuthal enlargement separately in each fixed reference ====
psms_enlarged_reference <- function(acoustics, evaluate, tolerance = 1e-12,
                                   axial = FALSE) {
  enlarged <- acoustics
  enlarged$m_max <- enlarged$m_max + if (axial) 0L else 16L
  enlarged$n_max <- enlarged$n_max + 16L
  first <- evaluate(enlarged)
  enlarged$n_max <- enlarged$n_max + 16L
  more_degrees <- evaluate(enlarged)
  expect_lt(max(Mod(more_degrees - first) / Mod(more_degrees)), tolerance)
  if (!axial) {
    enlarged$m_max <- enlarged$m_max + 16L
    more_orders <- evaluate(enlarged)
    expect_lt(max(Mod(more_orders - more_degrees) / Mod(more_orders)), tolerance)
  } else {
    more_orders <- more_degrees
  }
  list(field = more_orders, acoustics = enlarged)
}

test_that("quad radial overflow at the R boundary is not misclassified as missing", {
  skip_if_not(acousticTS:::.quad_precision_available())
  xi <- 0.07 / sqrt(0.07^2 - 0.01^2)
  # Degree 450 has finite native quad components beyond the R double range =====
  second <- Rmn(0, c(448, 450), 120, xi, kind = 2, precision = "quad")
  scalar <- Rmn(0, 450, 120, xi, kind = 2, precision = "quad")
  expect_false(anyNA(unlist(second)))
  expect_false(anyNA(unlist(scalar)))
  expect_true(any(is.infinite(c(second$value[2], second$derivative[2]))))
  expect_equal(scalar$value, second$value[2])
  expect_equal(scalar$derivative, second$derivative[2])
  # Incoming and outgoing waves preserve the second-kind signs on conversion ===
  for (kind in c(3, 4)) {
    radial <- Rmn(0, c(448, 450), 120, xi, kind = kind, precision = "quad")
    sign <- if (kind == 3) 1 else -1
    expect_false(anyNA(unlist(radial)))
    expect_equal(Im(radial$value), sign * second$value)
    expect_equal(Im(radial$derivative), sign * second$derivative)
  }
})

test_that("PSMS selects coefficient overlaps for full gas and liquid", {
  object <- fls_generate(prolate_spheroid(length_body = 0.14, radius_body = 0.01),
    density_body = 1028.9, sound_speed_body = 1480.3)
  initialize <- function(...) acousticTS:::psms_initialize(object, 70000, ...)
  liquid <- initialize()
  expect_identical(liquid@model_parameters$PSMS$parameters$Amn_method,
    "Amn_fluid_coefficients")
  gas <- initialize(boundary = "gas_filled")
  expect_identical(gas@model_parameters$PSMS$parameters$Amn_method, "Amn_fluid_gas")
  simplified <- initialize(simplify_Amn = TRUE)
  expect_identical(simplified@model_parameters$PSMS$parameters$Amn_method,
    "Amn_fluid_simplify")
  expect_true(acousticTS:::.psms_is_full_fluid_method("Amn_fluid_coefficients"))
  expect_identical(acousticTS:::.psms_Amn_method("liquid_filled", 1, "quad"),
    "Amn_fluid_simplify")
  expect_identical(acousticTS:::.psms_Amn_method("liquid_filled", 0, "quad"),
    "Amn_fluid_coefficients")
  expect_identical(acousticTS:::.psms_Amn_method("liquid_filled", NA, "quad"),
    NA)
  if (acousticTS:::.quad_precision_available()) {
    coefficients <- initialize(precision = "quad")
    expect_identical(coefficients@model_parameters$PSMS$parameters$Amn_method,
      "Amn_fluid_coefficients")
    simplified_quad <- initialize(precision = "quad", simplify_Amn = TRUE)
    expect_identical(simplified_quad@model_parameters$PSMS$parameters$Amn_method,
      "Amn_fluid_simplify")
  }
})

test_that("full liquid PSMS retains the converged axial field", {
  object <- fls_generate(
    prolate_spheroid(length_body = 0.14, radius_body = 0.01),
    density_body = 1028.9, sound_speed_body = 1480.3, theta_body = 0)
  # AcousticScattering 941d077 with degrees 33, 49 and 65 agrees at this axis ==
  reference <- complex(real = -1.8480331274694501e-6,
    imaginary = 4.445317502670616e-8)
  precisions <- c("double", if (acousticTS:::.quad_precision_available()) "quad")
  for (precision in precisions) {
    result <- target_strength(object, 70000, "psms", boundary = "liquid_filled",
      precision = precision, density_sw = 1026.8, sound_speed_sw = 1477.3)
    inputs <- result@model_parameters$PSMS
    wave_number <- inputs$parameters$acoustics$k_sw
    amplitude <- -2i / wave_number * result@model$PSMS$f_bs
    tolerance <- if (precision == "quad") 1e-10 else 1e-8
    expect_lt(Mod(amplitude - reference) / Mod(reference), tolerance)
    if (precision == "quad") {
      # Check modal convergence independently of agreement at matching limits ==
      enlarged <- inputs$parameters$acoustics
      enlarged$m_max <- enlarged$m_max + 16L
      enlarged$n_max <- enlarged$n_max + 16L
      field <- acousticTS:::.prolate_spheroidal_kernels_fixed(
        acoustics = enlarged, body = inputs$body, medium = inputs$medium,
        boundary_method = inputs$parameters$Amn_method,
        n_integration = 96L, precision = precision, adaptive = FALSE)
      expect_lt(Mod(result@model$PSMS$f_bs - field) / Mod(field), 1e-10)
    }
  }
})

test_that("quad liquid coefficient overlaps preserve the full quadrature field", {
  skip_if_not(acousticTS:::.quad_precision_available())
  for (angle in c(90, 37)) {
    object <- fls_generate(prolate_spheroid(length_body = 0.14, radius_body = 0.01),
      density_body = 1028.9, sound_speed_body = 1480.3, theta_body = angle * pi / 180)
    evaluate <- function(n_integration = 96L, adaptive = FALSE) {
      target_strength(object, c(70000, 200000), "psms", boundary = "liquid_filled",
        precision = "quad",
        n_integration = n_integration, adaptive = adaptive,
        density_sw = 1026.8, sound_speed_sw = 1477.3)
    }
    coefficients <- evaluate()
    inputs <- coefficients@model_parameters$PSMS
    evaluate_reference <- function(acoustics, order = 192L) {
      acousticTS:::.prolate_spheroidal_kernels_fixed(
        acoustics = acoustics, body = inputs$body, medium = inputs$medium,
        boundary_method = "Amn_fluid", n_integration = order,
        precision = "quad", adaptive = FALSE)
    }
    reference <- psms_enlarged_reference(inputs$parameters$acoustics,
      evaluate_reference)
    more_nodes <- evaluate_reference(reference$acoustics, 256L)
    expect_lt(max(Mod(more_nodes - reference$field) / Mod(more_nodes)), 1e-12)
    expect_lt(max(Mod(coefficients@model$PSMS$f_bs - reference$field) /
      Mod(reference$field)), 1e-12)
    expect_identical(coefficients@model_parameters$PSMS$parameters$Amn_method,
      "Amn_fluid_coefficients")
    # Coefficient overlaps retain the same modes without angular quadrature ====
    without_quadrature <- evaluate(n_integration = 8L, adaptive = TRUE)
    expect_identical(without_quadrature@model$PSMS$f_bs, coefficients@model$PSMS$f_bs)
    expect_identical(without_quadrature@model_parameters$PSMS$parameters$acoustics,
      coefficients@model_parameters$PSMS$parameters$acoustics)
  }
})

test_that("quad liquid overlaps agree across aspect ratios and material contrasts", {
  skip_if_not(acousticTS:::.quad_precision_available())
  fixtures <- list(
    list(length = 0.04, radius = 0.015, density = 1000, speed = 1400, angle = 90),
    list(length = 0.3, radius = 0.004, density = 1200, speed = 1600, angle = 37)
  )
  for (fixture in fixtures) {
    object <- fls_generate(
      prolate_spheroid(length_body = fixture$length, radius_body = fixture$radius),
      density_body = fixture$density, sound_speed_body = fixture$speed,
      theta_body = fixture$angle * pi / 180)
    fit <- target_strength(object, c(12000, 38000), "psms",
      boundary = "liquid_filled", precision = "quad",
      density_sw = 1026.8, sound_speed_sw = 1477.3)
    inputs <- fit@model_parameters$PSMS
    expect_true(all(is.finite(fit@model$PSMS$f_bs)))
    evaluate_reference <- function(acoustics, order = 192L) {
      acousticTS:::.prolate_spheroidal_kernels_fixed(
        acoustics = acoustics, body = inputs$body, medium = inputs$medium,
        boundary_method = "Amn_fluid", n_integration = order,
        precision = "quad", adaptive = FALSE)
    }
    reference <- psms_enlarged_reference(inputs$parameters$acoustics,
      evaluate_reference)
    more_nodes <- evaluate_reference(reference$acoustics, 256L)
    expect_lt(max(Mod(more_nodes - reference$field) / Mod(more_nodes)), 1e-12)
    expect_lt(max(Mod(fit@model$PSMS$f_bs - reference$field) /
      Mod(reference$field)), 1e-12)
  }
})

test_that("full gas spheroids agree with the independent native reference", {
  object <- gas_generate(
    shape = prolate_spheroid(length_body = 0.2, radius_body = 0.01),
    density_fluid = 1.24, sound_speed_fluid = 345
  )
  # Prol_Spheroid vectorized 995929a with precision_fbs = 1e-8 =================
  reference <- c(-26.549891089759697, -24.53103140836401)
  precisions <- c("double", if (acousticTS:::.quad_precision_available()) "quad")
  for (precision in precisions) {
    selected <- 1:2
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
  # Quad M=64/80, N=320/400 and 512/640-node enlargements agree ================
  result <- target_strength(
    object, 200e3, "psms", boundary = "gas_filled", precision = "quad",
    density_sw = 1026.8, sound_speed_sw = 1477.3
  )
  expect_equal(result@model$PSMS$TS, -29.0909526869064, tolerance = 1e-10)

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

test_that("full fluid refinement accepts equal initial order and degree limits", {
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
      # Only automatically refined methods require enlarged fixed references ===
      if (method %in% c("Amn_fixed_rigid", "Amn_pressure_release", "Amn_fluid_gas")) {
        evaluate_reference <- function(limits, order = 128L) {
          acousticTS:::prolate_spheroid_tmatrix_cpp(
            limits, body, medium, gauss_legendre(order, -1, 1), "quad", method
          )$f_scat
        }
        reference <- psms_enlarged_reference(acoustics, evaluate_reference,
          axial = angle == 0)
        retained <- reference$field
        if (method == "Amn_fluid_gas") {
          more_nodes <- evaluate_reference(reference$acoustics, 192L)
          expect_lt(Mod(more_nodes - retained) / Mod(more_nodes), 1e-12)
        }
      } else {
        retained <- acousticTS:::prolate_spheroid_tmatrix_cpp(
          acoustics, body, medium, quadrature, "quad", method
        )$f_scat
      }
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
  # Resolve the larger transverse incident support before checking stability ===
  reference_limits <- acoustics
  reference_limits$m_max <- 36L
  reference_limits$n_max <- 51L
  reference <- psms_enlarged_reference(reference_limits, function(limits) {
    acousticTS:::prolate_spheroid_tmatrix_cpp(
      limits, body, medium, quadrature, "double", "Amn_pressure_release"
    )$f_scat
  }, tolerance = 1e-13)
  expect_equal(actual, reference$field, tolerance = 1e-13)
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
    gas_promoted <- acousticTS:::.psms_promote_full_fluid_modal_limits(
      model_params = acousticTS:::.psms_complete_acoustics(
        model_params = gas_model_params,
        body_params = body_params,
        scatterer_shape = shape
      ),
      body_params = body_params,
      boundary = "gas_filled"
    )
    liquid_unchanged <- acousticTS:::.psms_promote_full_fluid_modal_limits(
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
