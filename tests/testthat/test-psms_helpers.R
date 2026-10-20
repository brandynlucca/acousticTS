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
    "not representable by finite spheroidal coordinates")
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

# Retained T-matrix evaluation supplies an independent fixed-support reference =
psms_fixed_reference <- function(acoustics, body, medium, boundary_method,
                                 n_integration, precision, adaptive = FALSE) {
  stopifnot(!adaptive)
  acousticTS:::prolate_spheroid_tmatrix_cpp(acoustics, body, medium,
    gauss_legendre(n_integration, -1, 1), precision, boundary_method)$f_scat
}

# Converge liquid references independently of adaptive stopping ================
psms_liquid_enlarged_reference <- function(acoustics, body, medium) {
  psms_enlarged_reference(acoustics, function(retained) {
    acousticTS:::prolate_spheroid_fbs(retained, body, medium,
      list(nodes = numeric(), weights = numeric()), precision = "quad",
      Amn_method = "Amn_fluid_coefficients", adaptive = FALSE)
  }, axial = body$theta_body %in% c(0, pi))$field
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

test_that("PSMS preserves liquid quadrature and gas coefficient overlaps", {
  object <- fls_generate(prolate_spheroid(length_body = 0.14, radius_body = 0.01),
    density_body = 1028.9, sound_speed_body = 1480.3)
  initialize <- function(...) acousticTS:::psms_initialize(object, 70000, ...)
  liquid <- initialize()
  expect_identical(liquid@model_parameters$PSMS$parameters$Amn_method,
    "Amn_fluid")
  gas <- initialize(boundary = "gas_filled")
  expect_identical(gas@model_parameters$PSMS$parameters$Amn_method, "Amn_fluid_gas")
  simplified <- initialize(simplify_Amn = TRUE)
  expect_identical(simplified@model_parameters$PSMS$parameters$Amn_method,
    "Amn_fluid_simplify")
  expect_true(acousticTS:::.psms_is_full_fluid_method("Amn_fluid_coefficients"))
  expect_identical(acousticTS:::.psms_Amn_method("liquid_filled", 1),
    "Amn_fluid_simplify")
  expect_identical(acousticTS:::.psms_Amn_method("liquid_filled", 0),
    "Amn_fluid")
  expect_identical(acousticTS:::.psms_Amn_method("liquid_filled", NA),
    NA)
  if (acousticTS:::.quad_precision_available()) {
    coefficients <- initialize(precision = "quad")
    expect_identical(coefficients@model_parameters$PSMS$parameters$Amn_method,
      "Amn_fluid")
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
      field <- psms_fixed_reference(
        acoustics = enlarged, body = inputs$body, medium = inputs$medium,
        boundary_method = inputs$parameters$Amn_method,
        n_integration = 96L, precision = precision, adaptive = FALSE)
      expect_lt(Mod(result@model$PSMS$f_bs - field) / Mod(field), 1e-10)
    }
  }
})

test_that("automatic liquid integration resolves the high-size overlap integrals", {
  skip_if_not(acousticTS:::.quad_precision_available())
  cases <- data.frame(h0 = c(90, 180, 120), angle = c(0, 0, 37),
    points = c(192L, 384L, 192L))
  for (index in seq_len(nrow(cases))) {
    selected <- cases[index, ]
    object <- fls_generate(
      prolate_spheroid(length_body = 0.14, radius_body = 0.01),
      density_body = 1028.9, sound_speed_body = 1480.3,
      theta_body = selected$angle * pi / 180)
    frequency <- selected$h0 * 1477.3 / (2 * pi * sqrt(0.07^2 - 0.01^2))
    evaluate <- function(adaptive, points = NULL) {
      target_strength(object, frequency, "psms", boundary = "liquid_filled",
        precision = "quad", adaptive = adaptive, n_integration = points,
        density_sw = 1026.8, sound_speed_sw = 1477.3)
    }
    reference <- evaluate(FALSE, selected$points)
    fixed <- evaluate(FALSE)
    adapted <- evaluate(TRUE)
    expect_identical(reference@model_parameters$PSMS$parameters$n_integration,
      selected$points)
    expect_true(is.na(fixed@model_parameters$PSMS$parameters$n_integration))
    expect_equal(fixed@model_parameters$PSMS$parameters$acoustics$n_max,
      reference@model_parameters$PSMS$parameters$acoustics$n_max)
    expect_lt(Mod(fixed@model$PSMS$f_bs - reference@model$PSMS$f_bs) /
      Mod(reference@model$PSMS$f_bs), 1e-12)
    # Check adaptive TS against independently enlarged modal support ===========
    inputs <- adapted@model_parameters$PSMS
    enlarged <- psms_liquid_enlarged_reference(inputs$parameters$acoustics,
      inputs$body, inputs$medium)
    expect_lte(max(abs(20 * log10(Mod(adapted@model$PSMS$f_bs) /
      Mod(enlarged)))), 1e-8)
  }
})

test_that("quad liquid quadrature and adaptive controls retain numerical agreement", {
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
    fixed <- evaluate()
    inputs <- fixed@model_parameters$PSMS
    evaluate_reference <- function(acoustics, order = 96L) {
      psms_fixed_reference(
        acoustics = acoustics, body = inputs$body, medium = inputs$medium,
        boundary_method = "Amn_fluid", n_integration = order,
        precision = "quad", adaptive = FALSE)
    }
    # Separate quadrature stability from agreement at prescribed modal limits ==
    reference <- evaluate_reference(inputs$parameters$acoustics)
    more_nodes <- evaluate_reference(inputs$parameters$acoustics, 128L)
    expect_lt(max(Mod(more_nodes - reference) / Mod(more_nodes)), 1e-12)
    expect_lt(max(Mod(fixed@model$PSMS$f_bs - reference) /
      Mod(reference)), 1e-12)
    expect_identical(fixed@model_parameters$PSMS$parameters$Amn_method,
      "Amn_fluid")
    # Explicit counts remain effective when adaptive quadrature is disabled ====
    denser <- evaluate(n_integration = 128L)
    expect_identical(denser@model_parameters$PSMS$parameters$n_integration, 128L)
    expect_lt(max(Mod(denser@model$PSMS$f_bs - fixed@model$PSMS$f_bs) /
      Mod(fixed@model$PSMS$f_bs)), 1e-12)
    # Adapt quadrature and outgoing tails while retaining full coupling ========
    expect_warning(adapted <- evaluate(n_integration = 8L, adaptive = TRUE),
      "'n_integration' is ignored when 'adaptive = TRUE'")
    expect_true(is.na(adapted@model_parameters$PSMS$parameters$n_integration))
    # Check adaptive TS against independently enlarged modal support ===========
    enlarged <- psms_liquid_enlarged_reference(inputs$parameters$acoustics,
      inputs$body, inputs$medium)
    expect_lte(max(abs(20 * log10(Mod(adapted@model$PSMS$f_bs) /
      Mod(enlarged)))), 1e-8)
    # Never reuse incomplete factors left by an adaptive solve =================
    repeated <- evaluate()
    expect_identical(repeated@model$PSMS$f_bs, fixed@model$PSMS$f_bs)
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
    evaluate_reference <- function(acoustics, order = 96L) {
      psms_fixed_reference(
        acoustics = acoustics, body = inputs$body, medium = inputs$medium,
        boundary_method = "Amn_fluid", n_integration = order,
        precision = "quad", adaptive = FALSE)
    }
    # Match modal support while checking quadrature stability independently ====
    reference <- evaluate_reference(inputs$parameters$acoustics)
    more_nodes <- evaluate_reference(inputs$parameters$acoustics, 128L)
    expect_lt(max(Mod(more_nodes - reference) / Mod(more_nodes)), 1e-12)
    expect_lt(max(Mod(fit@model$PSMS$f_bs - reference) /
      Mod(reference)), 1e-12)
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
  initialized <- acousticTS:::psms_initialize(
    object, 200e3, boundary = "gas_filled", precision = "quad",
    density_sw = 1026.8, sound_speed_sw = 1477.3
  )
  # Evaluate the independently verified converged support explicitly ===========
  inputs <- initialized@model_parameters$PSMS
  acoustics <- inputs$parameters$acoustics
  acoustics$m_max <- 64L
  acoustics$n_max <- 320L
  field <- acousticTS:::prolate_spheroid_fbs(
    acoustics, inputs$body, inputs$medium,
    list(nodes = numeric(), weights = numeric()), "quad", "Amn_fluid_gas",
    adaptive = FALSE
  )
  target <- 20 * log10(Mod(-2i / acoustics$k_sw * field))
  expect_equal(target, -29.0909526869064, tolerance = 1e-10)

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


test_that("quad full gas recovers the converged axial field at h0 = 120", {
  skip_if_not(acousticTS:::.quad_precision_available())
  semimajor <- 0.07
  semiminor <- 0.01
  frequency <- 120 * 1477.3 / (2 * pi * sqrt(semimajor^2 - semiminor^2))
  wave_number <- 2 * pi * frequency / 1477.3
  object <- gas_generate(
    prolate_spheroid(length_body = 2 * semimajor, radius_body = semiminor),
    density_fluid = 1.24, sound_speed_fluid = 345, theta_body = 0
  )
  # Native quad N=621 and N=776 agree in the complete complex field ============
  reference_field <- complex(
    real = 0.182688116344343582053401185024321438,
    imaginary = 0.657194044775623801719392956161668946
  )
  reference_amplitude <- -2i / wave_number * reference_field
  for (adaptive in c(FALSE, TRUE)) {
    result <- target_strength(
      object, frequency, "psms", boundary = "gas_filled", precision = "quad",
      simplify_Amn = FALSE, adaptive = adaptive,
      density_sw = 1026.8, sound_speed_sw = 1477.3
    )
    # Convert the returned modal sum to the physical scattering amplitude ======
    amplitude <- -2i / wave_number * result@model$PSMS$f_bs
    expect_length(amplitude, 1L)
    expect_true(all(is.finite(amplitude)))
    expect_lt(Mod(amplitude - reference_amplitude) / Mod(reference_amplitude),
      1e-12)
  }
})

test_that("full fluid solves retain equal prescribed order and degree limits", {
  skip_if_not(acousticTS:::.quad_precision_available())
  for (boundary in c("liquid_filled", "gas_filled")) {
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
    # Match the prescribed support in the independent retained T-matrix ========
    reference <- psms_fixed_reference(
      acoustics, model$body, model$medium, model$parameters$Amn_method,
      n_integration = if (boundary == "gas_filled") 192L else 96L,
      precision = "quad"
    )
    if (boundary == "gas_filled") {
      denser <- psms_fixed_reference(
        acoustics, model$body, model$medium, model$parameters$Amn_method,
        n_integration = 256L, precision = "quad"
      )
      expect_lt(Mod(denser - reference) / Mod(denser), 1e-12)
    }
    actual <- 20 * log10(Mod(-2i / acoustics$k_sw * field))
    expected <- 20 * log10(Mod(-2i / acoustics$k_sw * reference))
    expect_equal(as.numeric(actual), as.numeric(expected), tolerance = 2e-9)
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
      # Compare the two implementations at identical prescribed modal limits ===
      retained <- acousticTS:::prolate_spheroid_tmatrix_cpp(
        acoustics, body, medium, quadrature, "quad", method
      )$f_scat
      if (method == "Amn_fluid_gas") {
        more_nodes <- acousticTS:::prolate_spheroid_tmatrix_cpp(
          acoustics, body, medium, gauss_legendre(192L, -1, 1), "quad", method
        )$f_scat
        expect_lt(Mod(more_nodes - retained) / Mod(more_nodes), 1e-12)
      }
      for (adaptive in c(FALSE, TRUE)) {
        actual <- acousticTS:::prolate_spheroid_fbs(
          acoustics, body, medium, quadrature, "quad", method,
          adaptive = adaptive
        )
        if (method == "Amn_fluid" && adaptive) {
          # Adaptive liquid selects its support against the TS target ==========
          enlarged <- psms_liquid_enlarged_reference(acoustics, body, medium)
          expect_lte(abs(20 * log10(Mod(actual) / Mod(enlarged))), 1e-8)
        } else {
          expect_equal(actual, retained, tolerance = 1e-7,
            info = paste(method, angle, adaptive))
        }
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
  # Hold modal support fixed to isolate the specified receiver direction =======
  reference <- acousticTS:::prolate_spheroid_tmatrix_cpp(
    acoustics, body, medium, quadrature, "double", "Amn_pressure_release"
  )$f_scat
  expect_equal(actual, reference, tolerance = 1e-13)
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

test_that("PSMS rejects unrepresentable fixed quadrature counts", {
  object <- fls_generate(prolate_spheroid(length_body = 0.14, radius_body = 0.01),
    density_body = 1028.9, sound_speed_body = 1480.3)
  for (points in c(Inf, -Inf, NaN, NA_real_, 0, 1.5,
    as.double(.Machine$integer.max) + 1, 1e100)) {
    expect_error(acousticTS:::psms_initialize(object, 12000,
      boundary = "liquid_filled", adaptive = FALSE, n_integration = points),
      "'n_integration' must be a single positive integer within R's integer range")
  }
})

test_that("liquid quadrature controls reject insufficient fixed integration", {
  object <- fls_generate(prolate_spheroid(semimajor_length = 0.07,
    semiminor_length = 0.01), density_body = 1028.9,
    sound_speed_body = 1480.3, theta_body = 37 * pi / 180)
  precisions <- c("double", if (acousticTS:::.quad_precision_available()) "quad")
  for (precision in precisions) {
    evaluate <- function(points, adaptive) {
      target_strength(object, 12000, "psms", boundary = "liquid_filled",
        precision = precision, n_integration = points, adaptive = adaptive,
        density_sw = 1026.8, sound_speed_sw = 1477.3)
    }
    expect_error(evaluate(1L, FALSE), "[Qq]uadrature")
    fixed <- evaluate(96L, FALSE)
    # Compare integration methods without adaptive modal-tail truncation =======
    inputs <- fixed@model_parameters$PSMS
    coefficients <- acousticTS:::prolate_spheroid_fbs(
      inputs$parameters$acoustics, inputs$body, inputs$medium,
      list(nodes = numeric(), weights = numeric()), precision,
      "Amn_fluid_coefficients", adaptive = FALSE
    )
    expect_lt(Mod(coefficients - fixed@model$PSMS$f_bs) /
      Mod(fixed@model$PSMS$f_bs), if (precision == "quad") 1e-12 else 1e-8)
    # Requested quadrature counts must not affect adaptive integration =========
    expect_warning(adapted <- evaluate(1L, TRUE),
      "'n_integration' is ignored when 'adaptive = TRUE'")
    expect_warning(adapted_again <- evaluate(96L, TRUE),
      "'n_integration' is ignored when 'adaptive = TRUE'")
    expect_true(all(is.finite(adapted@model$PSMS$f_bs)))
    expect_true(is.na(adapted@model_parameters$PSMS$parameters$n_integration))
    expect_true(is.na(adapted_again@model_parameters$PSMS$parameters$n_integration))
    expect_identical(adapted@model$PSMS$f_bs, adapted_again@model$PSMS$f_bs)
  }
})


test_that("adaptive quad liquid support meets TS and complex targets", {
  skip_if_not(acousticTS:::.quad_precision_available())
  semimajor <- 0.01
  semiminor <- 0.75 * semimajor
  sound_speed <- 1477.3
  frequency <- 5 * sound_speed / (2 * pi * sqrt(semimajor^2 - semiminor^2))
  object <- fls_generate(prolate_spheroid(semimajor_length = semimajor,
    semiminor_length = semiminor), density_body = 1028.9,
    sound_speed_body = 1480.3, theta_body = 89.999 * pi / 180)
  optimized <- target_strength(object, frequency, "psms",
    boundary = "liquid_filled", precision = "quad", adaptive = TRUE,
    density_sw = 1026.8, sound_speed_sw = sound_speed)
  inputs <- optimized@model_parameters$PSMS
  # Enlarge both support directions without adaptive stopping ==================
  enlarged <- inputs$parameters$acoustics
  enlarged$m_max <- pmax(enlarged$m_max + 8L, 24L)
  enlarged$n_max <- pmax(enlarged$n_max + 16L, 36L, enlarged$m_max)
  reference <- acousticTS:::prolate_spheroid_fbs(
    enlarged, inputs$body, inputs$medium,
    list(nodes = numeric(), weights = numeric()), precision = "quad",
    Amn_method = "Amn_fluid_coefficients", adaptive = FALSE)
  reference_ts <- 20 * log10(Mod(-2i / enlarged$k_sw * reference))
  expect_true(all(is.finite(optimized@model$PSMS$TS)))
  expect_lte(max(abs(optimized@model$PSMS$TS - reference_ts)), 1e-8)
  optimized_field <- optimized@model$PSMS$f_bs
  expect_true(all(Mod(optimized_field - reference) <= 1e-18 +
    1e-14 * Mod(reference)))
})

test_that("adaptive quad liquid returns the exact matched-medium solution", {
  skip_if_not(acousticTS:::.quad_precision_available())
  object <- fls_generate(prolate_spheroid(semimajor_length = 0.01,
    semiminor_length = 0.0075), density_body = 1026.8,
    sound_speed_body = 1477.3, theta_body = 0.37)
  matched <- target_strength(object, 38000, "psms",
    boundary = "liquid_filled", precision = "quad", adaptive = TRUE,
    density_sw = 1026.8, sound_speed_sw = 1477.3)
  expect_true(all(matched@model$PSMS$f_bs == 0i))
  expect_identical(matched@model$PSMS$TS, -Inf)
})


test_that("adaptive quad hard boundaries match enlarged-support references", {
  skip_if_not(acousticTS:::.quad_precision_available())
  # Quad references use independently enlarged degree and order support ========
  cases <- data.frame(
    boundary = rep(c("fixed_rigid", "pressure_release"), 4L),
    semimajor_m = rep(c(0.001, 0.01, 0.001, 0.01), each = 2L),
    aspect_ratio = rep(c(100, 5, 1 / 0.9, 1 / 0.75),
      each = 2L),
    incidence_deg = rep(c(0, 45, 0.0001, 89.999), each = 2L),
    frequency = rep(c(0x1.ae89d283d4291p+24,
      0x1.d4aff795ee44ap+18,
      0x1.07612b3284a0ep+19,
      0x1.5b22dad00543ep+17), each = 2L),
    reference_real = c(6.37998830762401526607560478301551386e-06,
      -0.00231054766363223068713298784191789936,
      -0.270303041864811095929197098942407886,
      0.330896130822711613907144017486289914,
      -0.359307021100057917449804245034194428,
      0.483552953730301897755781437363595613,
      -1.57653241517880346146415797755356448,
      1.82190537750399140983966870507909978),
    reference_imaginary = c(0.000783301731615198875682995259801178308,
      -0.0633636064825864037530188103511905989,
      0.0565510431422121736774288513736311321,
      0.364658434424076755911778314881162379,
      0.250208551116167601850576350633109049,
      0.174753789240174688153480489535039548,
      1.07772791206368720569079611274673828,
      -0.557421921909752217543144402163310505)
  )
  for (index in seq_len(nrow(cases))) {
    case <- cases[index, ]
    object <- fls_generate(prolate_spheroid(length_body = 2 * case$semimajor_m,
      radius_body = case$semimajor_m / case$aspect_ratio),
      density_body = 1028.9, sound_speed_body = 1480.3,
      theta_body = case$incidence_deg * pi / 180)
    optimized <- target_strength(object, case$frequency, "psms",
      boundary = case$boundary, precision = "quad", adaptive = TRUE,
      density_sw = 1026.8, sound_speed_sw = 1477.3, phi_body = pi)
    expected <- complex(real = case$reference_real,
      imaginary = case$reference_imaginary)
    actual <- optimized@model$PSMS$f_bs
    wavenumber <- optimized@model_parameters$PSMS$parameters$acoustics$k_sw
    expect_true(all(is.finite(actual)))
    expect_lte(Mod(actual - expected), 1e-18 + 1e-14 * Mod(expected),
      label = paste(case$boundary, "field error at", case$incidence_deg, "degrees"))
    expect_lte(abs(optimized@model$PSMS$TS -
      20 * log10(Mod(-2i / wavenumber * expected))), 1e-8)
  }
})

test_that("adaptive quad PSMS preserves frequency vectors and scalar results", {
  skip_if_not(acousticTS:::.quad_precision_available())
  object <- gas_generate(
    prolate_spheroid(length_body = 0.02, radius_body = 0.002),
    density_fluid = 1.24, sound_speed_fluid = 345,
    theta_body = 37.19 * pi / 180
  )
  original <- object
  frequency <- c(2, 1, 2) * 1477.3 /
    (2 * pi * sqrt(0.01^2 - 0.002^2))
  for (boundary in c("gas_filled", "fixed_rigid", "pressure_release")) {
    evaluate <- function(values) {
      target_strength(object, values, "psms", boundary = boundary,
        adaptive = TRUE, precision = "quad", simplify_Amn = FALSE,
        density_sw = 1026.8, sound_speed_sw = 1477.3)
    }
    vector_result <- evaluate(frequency)
    scalar_results <- lapply(frequency, evaluate)
    expected_field <- vapply(scalar_results,
      function(result) result@model$PSMS$f_bs, complex(1))
    expected_ts <- vapply(scalar_results,
      function(result) result@model$PSMS$TS, numeric(1))
    expect_identical(vector_result@model$PSMS$frequency, frequency)
    expect_true(all(Mod(vector_result@model$PSMS$f_bs - expected_field) <=
      1e-18 + 1e-14 * Mod(expected_field)))
    expect_lte(max(abs(vector_result@model$PSMS$TS - expected_ts)), 1e-8)
    expect_identical(vector_result@model_parameters$PSMS$parameters$adaptive,
      TRUE)
    expect_identical(vector_result@model$PSMS$f_bs[1],
      vector_result@model$PSMS$f_bs[3])
  }
  expect_identical(object, original)
})

test_that("adaptive quad PSMS preserves simulation grids and generated inputs", {
  skip_if_not(acousticTS:::.quad_precision_available())
  withr::local_seed(1379)
  object <- gas_generate(
    prolate_spheroid(length_body = 0.02, radius_body = 0.002),
    density_fluid = 1.24, sound_speed_fluid = 345
  )
  original <- object
  angle_draws <- 0L
  speed_draws <- 0L
  frequency <- c(30000, 15000)
  parameters <- list(
    boundary = c("gas_filled", "pressure_release"),
    length_body = c(0.02, 0.03),
    radius_body = 0.002,
    theta = function() {
      angle_draws <<- angle_draws + 1L
      runif(1, 0.1, 1.3)
    },
    sound_speed_sw = function() {
      speed_draws <<- speed_draws + 1L
      runif(1, 1470, 1485)
    },
    density_sw = 1026.8,
    adaptive = TRUE,
    precision = "quad",
    simplify_Amn = FALSE
  )
  actual <- simulate_ts(object, frequency, "psms", parameters = parameters,
    parallel = FALSE, verbose = FALSE)$PSMS
  expect_identical(angle_draws, 4L)
  expect_identical(speed_draws, 4L)
  expect_identical(actual$frequency, rep(frequency, 4L))
  expect_equal(nrow(actual), 8L)
  expected_grid <- expand.grid(boundary = parameters$boundary,
    length_body = parameters$length_body, stringsAsFactors = FALSE)
  cells <- seq.int(1L, nrow(actual), by = length(frequency))
  expect_identical(actual$boundary[cells], expected_grid$boundary)
  expect_identical(actual$length_body[cells], expected_grid$length_body)
  expect_identical(actual$realization, rep(1L, 8L))
  expect_identical(actual$radius_body, rep(0.002, 8L))
  expect_true(all(actual$adaptive))
  expect_true(all(actual$precision == "quad"))
  expect_length(unique(actual$theta), 4L)
  expect_length(unique(actual$sound_speed_sw), 4L)
  # Reconstruct each resolved draw through independent public scalar calls =====
  for (index in seq_len(nrow(actual))) {
    row <- actual[index, ]
    reference_object <- gas_generate(
      prolate_spheroid(length_body = row$length_body,
        radius_body = row$radius_body),
      density_fluid = 1.24, sound_speed_fluid = 345,
      theta_body = row$theta
    )
    reference <- target_strength(reference_object, row$frequency, "psms",
      boundary = row$boundary, adaptive = TRUE, precision = "quad",
      simplify_Amn = FALSE, density_sw = row$density_sw,
      sound_speed_sw = row$sound_speed_sw)@model$PSMS
    expect_lte(Mod(row$f_bs - reference$f_bs),
      1e-18 + 1e-14 * Mod(reference$f_bs))
    expect_lte(abs(row$TS - reference$TS), 1e-8)
  }
  expect_identical(object, original)
})

test_that("adaptive quad gas agrees with independently enlarged references", {
  skip_if_not(acousticTS:::.quad_precision_available())
  # Two successive joint enlargements confirmed each native quad reference =====
  # Final (M,N,Ni): (11,66,66), (11,66,66), (16,241,241), (21,218,218) =========
  # The final axial Robin reference uses (M,N,Ni) = (0,931,1214) ===============
  cases <- data.frame(
    semimajor_m = c(0.001, 0.001, 0.001, 0.01, 0.001),
    frequency = c(0x1.1f068c57e2c6p+20, 0x1.1f068c57e2c6p+20,
      0x1.1f068c57e2c6p+22, 0x1.cb3dad596ad66p+18,
      0x1.ae89d283d4291p+24),
    incidence = c(0, 0x1.0c152382d7365p-1, 0x1.0c152382d7365p-1,
      0x1.921fb54442d18p+0, 0),
    chi_sw = c(0x1.4p+2, 0x1.4p+2, 0x1.4p+4,
      0x1.4000000000001p+4, 0x1.e000000000001p+6),
    chi_body = c(0x1.568ff420a636fp+4, 0x1.568ff420a636fp+4,
      0x1.568ff420a636fp+6, 0x1.568ff420a636fp+6,
      0x1.00ebf7187ca93p+9),
    reference_real = c(-0.0149132796747873439902380674882470732,
      -0.0330585259897998436184308747532578063,
      -0.0292157432844645998068917828336100188,
      -2.45013905130522499746404858835601201,
      -0.00321730834899562075031858529974459608),
    reference_imaginary = c(0.028068678398639622227980460933428143,
      -0.0465735001884262774735069284435437581,
      0.0128533562928119153046266375917432826,
      -2.81634227872851779357279836203608593,
      -0.0637679656318616669605227867913564312)
  )
  for (index in seq_len(nrow(cases))) {
    case <- cases[index, ]
    object <- gas_generate(
      prolate_spheroid(length_body = 2 * case$semimajor_m,
        radius_body = case$semimajor_m / 100),
      density_fluid = 1.24, sound_speed_fluid = 345,
      theta_body = case$incidence
    )
    optimized <- target_strength(object, case$frequency, "psms",
      boundary = "gas_filled", adaptive = TRUE, precision = "quad",
      simplify_Amn = FALSE, density_sw = 1026.8,
      sound_speed_sw = 1477.3, phi_body = pi)
    inputs <- optimized@model_parameters$PSMS
    expect_identical(inputs$parameters$acoustics$chi_sw, case$chi_sw)
    expect_identical(inputs$parameters$acoustics$chi_body, case$chi_body)
    expect_identical(inputs$body$xi, 0x1.000346ec78e6bp+0)
    reference <- complex(real = case$reference_real,
      imaginary = case$reference_imaginary)
    field <- optimized@model$PSMS$f_bs
    wavenumber <- inputs$parameters$acoustics$k_sw
    expect_true(all(is.finite(field)))
    expect_lte(Mod(field - reference), 1e-18 + 1e-14 * Mod(reference))
    expect_lte(abs(optimized@model$PSMS$TS -
      20 * log10(Mod(-2i / wavenumber * reference))), 1e-8)
  }
})
