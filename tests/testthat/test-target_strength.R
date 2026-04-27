library(acousticTS)

test_that("target_strength function works correctly", {
  # Test target_strength with calibration sphere
  cal_obj <- cal_generate()
  frequency <- c(38e3, 120e3)

  # Calculate target strength
  cal_with_ts <- target_strength(
    object = cal_obj,
    frequency = frequency,
    model = "calibration"
  )

  # Check that object is still CAL class
  expect_s4_class(cal_with_ts, "CAL")

  # Check that model slot was populated
  expect_true("calibration" %in% names(cal_with_ts@model))
  expect_s3_class(cal_with_ts@model$calibration, "data.frame")

  # Check that frequency was set correctly
  expect_equal(cal_with_ts@model$calibration$frequency, frequency)

  # Check that TS values were calculated (should be numeric, not NA)
  expect_type(cal_with_ts@model$calibration$TS, "double")
  expect_true(all(!is.na(cal_with_ts@model$calibration$TS)))
})

test_that("calibration sphere modal sum converges at high ka", {
  cal_obj <- cal_generate(material = "WC", diameter = 38.1e-3)

  cal_with_ts <- target_strength(
    object = cal_obj,
    frequency = 350e3,
    model = "calibration",
    sound_speed_sw = 1477.3,
    density_sw = 1026.8
  )

  expect_equal(
    cal_with_ts@model$calibration$TS,
    -44.2348898912626,
    tolerance = 1e-10
  )
})

test_that("calibration sphere can use the legacy fixed modal cutoff", {
  cal_obj <- cal_generate(material = "WC", diameter = 38.1e-3)

  cal_with_ts <- target_strength(
    object = cal_obj,
    frequency = 350e3,
    model = "calibration",
    sound_speed_sw = 1477.3,
    density_sw = 1026.8,
    adaptive = FALSE
  )

  expect_equal(
    cal_with_ts@model$calibration$TS,
    -44.2348180844783,
    tolerance = 1e-10
  )
})

test_that("target_strength works with different scatterer types", {
  # Test with FLS object and DWBA/SDWBA/KRM model
  data(krill, package = "acousticTS")
  data(sardine, package = "acousticTS")
  data(cod, package = "acousticTS")

  fls_obj <- fls_generate(
    x = krill@body$rpos[1, ],
    y = krill@body$rpos[2, ],
    z = krill@body$rpos[3, ],
    radius_body = krill@body$radius,
    g_body = krill@body$g,
    h_body = krill@body$h,
    radius_curvature_ratio = 3.0,
    theta_body = krill@body$theta
  )
  frequency <- seq(1e3, 10e3, 1e3)

  expect_warning(
    fls_with_ts <- target_strength(
      object = fls_obj,
      frequency = frequency,
      model = c("SDWBA", "SDWBA_curved", "DWBA", "KRM")
    ),
    "deprecated"
  )

  expect_s4_class(fls_with_ts, "FLS")
  expect_true(
    all(c("DWBA", "SDWBA", "SDWBA_curved", "KRM")
    %in% names(fls_with_ts@model))
  )

  # Test SBF object with KRM
  sardine_with_ts <- target_strength(
    object = sardine,
    frequency = frequency,
    model = "KRM"
  )

  expect_s4_class(sardine_with_ts, "SBF")
  expect_true("KRM" %in% names(sardine_with_ts@model))

  cod_with_ts <- target_strength(
    object = cod,
    frequency = frequency,
    model = "KRM"
  )

  expect_s4_class(cod_with_ts, "SBF")
  expect_true("KRM" %in% names(cod_with_ts@model))
  expect_equal(
    cod_with_ts@model_parameters$KRM$parameters$ns_b,
    ncol(cod@body$rpos)
  )
  expect_equal(
    cod_with_ts@model_parameters$KRM$parameters$ns_sb,
    ncol(cod@bladder$rpos)
  )

  # Test with GAS object and MSS Anderson (1950) model
  gas_obj <- gas_generate(shape = sphere(radius_body = 1e-3, n_segments = 80))

  gas_with_ts <- target_strength(
    object = gas_obj,
    frequency = frequency,
    model = "sphms"
  )

  expect_s4_class(gas_with_ts, "GAS")
  expect_true("SPHMS" %in% names(gas_with_ts@model))

  # Test ESS object with HP and MSS model
  ess_obj <- ess_generate(
    shape = sphere(radius_body = 10e-3, n_segments = 80),
    radius_shell = 10e-3, # 10 mm outer radius
    shell_thickness = 1e-3, # 1 mm shell thickness
    sound_speed_shell = 3750, # Shell sound speed (m/s)
    sound_speed_fluid = 1575, # Internal fluid sound speed (m/s)
    density_shell = 2565, # Shell density (kg/m³)
    density_fluid = 1077.3, # Internal fluid density (kg/m³)
    K = 70e9, # Bulk modulus (Pa)
    nu = 0.32 # Poisson's ratio
  )

  ess_with_ts <- target_strength(
    object = ess_obj,
    frequency = frequency,
    model = c("hpa", "essms")
  )

  expect_s4_class(ess_with_ts, "ESS")
  expect_true(
    all(c("ESSMS", "HPA")
    %in% names(ess_with_ts@model))
  )
})

test_that("target_strength can update existing models", {
  # Test that target_strength can be called multiple times
  cal_obj <- cal_generate()
  frequency1 <- c(38e3)
  frequency2 <- c(120e3, 200e3)

  # First calculation
  cal_ts1 <- target_strength(cal_obj, frequency1, "calibration")
  expect_equal(nrow(cal_ts1@model$calibration), length(frequency1))

  # Second calculation (should update/replace)
  cal_ts2 <- target_strength(cal_ts1, frequency2, "calibration")
  expect_equal(nrow(cal_ts2@model$calibration), length(frequency2))
  expect_equal(cal_ts2@model$calibration$frequency, frequency2)
})

test_that("target_strength handles multiple models", {
  # Test with FLS that can use multiple models
  data(krill, package = "acousticTS")
  fls_obj <- fls_generate(
    x = krill@body$rpos[1, ],
    y = krill@body$rpos[2, ],
    z = krill@body$rpos[3, ],
    radius_body = krill@body$radius,
    g_body = krill@body$g,
    h_body = krill@body$h,
    radius_curvature_ratio = 3.3,
    theta_body = krill@body$theta
  )
  frequency <- c(120e3)

  # Calculate with DWBA
  fls_dwba <- target_strength(fls_obj, frequency, "DWBA")

  # Should have the DWBA
  expect_true("DWBA" %in% names(fls_dwba@model))
})

test_that("target_strength applies shared and model-specific arguments cleanly", {
  data(krill, package = "acousticTS")

  fls_obj <- fls_generate(
    x = krill@body$rpos[1, ],
    y = krill@body$rpos[2, ],
    z = krill@body$rpos[3, ],
    radius_body = krill@body$radius,
    g_body = krill@body$g,
    h_body = krill@body$h,
    radius_curvature_ratio = 3.3,
    theta_body = krill@body$theta
  )

  out <- target_strength(
    object = fls_obj,
    frequency = 120e3,
    model = c("dwba", "sdwba"),
    density_sw = 1026,
    sound_speed_sw = 1478,
    model_args = list(
      sdwba = c(
        n_iterations = 5,
        n_segments_init = 14,
        phase_sd_init = 0.77,
        length_init = 38.35e-3,
        frequency_init = 120e3
      )
    )
  )

  expect_true(all(c("DWBA", "SDWBA") %in% names(out@model)))
  expect_equal(out@model_parameters$DWBA$medium$density, 1026)
  expect_equal(out@model_parameters$SDWBA$medium$density, 1026)
  expect_equal(
    out@model_parameters$SDWBA$parameters[[1]]$meta_params$p0,
    0.77
  )
  expect_equal(
    out@model_parameters$SDWBA$parameters[[1]]$meta_params$n_iterations,
    5
  )
})

test_that("target_strength lets model_args override shared arguments", {
  data(krill, package = "acousticTS")

  fls_obj <- fls_generate(
    x = krill@body$rpos[1, ],
    y = krill@body$rpos[2, ],
    z = krill@body$rpos[3, ],
    radius_body = krill@body$radius,
    g_body = krill@body$g,
    h_body = krill@body$h,
    radius_curvature_ratio = 3.3,
    theta_body = krill@body$theta
  )

  out <- target_strength(
    object = fls_obj,
    frequency = 120e3,
    model = c("dwba", "sdwba"),
    density_sw = 1026,
    sound_speed_sw = 1478,
    model_args = list(
      SDWBA = list(
        density_sw = 1027,
        n_iterations = 5,
        n_segments_init = 14,
        phase_sd_init = 0.77,
        length_init = 38.35e-3,
        frequency_init = 120e3
      )
    )
  )

  expect_equal(out@model_parameters$DWBA$medium$density, 1026)
  expect_equal(out@model_parameters$SDWBA$medium$density, 1027)
})

test_that("target_strength preserves object metadata", {
  # Test that metadata is preserved after target_strength calculation
  cal_obj <- cal_generate(ID = "TestSphere")
  original_metadata <- cal_obj@metadata

  cal_with_ts <- target_strength(cal_obj, 38e3, "calibration")

  # Metadata should be unchanged
  expect_equal(cal_with_ts@metadata, original_metadata)
  expect_equal(cal_with_ts@metadata$ID, "TestSphere")
})


test_that("Expected error states for target_strength", {
  data(krill, package = "acousticTS")

  # Test missing 'object'
  expect_error(
    target_strength(frequency = 120e3, model = "DWBA"),
    "Scattering object \\(\\'object'\\)\\ is required"
  )

  # Test missing 'frequency'
  expect_error(
    target_strength(object = krill, model = "DWBA"),
    "Frequency \\(\\Hz\\)\\ \\(\\'frequency'\\)\\ is required"
  )

  # Test missing 'model'
  expect_error(
    target_strength(object = krill, frequency = 120e3),
    "Target strength model \\(\\'model'\\)\\ is required"
  )

  # Test invalid 'model'
  expect_error(
    target_strength(object = krill, frequency = 120e3, model = "invalid"),
    "Unknown target strength model 'invalid'"
  )

  # Test missing initialize-model function pairing
  expect_error(
    target_strength(object = krill, frequency = 120e3, model = "spoof"),
    "Unknown target strength model 'spoof'"
  )

  expect_error(
    target_strength(
      object = krill,
      frequency = 120e3,
      model = "dwba",
      model_args = list(sdwba = list(phase_sd_init = 0.77))
    ),
    "'model_args' contains entries for model\\(s\\) not requested in 'model'"
  )

  expect_error(
    target_strength(
      object = krill,
      frequency = 120e3,
      model = "dwba",
      model_args = list(dwba = 1:3)
    ),
    "must be either a named list or a named atomic vector"
  )
})

test_that("target_strength verbosity works as intended", {
  data(krill, package = "acousticTS")
  data(sardine, package = "acousticTS")

  # Compute TS
  m1 <- capture_output_lines(
    target_strength(
      krill,
      frequency = 120e3, model = "DWBA", verbose = TRUE
    )
  )
  m2 <- capture_output_lines(
    target_strength(
      sardine,
      frequency = 120e3, model = "KRM", verbose = TRUE
    )
  )

  # Check first line: initialization confirmation
  expect_equal(
    m1[1],
    paste0(
      "DWBA model for FLS-object: Antarctic Euphausia superba ",
      "(McGehee et al., 1998) initialized."
    )
  )
  expect_equal(
    m2[1],
    paste0(
      "KRM model for SBF-object: Sardinops sagax caerulea ",
      "(Conti and Demer, 2003) initialized."
    )
  )

  # Check second line: model start
  expect_equal(
    m1[3],
    paste0(
      "Beginning TS modeling via DWBA model for FLS-object: Antarctic ",
      "Euphausia superba (McGehee et al., 1998) "
    )
  )
  expect_equal(
    m2[3],
    paste0(
      "Beginning TS modeling via KRM model for SBF-object: Sardinops sagax ",
      "caerulea (Conti and Demer, 2003) "
    )
  )

  # Check third line: model conclusion/success
  expect_equal(
    m1[4],
    paste0(
      "DWBA TS model predictions for FLS-object: Antarctic Euphausia ",
      "superba (McGehee et al., 1998) complete."
    )
  )
  expect_equal(
    m2[4],
    paste0(
      "KRM TS model predictions for SBF-object: Sardinops sagax ",
      "caerulea (Conti and Demer, 2003) complete."
    )
  )


  w <- capture_warnings(brake(krill, radius_curvature = 20e-3))
  expect_match(w, ".*Arc angle per segment", all = FALSE)
  expect_match(w, ".*One or more body segments", all = FALSE)
})
