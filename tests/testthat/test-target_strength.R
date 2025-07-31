test_that("target_strength function works correctly", {
  library(acousticTS)
  
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

test_that("target_strength works with different scatterer types", {
  library(acousticTS)
  
  # Test with FLS object and DWBA model
  data(krill, package="acousticTS")
  fls_obj <- fls_generate(
    x=krill@body$rpos[1, ],
    y=krill@body$rpos[2, ],
    z=krill@body$rpos[3, ],
    radius_body=krill@body$radius,
    g_body=krill@body$g,
    h_body=krill@body$h,
    radius_curvature_ratio=krill@body$radius_curvature_ratio,
    theta_body=krill@body$theta
  )
  frequency <- seq(1e3, 100e3, 1e3)
  
  fls_with_ts <- target_strength(
    object = fls_obj,
    frequency = frequency,
    model = "DWBA"
  )
  
  expect_s4_class(fls_with_ts, "FLS")
  expect_true("DWBA" %in% names(fls_with_ts@model))
  
  # Test with GAS object and MSS Anderson (1950) model
  gas_obj <- gas_generate(radius_body=1e-3)
  
  gas_with_ts <- target_strength(
    object = gas_obj,
    frequency = frequency,
    model = "mss_anderson"
  )
  
  expect_s4_class(gas_with_ts, "GAS")
  expect_true("MSS_anderson" %in% names(gas_with_ts@model))
})

test_that("target_strength can update existing models", {
  library(acousticTS)
  
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
  library(acousticTS)
  
  # Test with FLS that can use multiple models
  data(krill, package="acousticTS")
  fls_obj <- fls_generate(
    x=krill@body$rpos[1, ],
    y=krill@body$rpos[2, ],
    z=krill@body$rpos[3, ],
    radius_body=krill@body$radius,
    g_body=krill@body$g,
    h_body=krill@body$h,
    radius_curvature_ratio=krill@body$radius_curvature_ratio,
    theta_body=krill@body$theta
  )
  frequency <- c(120e3)
  
  # Calculate with DWBA
  fls_dwba <- target_strength(fls_obj, frequency, "DWBA")
  
  # Calculate with DCM (should add to existing models)
  fls_both <- target_strength(fls_dwba, frequency, "DCM")
  
  # Should have both models
  expect_true("DWBA" %in% names(fls_both@model))
  expect_true("DCM" %in% names(fls_both@model))
})

test_that("target_strength preserves object metadata", {
  library(acousticTS)
  
  # Test that metadata is preserved after target_strength calculation
  cal_obj <- cal_generate(ID = "TestSphere")
  original_metadata <- cal_obj@metadata
  
  cal_with_ts <- target_strength(cal_obj, 38e3, "calibration")
  
  # Metadata should be unchanged
  expect_equal(cal_with_ts@metadata, original_metadata)
  expect_equal(cal_with_ts@metadata$ID, "TestSphere")
})
