test_that("calibration_initialize works correctly", {
  library(acousticTS)
  
  # Test calibration_initialize
  cal_obj <- cal_generate()
  frequency <- c(38e3, 70e3, 120e3)
  
  # Initialize calibration model
  cal_initialized <- calibration_initialize(cal_obj, frequency = frequency)
  
  # Check that object is still CAL class
  expect_s4_class(cal_initialized, "CAL")
  
  # Check that model parameters slot was updated
  expect_type(cal_initialized@model_parameters, "list")
  expect_true("calibration" %in% names(cal_initialized@model_parameters))
  
  # Check that model slot was updated
  expect_type(cal_initialized@model, "list")
  expect_true("calibration" %in% names(cal_initialized@model))
  
  # Check that frequency was properly set
  expect_equal(cal_initialized@model$calibration$frequency, frequency)
})

test_that("mss_anderson_initialize works correctly", {
  library(acousticTS)
  
  # Test mss_anderson_initialize for GAS objects
  gas_obj <- gas_generate(radius_body=1)
  frequency <- c(38e3, 120e3)
  
  # Initialize MSS Anderson model
  gas_initialized <- mss_anderson_initialize(gas_obj, frequency = frequency)
  
  # Check that object is still GAS class
  expect_s4_class(gas_initialized, "GAS")
  
  # Check that model parameters slot was updated
  expect_type(gas_initialized@model_parameters, "list")
  expect_true("MSS_anderson" %in% names(gas_initialized@model_parameters))
  
  # Check that model slot was updated
  expect_type(gas_initialized@model, "list")
  expect_true("MSS_anderson" %in% names(gas_initialized@model))
  
  # Check that frequency was properly set
  expect_equal(gas_initialized@model$MSS_anderson$frequency, frequency)
})

test_that("mss_goodman_stern_initialize works correctly", {
  library(acousticTS)
  
  # Test mss_goodman_stern_initialize for ESS objects
  ess_obj <- ess_generate(radius_shell=1)
  frequency <- c(70e3, 200e3)
  
  # Initialize MSS Goodman-Stern model
  ess_initialized <- mss_goodman_stern_initialize(ess_obj, 
                                                  frequency = frequency)
  
  # Check that object is still ESS class
  expect_s4_class(ess_initialized, "ESS")
  
  # Check that model parameters slot was updated
  expect_type(ess_initialized@model_parameters, "list")
  expect_true("MSS_goodman_stern" %in% names(ess_initialized@model_parameters))
  
  # Check that model slot was updated
  expect_type(ess_initialized@model, "list")
  expect_true("MSS_goodman_stern" %in% names(ess_initialized@model))
})

test_that("dwba_initialize works correctly", {
  library(acousticTS)
  
  # Test dwba_initialize for FLS objects
  fls_obj <- fls_generate(x=1, y=1, z=1, radius_body=1, g_body=1, h_body=1)
  frequency <- c(38e3, 120e3)
  
  # Initialize DWBA model
  fls_initialized <- dwba_initialize(fls_obj, frequency = frequency)
  
  # Check that object is still FLS class
  expect_s4_class(fls_initialized, "FLS")
  
  # Check that model parameters slot was updated
  expect_type(fls_initialized@model_parameters, "list")
  expect_true("DWBA" %in% names(fls_initialized@model_parameters))
  
  # Check that model slot was updated
  expect_type(fls_initialized@model, "list")
  expect_true("DWBA" %in% names(fls_initialized@model))
  
  # Check that frequency was properly set
  expect_equal(fls_initialized@model$DWBA$frequency, frequency)
})

test_that("dwba_curved_initialize works correctly", {
  library(acousticTS)
  
  # Test dwba_curved_initialize for FLS objects
  fls_obj <- fls_generate(
    x=seq(0, 1, length.out=11), 
    y=rep(0, 11), 
    z=rep(0, 11), 
    radius_body=rep(0.5, 11), 
    g_body=1, h_body=1,
    radius_curvature_ratio=3
  )
  frequency <- c(120e3)
  
  # Initialize DWBA curved model
  fls_initialized <- dwba_curved_initialize(fls_obj, frequency = frequency)
  
  # Check that object is still FLS class
  expect_s4_class(fls_initialized, "FLS")
  
  # Check that model parameters slot was updated
  expect_type(fls_initialized@model_parameters, "list")
  expect_true("DWBA_curved" %in% names(fls_initialized@model_parameters))
  
  # Check that model slot was updated
  expect_type(fls_initialized@model, "list")
  expect_true("DWBA_curved" %in% names(fls_initialized@model))
})

test_that("dcm_initialize works correctly", {
  library(acousticTS)
  
  # Test dcm_initialize for FLS objects
  fls_obj <- fls_generate(x=1, y=1, z=1, radius_body=1, g_body=1, h_body=1)
  frequency <- c(38e3, 120e3)
  
  # Initialize DCM model
  fls_initialized <- dcm_initialize(fls_obj, frequency = frequency)
  
  # Check that object is still FLS class
  expect_s4_class(fls_initialized, "FLS")
  
  # Check that model parameters slot was updated
  expect_type(fls_initialized@model_parameters, "list")
  expect_true("DCM" %in% names(fls_initialized@model_parameters))
  
  # Check that model slot was updated
  expect_type(fls_initialized@model, "list")
  expect_true("DCM" %in% names(fls_initialized@model))
})

test_that("sdwba_initialize works correctly", {
  library(acousticTS)
  
  # Test sdwba_initialize for FLS objects
  fls_obj <- fls_generate(
    x=seq(0, 1, length.out=11), 
    y=rep(0, 11), 
    z=rep(0, 11), 
    radius_body=rep(0.5, 11), 
    g_body=1, h_body=1,
    radius_curvature_ratio=3
  )
  frequency <- c(120e3, 200e3)
  
  # Initialize SDWBA model
  fls_initialized <- sdwba_initialize(fls_obj, frequency = frequency)
  
  # Check that object is still FLS class
  expect_s4_class(fls_initialized, "FLS")
  
  # Check that model parameters slot was updated
  expect_type(fls_initialized@model_parameters, "list")
  expect_true("SDWBA" %in% names(fls_initialized@model_parameters))
  
  # Check that model slot was updated
  expect_type(fls_initialized@model, "list")
  expect_true("SDWBA" %in% names(fls_initialized@model))
})

test_that("sdwba_curved_initialize works correctly", {
  library(acousticTS)
  
  # Test sdwba_curved_initialize for FLS objects
  fls_obj <- fls_generate(
    x=seq(0, 1, length.out=11), 
    y=rep(0, 11), 
    z=rep(0, 11), 
    radius_body=rep(0.5, 11), 
    g_body=1, h_body=1,
    radius_curvature_ratio=3
  )
  frequency <- c(120e3)
  
  # Initialize SDWBA curved model
  fls_initialized <- sdwba_curved_initialize(fls_obj, frequency = frequency)
  
  # Check that object is still FLS class
  expect_s4_class(fls_initialized, "FLS")
  
  # Check that model parameters slot was updated
  expect_type(fls_initialized@model_parameters, "list")
  expect_true("SDWBA_curved" %in% names(fls_initialized@model_parameters))
  
  # Check that model slot was updated
  expect_type(fls_initialized@model, "list")
  expect_true("SDWBA_curved" %in% names(fls_initialized@model))
})

test_that("krm_initialize works correctly", {
  library(acousticTS)
  
  # Test krm_initialize for SBF objects
  # SBF should have both body and bladder
  x_body <- seq(0, 0.05, length.out = 6)
  w_body <- rep(0.005, 6)
  zU_body <- rep(0.0025, 6)
  zL_body <- rep(-0.0025, 6)
  x_bladder <- seq(0.01, 0.04, length.out = 4)
  w_bladder <- rep(0.002, 4)
  zU_bladder <- rep(0.001, 4)
  zL_bladder <- rep(-0.001, 4)
  
  sbf_obj <- sbf_generate(
    x_body = x_body, w_body = w_body, zU_body = zU_body, zL_body = zL_body,
    x_bladder = x_bladder, w_bladder = w_bladder,
    zU_bladder = zU_bladder, zL_bladder = zL_bladder,
    sound_speed_body = 1570, sound_speed_bladder = 345,
    density_body = 1070, density_bladder = 1.24
  )
  frequency <- c(38e3, 70e3)
  
  # Initialize KRM model
  sbf_initialized <- krm_initialize(sbf_obj, frequency = frequency)
  
  # Check that object is still SBF class
  expect_s4_class(sbf_initialized, "SBF")
  
  # Check that model parameters slot was updated
  expect_type(sbf_initialized@model_parameters, "list")
  expect_true("KRM" %in% names(sbf_initialized@model_parameters))
  
  # Check that model slot was updated
  expect_type(sbf_initialized@model, "list")
  expect_true("KRM" %in% names(sbf_initialized@model))
})

test_that("high_pass_stanton_initialize works correctly", {
  library(acousticTS)
  
  # Test high_pass_stanton_initialize for ESS objects
  ess_obj <- ess_generate(radius_shell=1)
  frequency <- c(120e3, 200e3)
  
  # Initialize High Pass Stanton model
  ess_initialized <- high_pass_stanton_initialize(ess_obj, 
                                                  frequency = frequency)
  
  # Check that object is still ESS class
  expect_s4_class(ess_initialized, "ESS")
  
  # Check that model parameters slot was updated
  expect_type(ess_initialized@model_parameters, "list")
  expect_true("high_pass_stanton" %in% names(ess_initialized@model_parameters))
  
  # Check that model slot was updated
  expect_type(ess_initialized@model, "list")
  expect_true("high_pass_stanton" %in% names(ess_initialized@model))
})

test_that("Initialization preserves object integrity", {
  library(acousticTS)
  
  # Test that initialization doesn't break the object structure
  cal_obj <- cal_generate()
  original_metadata <- cal_obj@metadata
  
  cal_initialized <- calibration_initialize(cal_obj, frequency = c(38e3))
  
  # Metadata should be preserved
  expect_equal(cal_initialized@metadata, original_metadata)
  
  # Object should still have all required slots
  expect_true("metadata" %in% slotNames(cal_initialized))
  expect_true("model_parameters" %in% slotNames(cal_initialized))
  expect_true("model" %in% slotNames(cal_initialized))
})

test_that("Initialization functions handle optional parameters", {
  library(acousticTS)
  
  # Test with optional medium parameters
  gas_obj <- gas_generate(radius_body=1)
  frequency <- c(38e3, 120e3)
  
  # Test with custom sound speed and density
  gas_initialized <- mss_anderson_initialize(
    gas_obj, 
    frequency = frequency,
    sound_speed_sw = 1480,  # Different from default
    density_sw = 1025       # Different from default
  )
  
  # Check that custom medium parameters were used
  medium_params <- gas_initialized@model_parameters$MSS_anderson$medium
  expect_equal(medium_params$sound_speed, 1480)
  expect_equal(medium_params$density, 1025)
})
