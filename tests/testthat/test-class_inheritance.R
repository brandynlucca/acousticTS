test_that("Scatterer class inheritance works correctly", {
  library(acousticTS)
  
  # Test that all scatterer classes inherit from Scatterer
  # Create instances of each scatterer type
  
  # Test CAL (calibration sphere) - inherits from ESS
  cal_obj <- cal_generate()
  expect_s4_class(cal_obj, "CAL")
  expect_s4_class(cal_obj, "ESS") # CAL inherits from ESS
  expect_s4_class(cal_obj, "Scatterer") # ESS inherits from Scatterer
  
  # Test ESS (elastic shelled sphere)  
  ess_obj <- ess_generate(radius_shell=1)
  expect_s4_class(ess_obj, "ESS")
  expect_s4_class(ess_obj, "Scatterer") # ESS inherits from Scatterer
  
  # Test GAS (gas-filled scatterer)
  gas_obj <- gas_generate(radius_body=1)
  expect_s4_class(gas_obj, "GAS")
  expect_s4_class(gas_obj, "Scatterer") # GAS inherits from Scatterer
  
  # Test FLS (fluid-like scatterer)
  fls_obj <- fls_generate(x=1, y=1, z=1, radius_body=1, g_body=1.00, 
                          h_body=1.00)
  expect_s4_class(fls_obj, "FLS")
  expect_s4_class(fls_obj, "Scatterer") # FLS inherits from Scatterer
  
  # Test SBF (swimbladdered fish) - inherits from GAS
  # Create a simple SBF object for testing
  x_body <- seq(0, 0.1, length.out = 11)
  w_body <- rep(0.01, 11)
  zU_body <- rep(0.005, 11)
  zL_body <- rep(-0.005, 11)
  x_bladder <- seq(0.02, 0.08, length.out = 7)
  w_bladder <- rep(0.005, 7)
  zU_bladder <- rep(0.003, 7)
  zL_bladder <- rep(-0.003, 7)
  
  sbf_obj <- sbf_generate(
    x_body = x_body, w_body = w_body, zU_body = zU_body, zL_body = zL_body,
    x_bladder = x_bladder, w_bladder = w_bladder, 
    zU_bladder = zU_bladder, zL_bladder = zL_bladder,
    sound_speed_body = 1570, sound_speed_bladder = 345,
    density_body = 1070, density_bladder = 1.24
  )
  
  expect_s4_class(sbf_obj, "SBF")
  expect_s4_class(sbf_obj, "GAS") # SBF inherits from GAS
  expect_s4_class(sbf_obj, "Scatterer") # GAS inherits from Scatterer
  
})

test_that("Scatterer objects have required slots", {
  library(acousticTS)
  
  # Test that Scatterer base class has required slots
  cal_obj <- cal_generate()
  
  # Check that slots can be accessed
  metadata <- cal_obj@metadata
  model_params <- cal_obj@model_parameters
  
  expect_type(metadata, "list")
  expect_type(model_params, "list")
  
  # Test with different scatterer types
  gas_obj <- gas_generate(radius_body=1)
  expect_true("metadata" %in% slotNames(gas_obj))
  expect_true("model_parameters" %in% slotNames(gas_obj))
  
  fls_obj <- fls_generate(x=1, y=1, z=1, radius_body=1, g_body=1, h_body=1)
  expect_true("metadata" %in% slotNames(fls_obj))
  expect_true("model_parameters" %in% slotNames(fls_obj))
})

test_that("Scatterer generation functions work", {
  library(acousticTS)
  
  # Test cal_generate
  cal_obj <- cal_generate()
  expect_s4_class(cal_obj, "CAL")
  expect_equal(cal_obj@metadata$ID, "Calibration sphere")  # Default ID
  
  # Test with custom ID
  cal_obj_custom <- cal_generate(ID = "TestSphere")
  expect_equal(cal_obj_custom@metadata$ID, "TestSphere")
  
  # Test ess_generate  
  ess_obj <- ess_generate(radius_shell=1)
  expect_s4_class(ess_obj, "ESS")
  expect_equal(ess_obj@metadata$ID, "UID")
  
  # Test gas_generate
  gas_obj <- gas_generate(radius_body=1)
  expect_s4_class(gas_obj, "GAS")
  expect_equal(gas_obj@metadata$ID, "UID")
  
  # Test fls_generate
  fls_obj <- fls_generate(x=1, y=1, z=1, radius_body=1, g_body=1, h_body=1)
  expect_s4_class(fls_obj, "FLS")
  expect_equal(fls_obj@metadata$ID, "UID")
})

test_that("Scatterer objects have proper structure", {
  library(acousticTS)
  
  # Test CAL structure
  cal_obj <- cal_generate()
  expect_true("body" %in% slotNames(cal_obj))
  expect_type(cal_obj@body, "list")
  
  # Test GAS structure  
  gas_obj <- gas_generate(radius_body=1)
  expect_true("body" %in% slotNames(gas_obj))
  expect_type(gas_obj@body, "list")
  
  # Test FLS structure
  fls_obj <- fls_generate(x=1, y=1, z=1, radius_body=1, g_body=1, h_body=1)
  expect_true("body" %in% slotNames(fls_obj))
  expect_type(fls_obj@body, "list")
  
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
  
  expect_true("body" %in% slotNames(sbf_obj))
  expect_true("bladder" %in% slotNames(sbf_obj))
  expect_type(sbf_obj@body, "list")
  expect_type(sbf_obj@bladder, "list")
})
