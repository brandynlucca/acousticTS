test_that("Generic show method works for all scatterer types", {
  library(acousticTS)
  
  # Suppression helper function
  suppressAllOutput <- function(expr) {
    suppressMessages(suppressWarnings(capture.output(expr, file = NULL)))
  }
  
  # Test show method for CAL objects
  cal_obj <- cal_generate()
  # Should not throw an error
  expect_error(suppressAllOutput(show(cal_obj)), NA)
  
  # Test show method for ESS objects
  ess_obj <- ess_generate(radius_shell=1)
  expect_error(suppressAllOutput(show(ess_obj)), NA)
  
  # Test show method for GAS objects
  gas_obj <- gas_generate(radius_body=1)
  expect_error(suppressAllOutput(show(gas_obj)), NA)
  
  # Test show method for FLS objects
  fls_obj <- fls_generate(x=1, y=1, z=1, radius_body=1, g_body=1, h_body=1)
  expect_error(suppressAllOutput(show(fls_obj)), NA)
  
  # Test show method for SBF objects
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
  expect_error(suppressAllOutput(show(sbf_obj)), NA)
})

test_that("Plot methods work for scatterer objects", {
  library(acousticTS)
  
  # Helper suppression function
  suppressPlot <- function(expr) {
    pdf(NULL)             # Open null graphics device
    on.exit(dev.off())    # Ensure it's closed after evaluation
    invisible(force(expr))
  }
  
  # Test plot method for CAL objects
  cal_obj <- cal_generate()
  expect_error(suppressPlot(plot(cal_obj)), NA)
  
  # Test plot method for FLS objects
  fls_obj <- fls_generate(x=c(1, 2, 3), 
                          y=rep(0, 3), 
                          z=rep(0, 3), 
                          radius_body=rep(3, 3), 
                          g_body=1, h_body=1)
  expect_error(suppressPlot(plot(fls_obj)), NA)
  
  # Test plot method for GAS objects
  gas_obj <- gas_generate(radius_body=1)
  expect_error(suppressPlot(plot(gas_obj)), NA)
  
  # Test plot method for SBF objects
  x_body <- seq(0, 0.04, length.out = 5)
  w_body <- rep(0.004, 5)
  zU_body <- rep(0.002, 5)
  zL_body <- rep(-0.002, 5)
  x_bladder <- seq(0.008, 0.032, length.out = 3)
  w_bladder <- rep(0.0015, 3)
  zU_bladder <- rep(0.001, 3)
  zL_bladder <- rep(-0.001, 3)
  
  sbf_obj <- sbf_generate(
    x_body = x_body, w_body = w_body, zU_body = zU_body, zL_body = zL_body,
    x_bladder = x_bladder, w_bladder = w_bladder,
    zU_bladder = zU_bladder, zL_bladder = zL_bladder,
    sound_speed_body = 1570, sound_speed_bladder = 345,
    density_body = 1070, density_bladder = 1.24
  )
  expect_error(suppressPlot(plot(sbf_obj)), NA)
})

test_that("Extract method works for scatterer objects", {
  library(acousticTS)
  
  # Test extract method for metadata
  cal_obj <- cal_generate()
  metadata <- extract(cal_obj, "metadata")
  expect_type(metadata, "list")
  expect_true("ID" %in% names(metadata))
  
  # Test extract method for shape_parameters
  shape_params <- extract(cal_obj, "shape_parameters")
  expect_type(shape_params, "list")
  
  # Test extract method for body
  body <- extract(cal_obj, "body")
  expect_type(body, "list")
  
  # Test with different scatterer types
  fls_obj <- fls_generate(x=1, y=1, z=1, radius_body=1, g_body=1, h_body=1)
  fls_metadata <- extract(fls_obj, "metadata")
  expect_type(fls_metadata, "list")
  
  gas_obj <- gas_generate(radius_body=1)
  gas_body <- extract(gas_obj, "body")
  expect_type(gas_body, "list")
})
