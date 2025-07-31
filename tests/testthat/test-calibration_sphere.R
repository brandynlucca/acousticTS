test_that("Compare calibration sphere model output at 38, 70, 120, and 200 kHz", {
  library( acousticTS )
  cal_sphere <- acousticTS::cal_generate( )
  # Class check
  expect_true(is(cal_sphere, "CAL"))
  # Parameterize model
  frequency <- c(38e3, 70e3, 120e3, 200e3)
  # Calculate TS; update original CAL object
  cal_sphere <- acousticTS::target_strength(object = cal_sphere,
                                            frequency = frequency,
                                            model = "calibration")
  # Calculate TS; store in a new CAL object
  cal_sphere_copy <- acousticTS::target_strength(object = cal_sphere,
                                                 frequency = frequency,
                                                 model = "calibration")
  # Class check
  expect_true(is(cal_sphere_copy, "CAL"))
  # Extract model results
  model_results <- acousticTS::extract(cal_sphere , "model")$calibration$TS
  model_results2 <- acousticTS::extract(cal_sphere_copy , "model")$calibration$TS
  ts_out <- model_results
  ts_out2 <- model_results2
  expect_equal( ts_out , ts_out2 )
  # Check output
  # Should be -42.42, -41.44, -39.52, and -39.05 dB at 38, 70, 120, and 200 kHz
  expect_equal( ts_out , 
                c( -42.4161095 , -41.4364324 , -39.5248183 , -39.0544814 ) ,
                tolerance = 1e-8 )
})
