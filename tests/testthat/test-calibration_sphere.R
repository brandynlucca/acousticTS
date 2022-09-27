test_that("Compare calibration sphere model output at 38, 70, 120, and 200 kHz", {
  cal_sphere <- cal_generate()
  # Class check
  expect_true(class(cal_sphere) == "CAL")
  # Parameterize model
  frequency <- c(38e3, 70e3, 120e3, 200e3)
  # Calculate TS; update original CAL object
  cal_sphere <- target_strength(object = cal_sphere,
                                frequency = frequency,
                                model = "calibration")
  # Calculate TS; store in a new CAL object
  cal_sphere_copy <- target_strength(object = cal_sphere,
                                     frequency = frequency,
                                     model = "calibration")
  # Class check
  expect_true(class(cal_sphere_copy) == "CAL")
  # Extract model results
  model_results <- extract(cal_sphere, "model")
  ts_out <- round(model_results$TS, 2)
  # Check output
  # Should be -42.42, -41.44, -39.52, and -39.05 dB at 38, 70, 120, and 200 kHz
  expect_equal(ts_out, c(-42.42, -41.44, -39.52, -39.05))
})
