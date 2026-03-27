library(acousticTS)

test_that(
  "Compare calibration sphere model output at 38, 70, 120, and 200 kHz",
  {

    # Create class
    cal_sphere <- acousticTS::cal_generate()
    # Class check
    expect_true(methods::is(cal_sphere, "CAL"))
    # Parameterize model
    frequency <- c(38e3, 70e3, 120e3, 200e3)
    # Calculate TS; update original CAL object
    cal_sphere <- target_strength(
      object = cal_sphere,
      frequency = frequency,
      model = "calibration"
    )
    # Calculate TS; store in a new CAL object
    cal_sphere_copy <- target_strength(
      object = cal_sphere,
      frequency = frequency,
      model = "calibration"
    )
    # Class check
    expect_true(methods::is(cal_sphere_copy, "CAL"))
    # Extract model results
    model_results <- acousticTS::extract(cal_sphere,
                                         "model")$calibration$TS
    model_results2 <- acousticTS::extract(cal_sphere_copy,
                                          "model")$calibration$TS
    ts_out <- model_results
    ts_out2 <- model_results2
    expect_equal(ts_out, ts_out2)
    # Check output
    # Should be -42.32, -41.07, -39.50, and -39.44 dB at 38, 70, 120, and
    # 200 kHz
    expect_equal(
      ts_out,
      c(-42.3296920811525,
        -41.07033000976835,
        -39.50264169247166,
        -39.43821198162192)
    )
  }
)
