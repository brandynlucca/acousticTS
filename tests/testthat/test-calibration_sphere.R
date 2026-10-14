library(acousticTS)

# Use scalar public special functions for the independent reference. ===========
.calibration_scalar_reference <- function(object, frequency, adaptive = TRUE) {
  body <- object@body
  x <- 2 * pi * frequency * body$radius / 1500
  xl <- x * 1500 / body$sound_speed_longitudinal
  xt <- x * 1500 / body$sound_speed_transversal
  term <- function(m) {
    sw <- -x * jsd(m, x) / js(m, x)
    longitudinal <- -xl * jsd(m, xl) / js(m, xl)
    transversal <- -xt * jsd(m, xt) / js(m, xt)
    beta <- -x * ysd(m, x) / ys(m, x)
    mm <- m * (m + 1)
    td <- mm - 1 - xt^2 / 2 + transversal
    numerator <- longitudinal / (longitudinal + 1) - mm / td
    denominator <- (mm - xt^2 / 2 + 2 * longitudinal) /
      (longitudinal + 1) - mm * (transversal + 1) / td
    phi <- xt^2 / 2 * numerator / denominator / (body$density / 1026)
    eta <- -js(m, x) / ys(m, x) * (phi + sw) / (phi + beta)
    cosine <- 1 / sqrt(1 + eta^2)
    sine <- eta * cosine
    (2 * m + 1) * as.numeric(Pn(m, cos(body$theta))) *
      sine * (1i * cosine - sine)
  }
  ml <- round(x) + 10L
  terms <- vapply(0:ml, term, complex(1))
  if (adaptive) {
    while (Mod(tail(terms, 1)) > 1e-10 && ml < 500L) {
      ml <- ml + 1L
      terms <- c(terms, term(ml))
    }
  }
  Mod(sum(terms)) * body$radius / x
}

test_that("calibration sequences preserve scalar modal solutions", {
  frequencies <- c(1000, 38000, 120000, 1000000)
  for (material in c("WC", "Cu", "Al")) {
    for (theta in c(0, 0.7, pi)) {
      object <- cal_generate(material = material, theta_sphere = theta)
      expected <- vapply(frequencies, function(f)
        .calibration_scalar_reference(object, f), numeric(1))
      actual <- target_strength(object, frequencies, "calibration",
        sound_speed_sw = 1500, density_sw = 1026)
      expect_equal(actual@model$calibration$f_bs, expected, tolerance = 1e-11)
      expect_equal(actual@model$calibration$TS, 20 * log10(expected),
        tolerance = 1e-7)
    }
  }
})

test_that(
  "Compare calibration sphere model output at 38, 70, 120, and 200 kHz",
  {
    # Create class =============================================================
    cal_sphere <- acousticTS::cal_generate()
    # Class check ==============================================================
    expect_true(methods::is(cal_sphere, "CAL"))
    # Parameterize model =======================================================
    frequency <- c(38e3, 70e3, 120e3, 200e3)
    # Calculate TS; update original CAL object =================================
    cal_sphere <- target_strength(
      object = cal_sphere,
      frequency = frequency,
      model = "calibration"
    )
    # Calculate TS; store in a new CAL object ==================================
    cal_sphere_copy <- target_strength(
      object = cal_sphere,
      frequency = frequency,
      model = "calibration"
    )
    # Class check ==============================================================
    expect_true(methods::is(cal_sphere_copy, "CAL"))
    # Extract model results ====================================================
    model_results <- acousticTS::extract(
      cal_sphere,
      "model"
    )$calibration$TS
    model_results2 <- acousticTS::extract(
      cal_sphere_copy,
      "model"
    )$calibration$TS
    ts_out <- model_results
    ts_out2 <- model_results2
    expect_equal(ts_out, ts_out2)
    # Check output =============================================================
    # Should be -42.32, -41.07, -39.50, and -39.44 dB at 38, 70, 120, and ======
    # 200 kHz ==================================================================
    expect_equal(
      ts_out,
      c(
        -42.3296920811525,
        -41.07033000976835,
        -39.50264169247166,
        -39.43821198162192
      )
    )
  }
)
