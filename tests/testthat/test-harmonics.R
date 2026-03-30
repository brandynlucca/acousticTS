library(acousticTS)

test_that("Legendre wrappers handle standard values and validation branches", {
  expect_equal(drop(Pn(0, 0.5)), 1)
  expect_equal(drop(Pn(1, 0.5)), 0.5)
  expect_equal(drop(Pn(2, 0.5)), -0.125, tolerance = 1e-10)
  expect_equal(drop(Pn(3, -1)), -1, tolerance = 1e-10)

  p_mat <- Pn(c(0, 1, 2.5), c(-1, 0, 1.5))
  expect_equal(dim(p_mat), c(3, 3))

  expect_equal(drop(Pndk(2, 0.5, 1)), 1.5, tolerance = 1e-8)
  expect_equal(drop(Pndk(2, 0.5, 2)), 3, tolerance = 1e-8)
  expect_equal(drop(Pndk(3, 1, 1)), 6, tolerance = 1e-8)
  expect_equal(drop(Pndk(2, 0.5, 3)), 0, tolerance = 1e-8)

  p_frac_fd <- (Pn(0.5, 0.25 + 1e-6) - Pn(0.5, 0.25 - 1e-6)) / (2e-6)
  expect_equal(drop(Pndk(0.5, 0.25, 1)), drop(p_frac_fd), tolerance = 1e-4)

  expect_error(Pn(1 + 1i, 0.5), "must be real numbers")
  expect_error(Pn(Inf, 0.5), "'n' must be finite.")
  expect_error(Pn(1, Inf), "must be a real, finite numeric")
  expect_error(Pndk(Inf, 0.5, 1), "'n' must be finite.")
  expect_error(Pndk(1, Inf, 1), "'x' must be a real, finite numeric.")
  expect_error(Pndk(2, 0.5, -1), "must be a non-negative integer")
})

test_that("Second-kind Legendre wrappers return expected values and warnings", {
  q0_expected <- 0.5 * log((1 + 0.5) / (1 - 0.5))
  q1_expected <- 0.5 * 0.5 * log((1 + 0.5) / (1 - 0.5)) - 1

  expect_equal(Re(drop(Qn(0, 0.5))), q0_expected, tolerance = 1e-8)
  expect_equal(Re(drop(Qn(1, 0.5))), q1_expected, tolerance = 1e-8)

  q_mat <- Qn(c(0, 1.5), c(0.25, 0.75))
  expect_equal(dim(q_mat), c(2, 2))

  q_complex <- drop(Qn(1, 2))
  expect_true(is.complex(q_complex))
  expect_gt(abs(Im(q_complex)), 0)

  fd <- (Qn(1, 0.5 + 1e-6) - Qn(1, 0.5 - 1e-6)) / (2e-6)
  expect_equal(drop(Qndk(1, 0.5, 1)), drop(fd), tolerance = 1e-4)

  fd2 <- (Qn(1, 0.5 + 1e-6) - 2 * Qn(1, 0.5) + Qn(1, 0.5 - 1e-6)) / (1e-6^2)
  expect_equal(drop(Qndk(1, 0.5, 2)), drop(fd2), tolerance = 1e-3)

  expect_warning(
    Qndk(1, 1, 1),
    "Derivatives near x = \\+/\\-1 may be inaccurate"
  )
  expect_error(Qn(1 + 1i, 0.5), "must be real numbers")
  expect_error(Qn(1, Inf), "'x' must be a real, finite numeric.")
  expect_error(Qn(Inf, 0.5), "'n' must be finite.")
  expect_error(acousticTS:::.validate_Qndk_inputs(Inf, 0.5, 1), "'n' must be finite.")
  expect_error(acousticTS:::.validate_Qndk_inputs(1, Inf, 1), "'x' must be a real, finite numeric.")
  expect_error(Qndk(1, 0.5 + 1i, 1), "must be real numbers")
  expect_error(Qndk(1, 0.5, 1.5), "must be a non-negative integer")
})
