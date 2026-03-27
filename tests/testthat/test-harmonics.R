library(acousticTS)

test_that("Legendre wrappers handle standard values and validation branches", {
  expect_equal(drop(Pn(0, 0.5)), 1)
  expect_equal(drop(Pn(1, 0.5)), 0.5)
  expect_equal(drop(Pn(2, 0.5)), -0.125, tolerance = 1e-10)

  expect_equal(drop(Pndk(2, 0.5, 1)), 1.5, tolerance = 1e-8)
  expect_equal(drop(Pndk(2, 0.5, 2)), 3, tolerance = 1e-8)

  expect_error(Pn(1 + 1i, 0.5), "must be real numbers")
  expect_error(Pn(1, Inf), "must be a real, finite numeric")
  expect_error(Pndk(2, 0.5, -1), "must be a non-negative integer")
})

test_that("Second-kind Legendre wrappers return expected values and warnings", {
  q0_expected <- 0.5 * log((1 + 0.5) / (1 - 0.5))
  q1_expected <- 0.5 * 0.5 * log((1 + 0.5) / (1 - 0.5)) - 1

  expect_equal(Re(drop(Qn(0, 0.5))), q0_expected, tolerance = 1e-8)
  expect_equal(Re(drop(Qn(1, 0.5))), q1_expected, tolerance = 1e-8)

  fd <- (Qn(1, 0.5 + 1e-6) - Qn(1, 0.5 - 1e-6)) / (2e-6)
  expect_equal(drop(Qndk(1, 0.5, 1)), drop(fd), tolerance = 1e-4)

  expect_warning(
    Qndk(1, 1, 1),
    "Derivatives near x = \\+/\\-1 may be inaccurate"
  )
  expect_error(Qn(Inf, 0.5), "'n' must be finite.")
  expect_error(Qndk(1, 0.5 + 1i, 1), "must be real numbers")
})
