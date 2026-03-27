library(acousticTS)

test_that("Mathematical utility functions work correctly", {

  # Test along_sum function
  rpos <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 2, ncol = 3)
  iterations <- 3
  result <- along_sum(rpos, iterations)
  expected <- rpos[, 1:(iterations - 1)] + rpos[, 2:iterations]
  expect_equal(result, expected)
  expect_error(
    along_sum(rpos, 1),
    "requires at least 2 columns"
  )

  # Test degrees conversion function
  expect_equal(degrees(pi), 180)
  expect_equal(degrees(pi / 2), 90)
  expect_equal(degrees(0), 0)
  expect_equal(degrees(2 * pi), 360)

  # Test radians conversion function
  expect_equal(radians(180), pi)
  expect_equal(radians(90), pi / 2)
  expect_equal(radians(0), 0)
  expect_equal(radians(360), 2 * pi)

  # Test vecnorm function (Euclidean norm)
  # Test single row
  values <- matrix(c(1, 2, 3), ncol = 3)
  expected_norm <- sqrt(1^2 + 2^2 + 3^2)
  expect_equal(vecnorm(values), expected_norm, tolerance = 1e-10)

  # Test multiple rows
  values_multi <- matrix(c(1, 2, 2, 3, 3, 4), nrow = 2, ncol = 3)
  expected_norms <- c(sqrt(1^2 + 2^2 + 3^2), sqrt(2^2 + 3^2 + 4^2))
  expect_equal(vecnorm(values_multi), expected_norms, tolerance = 1e-10)

  expect_equal(neumann(0:3), c(1, 2, 2, 2))
  expect_error(
    neumann(-1),
    "must be a non-negative integer"
  )
  expect_error(
    neumann(c(0, 1.5)),
    "All values in vector 'x' must be non-negative integers."
  )
})

test_that("Complex integration functions work", {

  # Test contour_integrate with a simple function
  # Define a simple complex function for testing
  simple_integral <- function(s, x, y) {
    complex(real = s, imaginary = s^2)
  }

  # Test the contour integration
  result <- contour_integrate(simple_integral, 1, 1)

  # The real part should be integral of s from 0 to 1 = 0.5
  # The imaginary part should be integral of s^2 from 0 to 1 = 1/3
  expect_equal(Re(result), 0.5, tolerance = 1e-6)
  expect_equal(Im(result), 1 / 3, tolerance = 1e-6)

  # Test that the result is complex
  expect_true(is.complex(result))

  phase_locked <- phase_integrate(
    x = 1,
    y = 1,
    n_iterations = 4,
    integral = simple_integral,
    phase_sd = 0
  )
  expect_equal(phase_locked, rep(result, 4))
})

test_that(
  "Gauss-Legendre quadrature validates n and supports the one-node rule", 
  {

    gl1 <- gauss_legendre(1)

    expect_equal(gl1$nodes, 0)
    expect_equal(gl1$weights, 2)

    expect_error(
      gauss_legendre(1.5),
      "n must be a positive integer"
    )
    expect_error(
      gauss_legendre(2, a = c(0, 1), b = 1),
      "a and b must be numeric scalars"
    )
    expect_error(
      gauss_legendre(2, a = 1, b = 1),
      "b must be greater than a"
    )
  }
)
