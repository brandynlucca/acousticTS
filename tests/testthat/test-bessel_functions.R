test_that("Bessel functions work correctly", {
  library(acousticTS)
  
  # Test jc (cylindrical Bessel function of the first kind)
  expect_equal(jc(0, 1), besselJ(1, 0))
  expect_equal(jc(1, 2), besselJ(2, 1))
  expect_equal(jc(0.5, 1.5), besselJ(1.5, 0.5))
  
  # Test yc (cylindrical Bessel function of the second kind)
  expect_equal(yc(0, 1), besselY(1, 0))
  expect_equal(yc(1, 2), besselY(2, 1))
  expect_equal(yc(0.5, 1.5), besselY(1.5, 0.5))
  
  # Test jcd (first derivative of jc)
  # Test edge case where n = 0 and l = 1
  expect_equal(jcd(1, 0), 0.5)
  # Test edge case where n = 0 and l != 1
  expect_equal(jcd(2, 0), 0.0)
  # Test normal case
  l <- 1
  n <- 2
  expected <- jc(l - 1, n) - (l / n) * jc(l, n)
  expect_equal(jcd(l, n), expected)
  
  # Test jcd with matrix input
  n_matrix <- matrix(c(1, 2, 3, 4), ncol = 2)
  result <- jcd(1, n_matrix)
  expect_true(is.matrix(result))
  expect_equal(ncol(result), 2)
  
  # Test ycd (first derivative of yc)
  # Test edge case where n = 0 and l = 1
  expect_equal(ycd(1, 0), 0.5)
  # Test edge case where n = 0 and l != 1
  expect_equal(ycd(2, 0), 0.0)
  # Test normal case
  l <- 1
  n <- 2
  expected <- yc(l - 1, n) - (l / n) * yc(l, n)
  expect_equal(ycd(l, n), expected)
  
  # Test jcdd (second derivative of jc)
  l <- 1
  n <- 2
  expected <- 0.25 * (jc(l - 2, n) - 2 * jc(l, n) + jc(l + 2, n))
  expect_equal(jcdd(l, n), expected)
})

test_that("Spherical Bessel functions work correctly", {
  library(acousticTS)
  
  # Test js (spherical Bessel function of the first kind)
  # js should be related to jc by: js(l, n) = sqrt(pi/(2*n)) * jc(l + 0.5, n)
  l <- 1
  n <- 2
  expected <- sqrt(pi / (2 * n)) * jc(l + 0.5, n)
  expect_equal(js(l, n), expected, tolerance = 1e-10)
  
  # Test ys (spherical Bessel function of the second kind)
  # ys should be related to yc by: ys(l, n) = sqrt(pi/(2*n)) * yc(l + 0.5, n)
  expected_ys <- sqrt(pi / (2 * n)) * yc(l + 0.5, n)
  expect_equal(ys(l, n), expected_ys, tolerance = 1e-10)
  
  # Test hs (spherical Hankel function)
  # hs should equal js + 1i * ys
  expected_hs <- js(l, n) + 1i * ys(l, n)
  expect_equal(hs(l, n), expected_hs, tolerance = 1e-10)
})

test_that("Hankel functions work correctly", {
  library(acousticTS)
  
  # Test hc (cylindrical Hankel function)
  # hc should equal jc + 1i * yc
  l <- 1
  n <- 2
  expected <- jc(l, n) + 1i * yc(l, n)
  expect_equal(hc(l, n), expected, tolerance = 1e-10)
  
  # Test hcd (first derivative of hc)
  # hcd should equal jcd + 1i * ycd
  expected_hcd <- jcd(l, n) + 1i * ycd(l, n)
  expect_equal(hcd(l, n), expected_hcd, tolerance = 1e-10)
})
