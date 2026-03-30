library(acousticTS)

test_that("Bessel functions work correctly", {
  # Test jc (cylindrical Bessel function of the first kind)
  expect_equal(Re(jc(0, 1)), besselJ(1, 0))
  expect_equal(Re(jc(1, 2)), besselJ(2, 1))
  expect_equal(Re(jc(0.5, 1.5)), besselJ(1.5, 0.5))

  # Test yc (cylindrical Bessel function of the second kind)
  expect_equal(Re(yc(0, 1)), besselY(1, 0))
  expect_equal(Re(yc(1, 2)), besselY(2, 1))
  expect_equal(Re(yc(0.5, 1.5)), besselY(1.5, 0.5))

  # Test matrix input for `yc`
  expect_equal(
    Re(t(yc(1, matrix(c(1, 2, 3))))),
    t(besselY(c(1, 2, 3), 1))
  )

  # Test jcd (first derivative of jc)
  # Test edge case where n = 0 and l = 1
  expect_equal(Re(jcd(1, 0)), 0.5)
  # Test edge case where n = 0 and l != 1
  expect_equal(Re(jcd(2, 0)), 0.0)
  # Test normal case
  l <- 1
  n <- 2
  expected <- jc(l - 1, n) - (l / n) * jc(l, n)
  expect_equal(jcd(l, n), expected)

  # Test jcd with matrix input
  n_l <- seq(1, 4, 1)
  n_z <- seq(1, 6, 1)
  result <- jcd(n_l, n_z)
  expect_true(is.matrix(result))
  expect_equal(dim(result), c(4, 6))

  # Test ycd (first derivative of yc)
  # Test edge case where n = 0 and l = 1
  expect_true(is.na(Re(ycd(1, 0))))
  # Test edge case where n = 0 and l != 1
  expect_true(is.na(Re(ycd(2, 0))))
  # Test normal case
  l <- 1
  n <- 2
  expected <- yc(l - 1, n) - (l / n) * yc(l, n)
  expect_equal(ycd(l, n), expected)

  # Test matrix input for `ycd`
  expect_equal(
    Re(as.vector(ycd(1, matrix(c(1, 2, 3))))),
    besselY(c(1, 2, 3), 0) - (1 / c(1, 2, 3)) * besselY(c(1, 2, 3), 1)
  )

  # Test jcdd (second derivative of jc)
  l <- 1
  n <- 2
  expected <- 0.25 * (jc(l - 2, n) - 2 * jc(l, n) + jc(l + 2, n))
  expect_equal(jcdd(l, n), expected)
})

test_that("Spherical Bessel functions work correctly", {
  # Test js (spherical Bessel function of the first kind)
  # js should be related to jc by: js(l, n) = sqrt(pi/(2*n)) * jc(l + 0.5, n)
  l <- 1
  n <- 2
  expected <- sqrt(pi / (2 * n)) * jc(l + 0.5, n)
  expect_equal(js(l, n), Re(expected), tolerance = 1e-10)

  # Test edge case
  # ---- Case: n == 0
  expect_equal(js(1, 0), 0)

  # Test error case
  # ---- Case: non-numeric 'l'
  expect_error(js("1", 1), "Inputs must be numeric or complex vectors.")
  # ---- Case: non-numeric 'n'
  expect_error(js(1, "1"), "Inputs must be numeric or complex vectors.")
  # ---- Case: non-numeric 'l' and 'n'
  expect_error(js("1", "1"), "Inputs must be numeric or complex vectors.")

  # Test jsd
  # ---- Case: Single value
  expect_equal(
    jsd(1, 1),
    js(0, 1) - (2) / 1 * js(1, 1)
  )
  # ---- Case: Vector[l]
  expect_equal(
    jsd(c(1, 2, 3), 1),
    js(c(1, 2, 3) - 1, 1) - (c(1, 2, 3) + 1) / 1 * js(c(1, 2, 3), 1)
  )
  # ---- Case: Vector[n]
  expect_equal(
    jsd(1, c(1, 2, 3)),
    js(0, c(1, 2, 3)) - (2) / c(1, 2, 3) * js(1, c(1, 2, 3))
  )

  # Test jsdd
  expect_equal(jsdd(1, 1), -0.177098575)

  # Test ys (spherical Bessel function of the second kind)
  # ys should be related to yc by: ys(l, n) = sqrt(pi/(2*n)) * yc(l + 0.5, n)
  expected_ys <- sqrt(pi / (2 * n)) * yc(l + 0.5, n)
  expect_equal(ys(l, n), Re(expected_ys), tolerance = 1e-10)

  # Test edge-cases
  # ---- Case: n < 0
  expect_equal(ys(1, -1), -sin(1) - cos(1))
  # --- Case: n == 0
  expect_equal(ys(1, 0), -Inf)

  # Test error case
  # ---- Case: non-numeric 'l'
  expect_error(ys("1", 1), "Inputs must be numeric or complex vectors.")
  # ---- Case: non-numeric 'n'
  expect_error(ys(1, "1"), "Inputs must be numeric or complex vectors.")
  # ---- Case: non-numeric 'l' and 'n'
  expect_error(ys("1", "1"), "Inputs must be numeric or complex vectors.")

  # Test ysd
  expect_equal(ysd(1, 1), 2 * sin(1) + cos(1))

  # Test ysdd
  expect_equal(ysdd(1, 1), -5 * sin(1) - 3 * cos(1))

  # Test hs (spherical Hankel function)
  # hs should equal js + 1i * ys
  expected_hs <- js(l, n) + 1i * ys(l, n)
  expect_equal(hs(l, n), expected_hs, tolerance = 1e-10)

  # Test hsd
  expect_equal(hsd(1, 1), (2 + 1i) * exp(1i))

  # Test higher order derivative functions (kth derivatives)
  expect_equal(jsdk(1, 1, 1), jsd(1, 1))
  expect_equal(jsdk(1, 1, 2), jsdd(1, 1))
  expect_equal(jsdk(0, 1, 1), jsd(0, 1))
  expect_equal(jsdk(0, 1, 2), jsdd(0, 1))
})

test_that("Spherical Bessel functions support complex arguments", {
  z <- 1 + 0.5i

  expect_equal(
    js(0, z),
    sin(z) / z,
    tolerance = 1e-12
  )
  expect_equal(
    js(1, z),
    sin(z) / z^2 - cos(z) / z,
    tolerance = 1e-12
  )
  expect_equal(
    ys(0, z),
    -cos(z) / z,
    tolerance = 1e-12
  )
  expect_equal(
    ys(1, z),
    -cos(z) / z^2 - sin(z) / z,
    tolerance = 1e-12
  )
  expect_equal(
    hs(1, z),
    js(1, z) + 1i * ys(1, z),
    tolerance = 1e-12
  )
  expect_equal(
    jsd(0, z),
    -js(1, z),
    tolerance = 1e-12
  )
  expect_equal(
    ysd(0, z),
    -ys(1, z),
    tolerance = 1e-12
  )
  expect_equal(
    hsd(0, z),
    -hs(1, z),
    tolerance = 1e-12
  )
})

test_that("Hankel functions work correctly", {
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

  # Test hcdd
  expect_equal(hcdd(1, 1), hc(1, 1) - hc(0, 1))
  # ---- Case: vector[n]
  expect_equal(
    hcdd(1, c(1, 2, 3)),
    c(
      hc(1, 1) - hc(0, 1),
      0.5 * (-hc(0, 2) - hc(1, 2)),
      1 / 9 * (-3 * hc(0, 3) - 7 * hc(1, 3))
    )
  )
  # ---- Case: matrix[n]
  expect_equal(
    as.vector(hcdd(1, matrix(c(1, 2, 3)))),
    c(
      hc(1, 1) - hc(0, 1),
      0.5 * (-hc(0, 2) - hc(1, 2)),
      1 / 9 * (-3 * hc(0, 3) - 7 * hc(1, 3))
    )
  )
  # ---- Case: n == 0
  expect_true(
    is.na(Re(hcdd(1, 0))) & is.na(Im(hcdd(1, 0)))
  )


  # Test hcdk
  # ---- Case: first derivative
  expect_equal(hcdk(1, 1, 1), hcd(1, 1))
  # ---- Case: second derivative
  expect_equal(hcdk(1, 1, 2), hcdd(1, 1))
  # ---- Case: vector[n]
  expect_equal(hcdk(1, c(1, 2, 3), 1), hcd(1, c(1, 2, 3)))
  # ---- Case: matrix[n]
  expect_equal(
    hcdk(1, matrix(c(1, 2, 3)), 1),
    hcd(1, matrix(c(1, 2, 3)))
  )
})
