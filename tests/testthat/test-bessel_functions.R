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
  
  # Test matrix input for `yc`
  expect_equal(
    t(yc(1, matrix(c(1, 2, 3)))),
    t(besselY(c(1, 2, 3), 1))
  )
  
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
  
  # Test edge case
  # ---- Case: n == 0
  expect_equal(js(1, 0), 0)
  
  # Test error case
  # ---- Case: non-numeric 'l'
  expect_error(js("1", 1), "Inputs must be numeric vectors.")
  # ---- Case: non-numeric 'n'
  expect_error(js(1, "1"), "Inputs must be numeric vectors.")
  # ---- Case: non-numeric 'l' and 'n'
  expect_error(js("1", "1"), "Inputs must be numeric vectors.")
  
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
  expect_equal(ys(l, n), expected_ys, tolerance = 1e-10)
  
  # Test edge-cases
  # ---- Case: n < 0
  expect_equal(ys(1, -1), -sin(1) - cos(1))
  # --- Case: n == 0
  expect_equal(ys(1, 0), -Inf)
  
  # Test error case
  # ---- Case: non-numeric 'l'
  expect_error(ys("1", 1), "Inputs must be numeric vectors.")
  # ---- Case: non-numeric 'n'
  expect_error(ys(1, "1"), "Inputs must be numeric vectors.")
  # ---- Case: non-numeric 'l' and 'n'
  expect_error(ys("1", "1"), "Inputs must be numeric vectors.")
  
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

test_that(
  "Bessel function cache for Goodman and Stern (1962) model works correctly", {
    library(acousticTS)
    
    # Calculate ka matrix
    # ---- Mock variables
    frequency <- seq(1e3, 10e3, 1e3)
    sound_speed_sw <- 1500
    sound_speed_fl <- 3000
    sound_speed_longitudinal <- 4000
    sound_speed_transversal <- 5000
    radius_shell <- 2
    radius_fluid <- 1
    m_limit <- 2
    m <- c(0, 1, 2)
    # ---- Compute
    ka_matrix <- calculate_ka_matrix(
      frequency, sound_speed_sw, sound_speed_fl, sound_speed_longitudinal, 
      sound_speed_transversal, radius_shell, radius_fluid
    )
    # ---- Compute
    ka_matrix_m <- lapply( rownames( ka_matrix ) , function( ka ) {
      modal_matrix( ka_matrix[ ka , ] , m_limit )
    } )
    names( ka_matrix_m ) <- rownames( ka_matrix )
    
    # Compute the cached Bessel functions
    cached_bessel <- .calculate_bessel_cache(ka_matrix_m, m)
    
    # Test list entries
    expect_equal(
      names(cached_bessel), 
      c("k1a_shell", "kLa_shell", "kTa_shell", "k1a_fluid", "kTa_fluid", 
        "kLa_fluid", "k3a_fluid")
      )
    
    # Test sums
    sums <- unname(
      vapply(cached_bessel, function(x) sum(unlist(x)), complex(1))
    )
    expected <- c(0.3454705398047 + 0.0295538507i, 
                  0.06375339922440937 + 0i, 
                  -0.10692964921272691 + 0i, 
                  0.0 + 0i, 
                  -11.894443889439168 + 0i, 
                  -3.5817720601741008 + 0i, 
                  0.46363940290204136 + 0i)
    # ---- Case: Reals
    expect_equal(Re(sums), Re(expected))
    # ---- Case: Imaginaries
    expect_equal(Im(sums), Im(expected))
  })
