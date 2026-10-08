test_that("angular batches preserve modal order and scalar values", {
  precisions <- c(
    "double",
    if (acousticTS:::.quad_precision_available()) "quad"
  )
  for (precision in precisions) {
    for (normalize in c(FALSE, TRUE)) {
      eta <- c(-0.7, 0, 0.4, 0.9)
      n <- c(4L, 1L, 3L, 1L)
      batch <- Smn(1, n, 2, eta, normalize, precision)
      for (field in c("value", "derivative")) {
        expected <- t(vapply(
          n,
          function(ni) {
            vapply(
              eta,
              function(e) Smn(1, ni, 2, e, normalize, precision)[[field]],
              numeric(1)
            )
          },
          numeric(length(eta))
        ))
        expect_equal(batch[[field]], expected, tolerance = 1e-10)
        expect_equal(
          Smn(1, n, 2, eta[1], normalize, precision)[[field]],
          expected[, 1],
          tolerance = 1e-10
        )
      }
      m <- c(1L, 0L, 1L, 0L)
      pair <- Smn(m, n, 2, 0.4, normalize, precision)
      pair_angles <- Smn(m, n, 2, eta, normalize, precision)
      for (field in c("value", "derivative")) {
        expect_equal(
          pair[[field]],
          vapply(
            seq_along(m),
            function(i) {
              Smn(m[i], n[i], 2, 0.4, normalize, precision)[[field]]
            },
            numeric(1)
          ),
          tolerance = 1e-10
        )
        expected_angles <- t(vapply(
          seq_along(m),
          function(i) {
            vapply(
              eta,
              function(e) Smn(m[i], n[i], 2, e, normalize, precision)[[field]],
              numeric(1)
            )
          },
          numeric(length(eta))
        ))
        expect_equal(pair_angles[[field]], expected_angles, tolerance = 1e-10)
      }
      m <- c(0L, 2L)
      n <- c(0L, 2L, 4L)
      outer <- Smn(m, n, 2, eta, normalize, precision)
      expect_length(outer, length(eta))
      for (k in seq_along(eta)) {
        scalar_eta <- Smn(m, n, 2, eta[k], normalize, precision)
        for (field in c("value", "derivative")) {
          expect_equal(
            outer[[k]][[field]],
            scalar_eta[[field]],
            tolerance = 1e-10
          )
          expect_true(is.nan(scalar_eta[[field]][2, 1]))
          for (i in seq_along(m)) {
            for (j in seq_along(n)) {
              if (n[j] >= m[i]) {
                expect_equal(
                  scalar_eta[[field]][i, j],
                  Smn(
                    m[i],
                    n[j],
                    2,
                    eta[k],
                    normalize,
                    precision
                  )[[field]],
                  tolerance = 1e-10
                )
              }
            }
          }
        }
      }
    }
  }
  expect_error(Smn(integer(), 1, 1, 0), "at least one")
  expect_error(Smn(0, integer(), 1, 0), "at least one")
  expect_error(Smn(0, 1, 1, numeric()), "at least one")
  expect_error(Smn(-1, 1, 1, 0), ">= 0")
  expect_error(Smn(0, -1, 1, 0), ">= 0")
  expect_error(Smn(0, 1, 1, 1.1), "angular prolate domain")
  expect_error(Smn(2, 1, 1, 0), "'n' must be >= 'm'", fixed = TRUE)
  expect_error(Smn(2, 0:3, 1, 0), "All 'n' values")
  expect_error(Smn(c(1, 2), c(2, 1), 1, 0), "pairwise")
})

test_that("radial batches agree with scalar modes and outgoing conventions", {
  precisions <- c(
    "double",
    if (acousticTS:::.quad_precision_available()) "quad"
  )
  for (precision in precisions) {
    for (kind in 1:4) {
      for (m in list(0L, c(0L, 1L), c(0L, 1L, 0L))) {
        n <- c(2L, 3L, 1L)
        result <- Rmn(m, n, 2, 1.7, kind, precision)
        for (field in c("value", "derivative")) {
          if (length(m) == 1) {
            expected <- vapply(
              n,
              function(ni) {
                Rmn(
                  m,
                  ni,
                  2,
                  1.7,
                  kind,
                  precision
                )[[field]]
              },
              complex(1)
            )
            expect_equal(result[[field]] + 0i, expected, tolerance = 1e-08)
          } else if (length(m) == length(n)) {
            expected <- vapply(
              seq_along(m),
              function(i) {
                Rmn(m[i], n[i], 2, 1.7, kind, precision)[[field]]
              },
              complex(1)
            )
            expect_equal(
              diag(result[[field]]) + 0i,
              expected,
              tolerance = 1e-08
            )
            expect_equal(
              result[[field]][row(result[[field]]) != col(result[[field]])],
              rep(if (kind <= 2) 0 else 0 + 0i, 6)
            )
          } else {
            expected <- t(vapply(
              m,
              function(mi) {
                vapply(
                  n,
                  function(ni) Rmn(mi, ni, 2, 1.7, kind, precision)[[field]],
                  complex(1)
                )
              },
              complex(length(n))
            ))
            expect_equal(result[[field]] + 0i, expected, tolerance = 1e-08)
          }
        }
      }
      expect_error(
        Rmn(2, 1, 2, 1.7, kind, precision),
        "'n' must be >= 'm'",
        fixed = TRUE
      )
      expect_error(Rmn(2, 0:3, 2, 1.7, kind, precision), "All 'n' values")
      expect_error(
        Rmn(c(1, 2), c(2, 1), 2, 1.7, kind, precision),
        "pairwise"
      )
    }
    for (field in c("value", "derivative")) {
      r1 <- Rmn(0, 0:3, 2, 1.7, 1, precision)[[field]]
      r2 <- Rmn(0, 0:3, 2, 1.7, 2, precision)[[field]]
      expect_equal(
        Rmn(0, 0:3, 2, 1.7, 3, precision)[[field]],
        r1 + 1i * r2,
        tolerance = 1e-09
      )
      expect_equal(
        Rmn(0, 0:3, 2, 1.7, 4, precision)[[field]],
        r1 - 1i * r2,
        tolerance = 1e-09
      )
    }
  }
})

test_that("radial Wronskians hold across size and focal distances", {
  # W[R1, R2] = 1 / (c * (xi^2 - 1)), independent of m and n.
  cases <- expand.grid(
    c = c(0.1, 1, 10, 50),
    xi = c(1.001, 1.1, 2, 5),
    m = c(0L, 1L, 3L)
  )
  for (i in seq_len(nrow(cases))) {
    case <- cases[i, ]
    n <- case$m + c(0L, 1L, 4L)
    r1 <- Rmn(case$m, n, case$c, case$xi, kind = 1)
    r2 <- Rmn(case$m, n, case$c, case$xi, kind = 2)
    wronskian <- r1$value * r2$derivative - r1$derivative * r2$value
    expect_equal(
      wronskian * case$c * (case$xi^2 - 1),
      rep(1, length(n)),
      tolerance = 2e-06,
      info = paste("case", i)
    )
  }
})

test_that("normalized angular solutions obey parity and orthonormality", {
  quad <- gauss_legendre(96, -1, 1)
  for (c in c(0.1, 1, 10, 50)) {
    for (m in c(0L, 1L, 3L)) {
      n <- m + 0:4
      angular <- Smn(m, n, c, quad$nodes, normalize = TRUE)
      gram <- (angular$value * rep(quad$weights, each = length(n))) %*%
        t(angular$value)
      expect_equal(gram, diag(length(n)), tolerance = 1e-09)
      positive <- Smn(m, n, c, 0.37, normalize = TRUE)
      negative <- Smn(m, n, c, -0.37, normalize = TRUE)
      expect_equal(
        negative$value,
        (-1)^(n - m) * positive$value,
        tolerance = 1e-10
      )
      expect_equal(
        negative$derivative,
        (-1)^(n - m + 1) * positive$derivative,
        tolerance = 1e-10
      )
    }
  }
})

test_that("spheroidal focal and polar limits match nearby regular values", {
  for (c in c(0.1, 2, 20)) {
    for (m in 0:3) {
      n <- m + 0:3
      polar <- Smn(m, n, c, 1)
      near <- Smn(m, n, c, 1 - 1e-10)
      if (m == 0) {
        expect_equal(polar$value, near$value, tolerance = 1e-07)
        h <- 1e-05
        derivative <- (
          3 * polar$value - 4 * Smn(m, n, c, 1 - h)$value +
            Smn(m, n, c, 1 - 2 * h)$value
        ) / (2 * h)
        expect_lt(max(abs(polar$derivative - derivative)), 1e-06)
      } else {
        expect_equal(polar$value, rep(0, length(n)), tolerance = 1e-12)
      }
      focal <- Rmn(m, n, c, 1, kind = 1)
      if (m == 0) {
        near <- Rmn(m, n, c, 1 + 1e-10, kind = 1)
        expect_equal(focal$value, near$value, tolerance = 1e-07)
        h <- 1e-06
        derivative <- (
          -3 * focal$value + 4 * Rmn(m, n, c, 1 + h, kind = 1)$value -
            Rmn(m, n, c, 1 + 2 * h, kind = 1)$value
        ) / (2 * h)
        expect_lt(max(abs(focal$derivative - derivative)), 1e-06)
      } else {
        expect_equal(focal$value, rep(0, length(n)), tolerance = 1e-12)
      }
    }
  }
})

test_that("radial Wronskians remain normalized for large size parameters", {
  cases <- expand.grid(
    c = c(100, 300),
    xi = c(1 + 1e-08, 1.01, 1.5),
    m = c(
      0L,
      5L,
      20L
    )
  )
  for (i in seq_len(nrow(cases))) {
    case <- cases[i, ]
    n <- case$m + c(0L, 1L, 8L)
    r1 <- Rmn(case$m, n, case$c, case$xi, kind = 1)
    r2 <- Rmn(case$m, n, case$c, case$xi, kind = 2)
    wronskian <- r1$value * r2$derivative - r1$derivative * r2$value
    expect_equal(
      wronskian * case$c * (case$xi^2 - 1),
      rep(1, length(n)),
      tolerance = 1e-05,
      info = paste("case", i)
    )
  }
})

test_that("radial solutions remain independent through the modal transition", {
  for (c in c(50, 100, 300)) {
    for (xi in c(1.01, 1.1, 1.5, 2, 5)) {
      n <- 0:ceiling(1.5 * c)
      r1 <- Rmn(0, n, c, xi, kind = 1)
      r2 <- Rmn(0, n, c, xi, kind = 2)
      wronskian <- r1$value * r2$derivative - r1$derivative * r2$value
      normalized_wronskian <- wronskian * c * (xi^2 - 1)
      expect_equal(
        normalized_wronskian,
        rep(1, length(n)),
        tolerance = 1e-06,
        info = paste("c =", c, "xi =", xi)
      )
    }
  }
})

test_that("radial solutions retain accuracy for high azimuthal orders", {
  cases <- rbind(
    expand.grid(
      c = c(50, 200, 500),
      xi = c(1.05, 1.25, 2),
      m = c(
        25L,
        50L
      )
    ),
    expand.grid(c = c(50, 200, 500), xi = c(1.25, 2), m = 100L)
  )
  for (i in seq_len(nrow(cases))) {
    case <- cases[i, ]
    n <- case$m + 0:ceiling(case$c)
    r1 <- Rmn(case$m, n, case$c, case$xi, kind = 1)
    r2 <- Rmn(case$m, n, case$c, case$xi, kind = 2)
    wronskian <- r1$value * r2$derivative - r1$derivative * r2$value
    normalized_wronskian <- wronskian * case$c * (case$xi^2 - 1)
    expect_equal(
      normalized_wronskian,
      rep(1, length(n)),
      tolerance = 1e-06,
      info = paste("case", i)
    )
  }
})

test_that("Spheroidal wave functions work correctly", {
  # Angular wave function, Smn
  expect_equal(
    Smn(2, 3, 1, 0.5)$value,
    5.650368053851631
  )
  expect_equal(
    Smn(2, 3, 1, 0.5)$derivative,
    3.454326321112444
  )
  expect_equal(
    Smn(2, 3, 1, 0.0)$value,
    0
  )
  expect_equal(
    Smn(2, 3, 1, 0.0)$derivative,
    15.2772229786631266542735101
  )
  expect_equal(
    Smn(0, 3, 1, 0.5)$value,
    -0.4302211279618986
  )
  expect_equal(
    Smn(0, 3, 1, 0.5)$derivative,
    0.4316245260227451
  )
})

test_that("Spheroidal wrappers validate arguments and expose radial kinds", {
  expect_error(
    Smn(0, 1.5, 1, 0.5),
    "'n' must be a real integer"
  )
  expect_error(
    Smn("bad", 1, 1, 0.5),
    "'m' must be a real integer"
  )
  expect_error(
    Smn(0, 1, 1, "bad"),
    "'eta' must be a real number"
  )
  expect_error(
    Smn(0, 1, c(1, 2), 0.5),
    "'c' must be a single, real number."
  )
  expect_error(
    Smn(0, 1, 1, 0.5, precision = "half"),
    "'precision' must either be 'double' \\(default\\) or 'quad'"
  )

  radial_first <- Rmn(m = 0, n = 1, c = 1, xi = 1.5, kind = 1)
  radial_third <- Rmn(m = 0, n = 1, c = 1, xi = 1.5, kind = 3)

  expect_true(is.list(radial_first))
  expect_true(all(c("value", "derivative") %in% names(radial_first)))
  expect_true(is.complex(radial_third$value))
  expect_error(
    Rmn(m = 0, n = Inf, c = 1, xi = 1.5),
    "'n' must be a real integer"
  )
  expect_error(
    Rmn(m = "bad", n = 1, c = 1, xi = 1.5),
    "'m' must be a real integer"
  )
  expect_error(
    Rmn(m = 0, n = 1, c = c(1, 2), xi = 1.5),
    "'c' must be a single, real number."
  )
  expect_error(
    Rmn(m = 0, n = 1, c = 1, xi = c(1.5, 2)),
    "'xi' must be a single, real number."
  )
  expect_error(
    Rmn(m = 0, n = 1, c = 1, xi = 1.5, precision = "half"),
    "'precision' must either be 'double' \\(default\\) or 'quad'"
  )
})
