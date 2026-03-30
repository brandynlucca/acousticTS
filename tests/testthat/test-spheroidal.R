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
