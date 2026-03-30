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

test_that("Smn_cpp covers shared-order, pairwise, outer-product, and quad paths", {
  eta_vec <- c(-0.5, 0, 0.5)

  shared <- acousticTS:::Smn_cpp(
    m = 1L,
    n = 1:3,
    c = 1.5,
    arg = eta_vec,
    normalize = TRUE,
    precision = "double"
  )
  expect_equal(dim(shared$value), c(3L, 3L))
  expect_equal(
    shared$value[2, ],
    vapply(
      eta_vec,
      function(eta_i) Smn(1, 2, 1.5, eta_i, normalize = TRUE)$value,
      numeric(1)
    ),
    tolerance = 1e-12
  )
  expect_equal(
    shared$derivative[3, ],
    vapply(
      eta_vec,
      function(eta_i) Smn(1, 3, 1.5, eta_i, normalize = TRUE)$derivative,
      numeric(1)
    ),
    tolerance = 1e-12
  )

  pairwise <- acousticTS:::Smn_cpp(
    m = c(0L, 1L),
    n = c(1L, 2L),
    c = 1.5,
    arg = 0.25,
    normalize = FALSE,
    precision = "double"
  )
  expect_equal(
    pairwise$value,
    c(
      Smn(0, 1, 1.5, 0.25)$value,
      Smn(1, 2, 1.5, 0.25)$value
    ),
    tolerance = 1e-12
  )
  expect_equal(
    pairwise$derivative,
    c(
      Smn(0, 1, 1.5, 0.25)$derivative,
      Smn(1, 2, 1.5, 0.25)$derivative
    ),
    tolerance = 1e-12
  )

  outer <- acousticTS:::Smn_cpp(
    m = c(0L, 1L),
    n = 1:3,
    c = 1.5,
    arg = c(-0.5, 0.5),
    normalize = FALSE,
    precision = "double"
  )
  expect_type(outer, "list")
  expect_length(outer, 2L)
  expect_equal(dim(outer[[1]]$value), c(2L, 3L))
  expect_equal(outer[[1]]$value[1, 1], Smn(0, 1, 1.5, -0.5)$value, tolerance = 1e-12)
  expect_equal(outer[[2]]$value[2, 2], Smn(1, 2, 1.5, 0.5)$value, tolerance = 1e-12)

  shared_double <- acousticTS:::Smn_cpp(
    m = 1L,
    n = 1:3,
    c = 1.5,
    arg = c(-0.5, 0.5),
    normalize = FALSE,
    precision = "double"
  )
  shared_quad_scalar <- Smn(1, 2, 1.5, 0.5, normalize = FALSE, precision = "quad")
  shared_double_scalar <- Smn(1, 2, 1.5, 0.5, normalize = FALSE, precision = "double")
  expect_equal(shared_quad_scalar$value, shared_double_scalar$value, tolerance = 1e-10)
  expect_equal(shared_quad_scalar$derivative, shared_double_scalar$derivative, tolerance = 1e-10)

  expect_error(
    acousticTS:::Smn_cpp(integer(0), 1L, 1.5, 0.5, FALSE, "double"),
    "must have at least one element"
  )
  expect_error(
    acousticTS:::Smn_cpp(0L, 1L, 1.5, numeric(0), FALSE, "double"),
    "must have at least one element"
  )
  expect_error(
    acousticTS:::Smn_cpp(-1L, 1L, 1.5, 0.5, FALSE, "double"),
    "must be >= 0"
  )
  expect_error(
    acousticTS:::Smn_cpp(0L, 1L, 1.5, 1.1, FALSE, "double"),
    "\\|eta\\| must be <= 1"
  )
  expect_error(
    acousticTS:::Smn_cpp(c(1L, 2L), c(0L, 1L), 1.5, 0.5, FALSE, "double"),
    "pairwise evaluation"
  )
  expect_error(
    acousticTS:::Smn_cpp(0L, 1L, 1.5, 0.5, FALSE, "half"),
    "'precision' must be 'double' or 'quad'"
  )
})

test_that("Rmn_cpp covers real and complex layout branches", {
  kind1 <- acousticTS:::Rmn_cpp(
    m = 0L,
    n = 1:3,
    c = 1.5,
    x1 = 1.2,
    kind = 1L,
    precision = "double"
  )
  kind2 <- acousticTS:::Rmn_cpp(
    m = 0L,
    n = 1:3,
    c = 1.5,
    x1 = 1.2,
    kind = 2L,
    precision = "double"
  )
  kind3 <- acousticTS:::Rmn_cpp(
    m = 0L,
    n = 1:3,
    c = 1.5,
    x1 = 1.2,
    kind = 3L,
    precision = "double"
  )
  kind4 <- acousticTS:::Rmn_cpp(
    m = c(0L, 1L),
    n = 1:3,
    c = 1.5,
    x1 = 1.2,
    kind = 4L,
    precision = "double"
  )

  expect_type(kind1$value, "double")
  expect_equal(length(kind1$value), 3L)
  expect_equal(kind2$value[2], Rmn(0, 2, 1.5, 1.2, kind = 2)$value, tolerance = 1e-12)
  expect_equal(kind3$value, kind1$value + 1i * kind2$value, tolerance = 1e-10)
  expect_equal(kind3$derivative, kind1$derivative + 1i * kind2$derivative, tolerance = 1e-10)
  expect_true(is.matrix(kind4$value))
  expect_equal(dim(kind4$value), c(2L, 3L))
  expect_equal(kind4$value[1, 2], Conj(kind3$value[2]), tolerance = 1e-10)

  pairwise_complex <- acousticTS:::Rmn_cpp(
    m = c(0L, 1L),
    n = c(1L, 2L),
    c = 1.5,
    x1 = 1.2,
    kind = 3L,
    precision = "double"
  )
  expect_equal(dim(pairwise_complex$value), c(2L, 2L))
  expect_equal(pairwise_complex$value[1, 1], Rmn(0, 1, 1.5, 1.2, kind = 3)$value, tolerance = 1e-10)
  expect_equal(pairwise_complex$value[2, 2], Rmn(1, 2, 1.5, 1.2, kind = 3)$value, tolerance = 1e-10)
  expect_equal(pairwise_complex$value[1, 2], 0 + 0i, tolerance = 1e-12)

  kind1_quad_scalar <- Rmn(
    m = 0L,
    n = 2L,
    c = 1.5,
    xi = 1.2,
    kind = 1L,
    precision = "quad"
  )
  kind1_double_scalar <- Rmn(
    m = 0L,
    n = 2L,
    c = 1.5,
    xi = 1.2,
    kind = 1L,
    precision = "double"
  )
  expect_equal(kind1_quad_scalar$value, kind1_double_scalar$value, tolerance = 1e-10)
  expect_equal(kind1_quad_scalar$derivative, kind1_double_scalar$derivative, tolerance = 1e-10)

  invalid_outer <- acousticTS:::Rmn_cpp(
    m = c(1L, 2L),
    n = c(0L, 3L, 4L),
    c = 1.5,
    x1 = 1.2,
    kind = 1L,
    precision = "double"
  )
  expect_true(is.nan(invalid_outer$value[1, 1]))
  expect_error(
    acousticTS:::Rmn_cpp(c(1L, 2L), c(0L, 1L), 1.5, 1.2, 1L, "double"),
    "pairwise evaluation"
  )
  expect_error(
    acousticTS:::Rmn_cpp(0L, 1L, 1.5, 1.2, 1L, "half"),
    "'precision' must be 'double' or 'quad'"
  )
})
