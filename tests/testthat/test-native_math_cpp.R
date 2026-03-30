library(acousticTS)

test_that("direct Bessel cpp exports cover matrix, pairwise, and k-zero branches", {
  z_pair <- as.complex(c(1, 2))
  nu_pair <- as.complex(c(0, 1))
  z_matrix <- as.complex(c(1, 2))
  nu_matrix <- as.complex(c(0, 1, 2))

  expect_length(acousticTS:::jc_cpp(z_pair, nu_pair), 2)
  expect_length(acousticTS:::yc_cpp(z_pair, nu_pair), 2)
  expect_length(acousticTS:::hc_cpp(z_pair, nu_pair), 2)

  jc_mat <- acousticTS:::jc_cpp(z_matrix, nu_matrix)
  yc_mat <- acousticTS:::yc_cpp(z_matrix, nu_matrix)
  hc_mat <- acousticTS:::hc_cpp(z_matrix, nu_matrix)
  expect_true(is.matrix(jc_mat))
  expect_true(is.matrix(yc_mat))
  expect_true(is.matrix(hc_mat))
  expect_equal(dim(jc_mat), c(3L, 2L))
  expect_equal(dim(yc_mat), c(3L, 2L))
  expect_equal(dim(hc_mat), c(3L, 2L))

  expect_equal(acousticTS:::jc_deriv_cpp(z_pair, nu_pair, 0L), acousticTS:::jc_cpp(z_pair, nu_pair))
  expect_equal(acousticTS:::yc_deriv_cpp(z_pair, nu_pair, 0L), acousticTS:::yc_cpp(z_pair, nu_pair))
  expect_equal(acousticTS:::hc_deriv_cpp(z_pair, nu_pair, 0L), acousticTS:::hc_cpp(z_pair, nu_pair))
  expect_true(is.matrix(acousticTS:::jc_deriv_cpp(z_matrix, nu_matrix, 3L)))
  expect_true(is.matrix(acousticTS:::yc_deriv_cpp(z_matrix, nu_matrix, 3L)))
  expect_true(is.matrix(acousticTS:::hc_deriv_cpp(z_matrix, nu_matrix, 3L)))

  l_pair <- 0:1
  l_matrix <- 0:2
  z_real_pair <- c(1, 2)
  z_real_matrix <- c(1, 2)

  expect_length(acousticTS:::js_cpp(l_pair, z_real_pair), 2)
  expect_length(acousticTS:::ys_cpp(l_pair, z_real_pair), 2)
  expect_length(acousticTS:::hs_cpp(l_pair, z_real_pair), 2)
  expect_true(is.matrix(acousticTS:::js_cpp(l_matrix, z_real_matrix)))
  expect_true(is.matrix(acousticTS:::ys_cpp(l_matrix, z_real_matrix)))
  expect_true(is.matrix(acousticTS:::hs_cpp(l_matrix, z_real_matrix)))

  expect_equal(acousticTS:::js_deriv_cpp(l_pair, z_real_pair, 0L), acousticTS:::js_cpp(l_pair, z_real_pair))
  expect_equal(acousticTS:::ys_deriv_cpp(l_pair, z_real_pair, 0L), acousticTS:::ys_cpp(l_pair, z_real_pair))
  expect_equal(acousticTS:::hs_deriv_cpp(l_pair, z_real_pair, 0L), acousticTS:::hs_cpp(l_pair, z_real_pair))
  expect_true(is.matrix(acousticTS:::js_deriv_cpp(l_matrix, z_real_matrix, 3L)))
  expect_true(is.matrix(acousticTS:::ys_deriv_cpp(l_matrix, z_real_matrix, 3L)))
  expect_true(is.matrix(acousticTS:::hs_deriv_cpp(l_matrix, z_real_matrix, 3L)))

  z_complex_pair <- as.complex(c(1 + 0.5i, 2 + 0.25i))
  z_complex_matrix <- as.complex(c(1 + 0.5i, 2 + 0.25i))

  expect_length(acousticTS:::js_complex_cpp(l_pair, z_complex_pair), 2)
  expect_length(acousticTS:::ys_complex_cpp(l_pair, z_complex_pair), 2)
  expect_length(acousticTS:::hs_complex_cpp(l_pair, z_complex_pair), 2)
  expect_true(is.matrix(acousticTS:::js_complex_cpp(l_matrix, z_complex_matrix)))
  expect_true(is.matrix(acousticTS:::ys_complex_cpp(l_matrix, z_complex_matrix)))
  expect_true(is.matrix(acousticTS:::hs_complex_cpp(l_matrix, z_complex_matrix)))

  expect_equal(
    acousticTS:::js_complex_deriv_cpp(l_pair, z_complex_pair, 0L),
    acousticTS:::js_complex_cpp(l_pair, z_complex_pair)
  )
  expect_equal(
    acousticTS:::ys_complex_deriv_cpp(l_pair, z_complex_pair, 0L),
    acousticTS:::ys_complex_cpp(l_pair, z_complex_pair)
  )
  expect_equal(
    acousticTS:::hs_complex_deriv_cpp(l_pair, z_complex_pair, 0L),
    acousticTS:::hs_complex_cpp(l_pair, z_complex_pair)
  )
  expect_true(is.matrix(acousticTS:::js_complex_deriv_cpp(l_matrix, z_complex_matrix, 3L)))
  expect_true(is.matrix(acousticTS:::ys_complex_deriv_cpp(l_matrix, z_complex_matrix, 3L)))
  expect_true(is.matrix(acousticTS:::hs_complex_deriv_cpp(l_matrix, z_complex_matrix, 3L)))
})

test_that("direct Legendre cpp exports cover integer, fractional, endpoint, and derivative branches", {
  pn <- acousticTS:::Pn_cpp(c(0, 2, 2.5), c(-0.5, 0.25, 1.5))
  expect_true(is.matrix(pn))
  expect_equal(dim(pn), c(3L, 3L))
  expect_equal(pn[1, ], rep(1, 3), tolerance = 1e-12)

  expect_error(
    acousticTS:::Pn_deriv_cpp(c(1, 2), c(0.1, 0.3), -1L),
    "Derivative order k must be non-negative"
  )
  expect_equal(
    acousticTS:::Pn_deriv_cpp(c(1, 2), c(0.1, 0.3), 0L),
    acousticTS:::Pn_cpp(c(1, 2), c(0.1, 0.3))
  )
  pn_deriv <- acousticTS:::Pn_deriv_cpp(c(2, 2.5), c(-1, 0.25, 1), 1L)
  expect_true(is.matrix(pn_deriv))
  expect_equal(dim(pn_deriv), c(2L, 3L))
  expect_warning(
    pn_frac_high <- acousticTS:::Pn_deriv_cpp(c(2.5), c(0.25), 3L),
    "Fractional order derivatives k>1 use finite differences"
  )
  expect_true(all(is.finite(pn_frac_high)))

  qn <- acousticTS:::Qn_cpp(c(2, 2.5), c(-0.5, 1, 2))
  expect_true(is.matrix(qn))
  expect_equal(dim(qn), c(2L, 3L))
  expect_true(all(is.infinite(Re(qn[, 2]))))
  expect_true(all(is.infinite(Im(qn[, 2]))))
  expect_true(any(Im(qn[, 3]) != 0))

  expect_error(
    acousticTS:::Qn_deriv_cpp(c(1, 2), c(0.1, 0.3), -1L),
    "Derivative order k must be non-negative"
  )
  expect_equal(
    acousticTS:::Qn_deriv_cpp(c(1, 2), c(0.1, 0.3), 0L),
    acousticTS:::Qn_cpp(c(1, 2), c(0.1, 0.3))
  )
  qn_deriv <- acousticTS:::Qn_deriv_cpp(c(2, 2.5), c(-0.5, 2), 2L)
  expect_true(is.matrix(qn_deriv))
  expect_equal(dim(qn_deriv), c(2L, 2L))
})

test_that("Goodman-Stern cpp exports cover old and new shell-boundary interfaces", {
  old_res <- acousticTS:::elastic_shell_boundary_conditions_old(
    k1a = c(0.7, 1.1),
    kLa_shell = c(0.9, 1.3),
    kTa_shell = c(0.6, 0.8),
    kLa_fluid = c(0.5, 0.7),
    kTa_fluid = c(0.4, 0.6),
    k3a_fluid = c(0.3, 0.5),
    m_vec = 0:2,
    lambda = 2.5,
    mu = 1.2,
    rho_ratio_sw = 0.9,
    rho_ratio_fl = 0.8
  )
  expect_true(is.matrix(old_res))
  expect_equal(dim(old_res), c(2L, 3L))

  ka_matrix <- rbind(
    k1a_shell = c(0.7, 1.1),
    kLa_shell = c(0.9, 1.3),
    kTa_shell = c(0.6, 0.8),
    unused_shell = c(0.2, 0.2),
    kTa_fluid = c(0.4, 0.6),
    kLa_fluid = c(0.5, 0.7),
    k3a_fluid = c(0.3, 0.5)
  )
  new_res <- acousticTS:::elastic_shell_boundary_conditions(
    ka_matrix = ka_matrix,
    m_limit = c(2L, 1L),
    lambda = 2.5,
    mu = 1.2,
    rho_ratio_sw = 0.9,
    rho_ratio_fl = 0.8
  )
  expect_true(is.matrix(new_res))
  expect_equal(dim(new_res), c(2L, 3L))
  expect_true(is.na(new_res[2, 3]))
})
