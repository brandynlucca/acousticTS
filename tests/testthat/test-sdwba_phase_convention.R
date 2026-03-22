test_that("sdwba_stochastic_summary follows the paper phase convention", {
  segment_integrals <- matrix(c(1 + 0i, 1 + 0i), nrow = 1)
  phase_sd <- 0.5
  n_iterations <- 5

  set.seed(123)
  out <- acousticTS:::sdwba_stochastic_summary(
    segment_integrals = segment_integrals,
    phase_sd = phase_sd,
    n_iterations = n_iterations
  )

  set.seed(123)
  rng_1 <- stats::rnorm(n_iterations, mean = 0, sd = 1)
  rng_2 <- stats::rnorm(n_iterations, mean = 0, sd = 1)
  manual_f_bs <- exp(1i * rng_1 * phase_sd) + exp(1i * rng_2 * phase_sd)
  manual_sigma <- Mod(manual_f_bs)^2

  expect_equal(out$f_bs, mean(manual_f_bs))
  expect_equal(out$sigma_bs, mean(manual_sigma))
  expect_equal(out$TS_mean, acousticTS::db(mean(manual_sigma)))
})
