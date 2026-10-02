.sdwba_test_body <- function() {
  fls_generate(
    x_body = seq(0, 0.03835, length.out = 15),
    y_body = rep(0, 15), z_body = rep(0, 15),
    radius_body = rep(0.001, 15),
    g_body = 1.0357, h_body = 1.0279,
    radius_curvature_ratio = 5
  )
}

test_that("SDWBA phase scaling follows each group's member frequencies", {
  object <- .sdwba_test_body()
  frequencies <- list(c(38e3, 120e3, 200e3), c(200e3, 38e3, 120e3, 200e3))
  for (curved in c(FALSE, TRUE)) {
    initialize <- if (curved) sdwba_curved_initialize else sdwba_initialize
    slot <- if (curved) "SDWBA_curved" else "SDWBA"
    for (frequency in frequencies) {
      batch <- suppressWarnings(initialize(object, frequency))
      groups <- batch@model_parameters[[slot]]$parameters
      for (group in groups) {
        for (f in group$acoustics$frequency) {
          single <- suppressWarnings(initialize(object, f))
          reference <- single@model_parameters[[slot]]$parameters[[1]]
          expect_equal(group$meta_params$phase_sd, reference$meta_params$phase_sd)
          expect_equal(group$n_segments, reference$n_segments)
        }
      }
    }
  }
})

test_that("SDWBA retains input order and duplicate frequencies", {
  object <- .sdwba_test_body()
  frequency <- c(38e3, 200e3, 120e3, 38e3, 300e3, 200e3)
  for (model in c("sdwba", "sdwba_curved")) {
    slot <- if (model == "sdwba") "SDWBA" else "SDWBA_curved"
    run <- function(f) {
      suppressWarnings(target_strength(
        object, frequency = f, model = model, phase_sd_init = 0,
        n_iterations = 2, sound_speed_sw = 1456, density_sw = 1025
      ))@model[[slot]]
    }
    batch <- run(frequency)
    singles <- do.call(rbind, lapply(frequency, run))
    expect_equal(batch$frequency, frequency)
    for (column in c("f_bs", "sigma_bs", "TS", "TS_sd")) {
      expect_equal(batch[[column]], singles[[column]], tolerance = 1e-12)
    }
  }
})

test_that("SDWBA resampling preserves nodes at unchanged resolution", {
  x <- c(0, 0.004, 0.02, 0.03)
  for (indices in list(1:4, 4:1)) {
    object <- fls_generate(
      x_body = x[indices], y_body = c(0, 0.002, 0.001, 0)[indices],
      z_body = c(0, 0.003, 0.001, 0)[indices],
      radius_body = c(0.0005, 0.001, 0.002, 0.0003)[indices],
      g_body = 1.0357, h_body = 1.0279
    )
    resampled <- sdwba_resample(object, 3)
    expect_equal(resampled@body$rpos, object@body$rpos)
    expect_equal(resampled@body$radius, object@body$radius)
  }
})

test_that("SDWBA resampling linearly interpolates position and radius", {
  object <- fls_generate(
    x_body = c(0, 0.01, 0.03), y_body = c(0, 0.002, 0),
    z_body = c(0, 0.003, 0.001), radius_body = c(0.0005, 0.002, 0.001),
    g_body = 1.0357, h_body = 1.0279
  )
  resampled <- sdwba_resample(object, 6)
  expect_equal(resampled@body$rpos[1, ], seq(0, 0.03, length.out = 7))
  expect_equal(resampled@body$radius,
               c(0.0005, 0.00125, 0.002, 0.00175, 0.0015, 0.00125, 0.001))
  expect_equal(unname(resampled@body$rpos[2, ]),
               c(0, 0.001, 0.002, 0.0015, 0.001, 0.0005, 0))
  expect_equal(unname(resampled@body$rpos[3, ]),
               c(0, 0.0015, 0.003, 0.0025, 0.002, 0.0015, 0.001))
  expect_equal(resampled@body$rpos["zU", ],
               resampled@body$rpos["z", ] + resampled@body$radius)
  expect_equal(resampled@body$rpos["zL", ],
               resampled@body$rpos["z", ] - resampled@body$radius)
  reversed <- object
  reversed@body$rpos <- object@body$rpos[, 3:1]
  reversed@body$radius <- rev(object@body$radius)
  reversed <- sdwba_resample(reversed, 6)
  expect_equal(unname(reversed@body$rpos), unname(resampled@body$rpos[, 7:1]))
  expect_equal(reversed@body$radius, rev(resampled@body$radius))
  # Two endpoint nodes must remain endpoints when coarsening/refining.
  coarse <- sdwba_resample(object, 1)
  expect_equal(coarse@body$radius, object@body$radius[c(1, 3)])
  expect_equal(sdwba_resample(coarse, 6)@body$radius,
               seq(0.0005, 0.001, length.out = 7))
})

test_that("Zero-phase SDWBA agrees with DWBA for an unchanged cylinder", {
  object <- .sdwba_test_body()
  for (theta in c(pi / 2, pi / 3)) {
    object@body$theta <- theta
    for (f in c(38e3, 120e3, 200e3)) {
      result <- target_strength(
        object, frequency = f, model = c("dwba", "sdwba"),
        phase_sd_init = 0, n_iterations = 3,
        sound_speed_sw = 1456, density_sw = 1025
      )@model
      expect_equal(result$SDWBA$f_bs, result$DWBA$f_bs, tolerance = 1e-10)
      expect_equal(result$SDWBA$TS, result$DWBA$TS, tolerance = 1e-10)
      expect_equal(result$SDWBA$TS_sd, 0)
    }
  }
})

test_that("SDWBA reports sample dispersion in realization TS", {
  set.seed(20261001)
  result <- sdwba_stochastic_summary(matrix(c(1, 1), nrow = 1), 0.7, 10000)
  set.seed(20261001)
  amplitude <- exp(1i * 0.7 * rnorm(10000)) + exp(1i * 0.7 * rnorm(10000))
  expect_equal(result$TS_sd, sd(10 * log10(Mod(amplitude)^2)))
  expect_equal(result$sigma_bs, mean(Mod(amplitude)^2))
  expect_equal(result$TS_mean, 10 * log10(mean(Mod(amplitude)^2)))
  expect_equal(sdwba_stochastic_summary(matrix(c(1, 1), 1), 0, 3)$TS_sd, 0)
  expect_true(is.na(sdwba_stochastic_summary(matrix(1, 1), 0, 1)$TS_sd))
  expect_identical(sdwba_stochastic_summary(matrix(0, 1), 0, 3)$TS_sd, NA_real_)
})

test_that("SDWBA stochastic power agrees with the Gaussian expectation", {
  amplitudes <- c(1 + 0.5i, -0.3 + 0.7i, 0.2 - 0.1i)
  phase_sd <- 0.7
  n <- 100000
  set.seed(42)
  result <- sdwba_stochastic_summary(matrix(amplitudes, nrow = 1), phase_sd, n)
  set.seed(42)
  draws <- matrix(rnorm(n * length(amplitudes)), nrow = n)
  power <- Mod(rowSums(sweep(exp(1i * phase_sd * draws), 2, amplitudes, `*`)))^2
  expected <- exp(-phase_sd^2) * Mod(sum(amplitudes))^2 +
    (1 - exp(-phase_sd^2)) * sum(Mod(amplitudes)^2)
  # Six estimated standard errors bound sampling noise, not solver error.
  expect_lt(abs(result$sigma_bs - expected), 6 * sd(power) / sqrt(n))
})

test_that("SDWBA preserves bent profiles in the zero-phase limit", {
  object <- .sdwba_test_body()
  object@body$rpos[3, ] <- 0.002 * sin(seq(0, pi, length.out = 15))
  object@body$rpos["zU", ] <- object@body$rpos[3, ] + object@body$radius
  object@body$rpos["zL", ] <- object@body$rpos[3, ] - object@body$radius
  object@body$theta <- pi / 3
  result <- target_strength(
    object, frequency = 120e3, model = c("dwba", "sdwba"),
    phase_sd_init = 0, n_iterations = 3
  )@model
  expect_equal(result$SDWBA$f_bs, result$DWBA$f_bs, tolerance = 1e-10)
  expect_equal(result$SDWBA$TS, result$DWBA$TS, tolerance = 1e-10)
  expect_equal(result$SDWBA$TS_sd, 0)
})

test_that("SDWBA resampling retains length metadata and validates its grid", {
  object <- .sdwba_test_body()
  object@shape_parameters$length <- 0.035
  expect_equal(sdwba_resample(object, 24)@shape_parameters$length, 0.035)
  for (n in list(0, -1, 1.5, NA_real_, Inf, c(1, 2))) {
    expect_error(sdwba_resample(object, n), "positive integer")
  }
  object@body$rpos[1, 4] <- object@body$rpos[1, 2]
  expect_error(sdwba_resample(object, 24), "strictly monotonic")
})
