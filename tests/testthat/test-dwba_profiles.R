library(acousticTS)

test_that("DWBA integrators validate geometry and preserve the axial limit", {
  profile <- rbind(x = c(0, 0.04), y = 0, z = 0, a = 0.003)
  k <- 2 * pi * c(12000, 38000) / 1500
  for (fun in list(
    acousticTS:::dwba_fbs_cpp,
    acousticTS:::dwba_segment_integrals_cpp
  )) {
    expect_error(fun(matrix(0, 3, 2), k, 0, 1, 0.02), "4 x n matrix")
    expect_error(fun(profile[, 1, drop = FALSE], k, 0, 1, 0.02), "at least two")
    expect_error(
      fun(profile, k, 0, 1, 0.02, subdivisions = 0),
      "positive integer"
    )
    expect_error(fun(profile, k, 0, 0, 0.02), "non-zero")
    expect_error(
      fun(profile, k, 0, 1, 0.02, rel_tol = 0, abs_tol = 0),
      "QUADPACK error"
    )
    expect_true(all(is.finite(fun(profile, k, 0, 1.03, 0.02))))
  }
  direct <- acousticTS:::dwba_fbs_cpp(profile, k, 0, 1.03, 0.02)
  segments <- acousticTS:::dwba_segment_integrals_cpp(
    profile,
    k,
    0,
    1.03,
    0.02
  )
  expect_equal(as.vector(direct), as.vector(segments), tolerance = 1e-10)
})

.make_dwba_arbitrary <- function(shape_name, n_segments) {
  n_nodes <- n_segments + 1

  if (shape_name == "sphere") {
    radius <- 0.01
    v <- seq(0, pi, length.out = n_nodes)
    x <- radius + radius * cos(v)
    a <- radius * sin(v)
  } else if (shape_name == "prolate_spheroid") {
    semi_major <- 0.07
    semi_minor <- 0.01
    v <- seq(0, pi, length.out = n_nodes)
    x <- semi_major + semi_major * cos(v)
    a <- semi_minor * sin(v)
  } else {
    x <- seq(0.07, 0, length.out = n_nodes)
    a <- rep(0.01, n_nodes)
  }

  fls_generate(
    shape = arbitrary(
      x_body = x,
      y_body = rep(0, length(x)),
      z_body = rep(0, length(x)),
      radius_body = a
    ),
    theta_body = pi / 2,
    density_body = 1028.9,
    sound_speed_body = 1480.3
  )
}

.make_dwba_canonical <- function(shape_name, n_segments) {
  switch(shape_name,
    sphere = fls_generate(
      shape = sphere(
        radius_body = 0.01,
        n_segments = n_segments
      ),
      theta_body = pi / 2,
      density_body = 1028.9,
      sound_speed_body = 1480.3
    ),
    prolate_spheroid = fls_generate(
      shape = prolate_spheroid(
        length_body = 0.14,
        radius_body = 0.01,
        n_segments = n_segments
      ),
      theta_body = pi / 2,
      density_body = 1028.9,
      sound_speed_body = 1480.3
    ),
    cylinder = fls_generate(
      shape = cylinder(
        length_body = 0.07,
        radius_body = 0.01,
        n_segments = n_segments
      ),
      theta_body = pi / 2,
      density_body = 1028.9,
      sound_speed_body = 1480.3
    )
  )
}

test_that(
  paste0(
    "Canonical DWBA shapes resolve to the same nodewise profiles as ",
    "equivalent arbitrary shapes"
  ),
  {
    freq <- c(120e3, 240e3)

    for (shape_name in c("sphere", "prolate_spheroid", "cylinder")) {
      n_segments <- if (shape_name == "cylinder") 120 else 100

      ts_canonical <- target_strength(
        .make_dwba_canonical(shape_name, n_segments),
        frequency = freq,
        model = "DWBA",
        sound_speed_sw = 1477.3,
        density_sw = 1026.8
      )@model$DWBA$TS

      ts_arbitrary <- target_strength(
        .make_dwba_arbitrary(shape_name, n_segments),
        frequency = freq,
        model = "DWBA",
        sound_speed_sw = 1477.3,
        density_sw = 1026.8
      )@model$DWBA$TS

      expect_equal(ts_canonical, ts_arbitrary, tolerance = 1e-10)
    }
  }
)

test_that("Canonical SDWBA matches the selected analytic profile", {
  cases <- list(
    deterministic = list(n_iterations = 1, phase_sd_init = 0),
    stochastic = list(n_iterations = 10, phase_sd_init = sqrt(2) / 32)
  )

  for (shape_name in c("sphere", "prolate_spheroid", "cylinder")) {
    body_length <- switch(
      shape_name,
      sphere = 0.02,
      prolate_spheroid = 0.14,
      cylinder = 0.07
    )
    n0 <- if (shape_name == "cylinder") {
      50
    } else {
      100
    }
    for (f in c(38e3, 240e3)) {
      n_selected <- max(
        n0,
        ceiling(n0 * f * body_length / (120e3 * 0.03835))
      )
      for (settings in cases) {
        run <- function(object) {
          set.seed(1)
          target_strength(
            object,
            frequency = f,
            model = "SDWBA",
            sound_speed_sw = 1477.4,
            density_sw = 1026.8,
            n_iterations = settings$n_iterations,
            n_segments_init = n0,
            phase_sd_init = settings$phase_sd_init,
            length_init = 0.03835,
            frequency_init = 120e3
          )@model$SDWBA
        }
        expect_equal(
          run(.make_dwba_canonical(shape_name, 1000)),
          run(.make_dwba_arbitrary(shape_name, n_selected)),
          tolerance = 1e-10,
          info = paste(shape_name, f, settings$phase_sd_init)
        )
      }
    }
  }
})

test_that("SDWBA analytic geometry and TS ignore constructor resolution", {
  frequency <- c(12e3, 38e3, 120e3, 400e3)
  for (shape_name in c("sphere", "prolate_spheroid", "cylinder")) {
    n0 <- if (shape_name == "cylinder") {
      50
    } else {
      100
    }
    objects <- lapply(c(20, 50, 51, 1000), function(n) {
      .make_dwba_canonical(shape_name, n)
    })
    run <- function(object) {
      set.seed(27)
      target_strength(
        object,
        frequency = frequency,
        model = "sdwba",
        n_segments_init = n0,
        phase_sd_init = sqrt(2) / 32,
        n_iterations = 10,
        sound_speed_sw = 1477.4,
        density_sw = 1026.8
      )
    }
    results <- lapply(objects, run)
    reference <- results[[1]]
    body_length <- reference@shape_parameters$length
    n_expected <- pmax(
      n0,
      ceiling(n0 * frequency * body_length / (120e3 * 0.03835))
    )
    groups <- reference@model_parameters$SDWBA$parameters
    for (group in groups) {
      expect_equal(group$n_segments, n_expected[group$input_indices[1]])
      expect_equal(ncol(group$body_params$rpos) - 1L, group$n_segments)
      expect_equal(
        group$meta_params$phase_sd,
        sqrt(2) / 32 * n0 * body_length / (group$n_segments * 0.03835)
      )
    }
    for (result in results[-1]) {
      expect_equal(result@model$SDWBA, reference@model$SDWBA, tolerance = 1e-12)
      actual <- result@model_parameters$SDWBA$parameters
      for (i in seq_along(groups)) {
        expect_equal(actual[[i]]$body_params$rpos, groups[[i]]$body_params$rpos)
        expect_equal(
          actual[[i]]$body_params$radius,
          groups[[i]]$body_params$radius
        )
        expect_equal(actual[[i]]$meta_params, groups[[i]]$meta_params)
      }
    }
  }
})

test_that("SDWBA analytic regeneration retains placement and axial direction", {
  for (shape_name in c("sphere", "prolate_spheroid", "cylinder")) {
    object <- .make_dwba_canonical(shape_name, 50)
    reference <- sdwba_resample(object, 51)
    object@body$rpos[1, ] <- object@body$rpos[1, ] + 0.123
    object@body$rpos <- object@body$rpos[,
      rev(seq_len(ncol(object@body$rpos))),
      drop = FALSE
    ]
    object@body$radius <- rev(object@body$radius)
    actual <- sdwba_resample(object, 51)
    expected <- reference@body$rpos[, 52:1]
    expected[1, ] <- expected[1, ] + 0.123
    expect_equal(actual@body$rpos, expected)
    expect_equal(actual@body$radius, rev(reference@body$radius))
    expect_equal(actual@shape_parameters$n_segments, 51)
  }
})

test_that("SDWBA preserves explicit curved geometry on a canonical shape", {
  object <- .make_dwba_canonical("cylinder", 50)
  z <- 0.002 * sin(seq(0, pi, length.out = 51))
  object@body$rpos[3, ] <- z
  object@body$rpos["zU", ] <- z + object@body$radius
  object@body$rpos["zL", ] <- z - object@body$radius
  explicit <- object
  explicit@shape_parameters$shape <- "Arbitrary"
  for (n in c(20, 50, 100)) {
    actual <- sdwba_resample(object, n)
    reference <- sdwba_resample(explicit, n)
    expect_equal(actual@body, reference@body)
    expect_equal(actual@shape_parameters$n_segments, n)
  }
})
