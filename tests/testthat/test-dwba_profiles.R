library(acousticTS)

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
  switch(
    shape_name,
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

test_that("Canonical DWBA shapes resolve to the same nodewise profiles as equivalent arbitrary shapes", {
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
})

test_that("Canonical SDWBA shapes match equivalent arbitrary shapes", {
  freq <- c(120e3, 240e3)

  cases <- list(
    deterministic = list(n_iterations = 1, phase_sd_init = 0),
    stochastic = list(n_iterations = 10, phase_sd_init = 0.044)
  )

  for (case_name in names(cases)) {
    settings <- cases[[case_name]]

    for (shape_name in c("sphere", "prolate_spheroid", "cylinder")) {
      n_segments <- if (shape_name == "cylinder") 120 else 100

      set.seed(1)
      ts_canonical <- target_strength(
        .make_dwba_canonical(shape_name, n_segments),
        frequency = freq,
        model = "SDWBA",
        sound_speed_sw = 1477.3,
        density_sw = 1026.8,
        n_iterations = settings$n_iterations,
        n_segments_init = 50,
        phase_sd_init = settings$phase_sd_init,
        length_init = 38.35e-3,
        frequency_init = 120e3
      )@model$SDWBA$TS

      set.seed(1)
      ts_arbitrary <- target_strength(
        .make_dwba_arbitrary(shape_name, n_segments),
        frequency = freq,
        model = "SDWBA",
        sound_speed_sw = 1477.3,
        density_sw = 1026.8,
        n_iterations = settings$n_iterations,
        n_segments_init = 50,
        phase_sd_init = settings$phase_sd_init,
        length_init = 38.35e-3,
        frequency_init = 120e3
      )@model$SDWBA$TS

      expect_equal(
        ts_canonical,
        ts_arbitrary,
        tolerance = 1e-10,
        info = paste(case_name, shape_name)
      )
    }
  }
})
