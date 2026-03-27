library(acousticTS)

.nonuniform_profile <- function(length_body, radius_body, n_points = 30L) {
  ang <- seq(-pi / 2, pi / 2, length.out = n_points)
  data.frame(
    x = (length_body / 2) * (sin(ang) + 1),
    w = 2 * radius_body * cos(ang),
    zU = radius_body * cos(ang),
    zL = -radius_body * cos(ang)
  )
}

test_that("KRM vectorized FLS matches per-frequency evaluation", {
  freqs <- c(12e3, 166e3, 314e3)

  for (shape in c("sphere", "prolate")) {
    prof <- switch(shape,
      sphere = .nonuniform_profile(0.02, 0.01),
      prolate = .nonuniform_profile(0.14, 0.01)
    )

    obj <- fls_generate(
      arbitrary(
        x_body = prof$x,
        w_body = prof$w,
        zU_body = prof$zU,
        zL_body = prof$zL
      ),
      density_body = 1028.9,
      sound_speed_body = 1480.3,
      theta_body = pi / 2
    )

    vec <- target_strength(
      obj,
      frequency = freqs,
      model = "KRM",
      density_sw = 1026.8,
      sound_speed_sw = 1477.3
    )

    loop_ts <- vapply(freqs, function(f) {
      target_strength(
        obj,
        frequency = f,
        model = "KRM",
        density_sw = 1026.8,
        sound_speed_sw = 1477.3
      )@model$KRM$TS
    }, numeric(1))

    expect_equal(vec@model$KRM$TS, loop_ts, tolerance = 1e-10)
  }
})

test_that("KRM vectorized SBF matches per-frequency evaluation", {
  freqs <- c(38e3, 120e3, 240e3)
  body <- .nonuniform_profile(0.14, 0.01, n_points = 30L)
  bladder <- .nonuniform_profile(0.06, 0.004, n_points = 20L)
  bladder$x <- bladder$x + 0.04

  obj <- sbf_generate(
    x_body = body$x,
    w_body = body$w,
    zU_body = body$zU,
    zL_body = body$zL,
    x_bladder = bladder$x,
    w_bladder = bladder$w,
    zU_bladder = bladder$zU,
    zL_bladder = bladder$zL,
    density_body = 1028.9,
    sound_speed_body = 1480.3,
    density_bladder = 1.24,
    sound_speed_bladder = 345,
    theta_body = pi / 2,
    theta_bladder = pi / 2
  )

  vec <- target_strength(
    obj,
    frequency = freqs,
    model = "KRM",
    density_sw = 1026.8,
    sound_speed_sw = 1477.3
  )

  loop_ts <- vapply(freqs, function(f) {
    target_strength(
      obj,
      frequency = f,
      model = "KRM",
      density_sw = 1026.8,
      sound_speed_sw = 1477.3
    )@model$KRM$TS
  }, numeric(1))

  expect_equal(vec@model$KRM$TS, loop_ts, tolerance = 1e-10)
})

test_that("KRM sardine body+bladder matches KRMr/NOAA reference frequencies", {
  data(sardine, package = "acousticTS")

  freqs <- c(12e3, 24e3, 60e3, 100e3, 120e3, 200e3, 300e3, 400e3)
  expected_ts <- c(
    -33.9816000823311,
    -35.9807232703388,
    -40.7713426107244,
    -46.0211785786631,
    -76.0749010872345,
    -39.1283077156313,
    -46.2737714255569,
    -42.8167934500733
  )

  out <- target_strength(
    sardine,
    frequency = freqs,
    model = "KRM",
    density_sw = 1030,
    sound_speed_sw = 1490
  )

  expect_equal(out@model$KRM$TS, expected_ts, tolerance = 1e-6)

  out_explicit <- target_strength(
    sardine,
    frequency = freqs,
    model = "KRM",
    density_sw = 1030,
    sound_speed_sw = 1490,
    krm_variant = "lowcontrast"
  )

  out_body <- target_strength(
    sardine,
    frequency = freqs,
    model = "KRM",
    density_sw = 1030,
    sound_speed_sw = 1490,
    krm_variant = "body_embedded"
  )

  expect_equal(out@model$KRM$TS, out_explicit@model$KRM$TS, tolerance = 1e-10)
  expect_true(any(abs(out@model$KRM$TS - out_body@model$KRM$TS) > 1))
})
