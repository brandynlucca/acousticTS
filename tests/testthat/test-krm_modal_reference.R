.krm_reference_fish <- function(x = seq(-0.05, 0.05, length.out = 101),
                                theta = pi / 2, dx = 0, dz = 0,
                                body_dx = 0, centerline = rep(0, length(x))) {
  sbf_generate(
    body_shape = arbitrary(
      x_body = x + body_dx, w_body = rep(0.02, length(x)),
      zU_body = rep(0.01, length(x)), zL_body = rep(-0.01, length(x))
    ),
    bladder_shape = arbitrary(
      x_body = x + dx, w_body = rep(0.002, length(x)),
      zU_body = centerline + dz + 0.001,
      zL_body = centerline + dz - 0.001
    ),
    density_body = 1070, sound_speed_body = 1575,
    density_bladder = 1.2, sound_speed_bladder = 345,
    theta_body = theta, theta_bladder = theta
  )
}

.krm_reference_result <- function(object, frequency, variant = "lowcontrast") {
  target_strength(
    object, frequency, "krm", density_sw = 1025, sound_speed_sw = 1500,
    krm_variant = variant
  )@model$KRM
}

# Clay (1991), Eqs. (7)-(12); Clay (1992), Eqs. (15)-(16), (A7).
.krm_published_modal <- function(frequency, theta, sound_speed = 1500,
                                 density = 1025) {
  k <- 2 * pi * frequency / sound_speed
  q <- k * 0.001
  qi <- 2 * pi * frequency / 345 * 0.001
  gh <- (1.2 / density) * (345 / sound_speed)
  numerator <- -besselJ(qi, 1) * besselY(q, 0) +
    gh * besselY(q, 1) * besselJ(qi, 0)
  denominator <- -besselJ(qi, 1) * besselJ(q, 0) +
    gh * besselJ(q, 1) * besselJ(qi, 0)
  b0 <- -denominator / (denominator + 1i * numerator)
  delta <- k * 0.1 * cos(theta)
  sinc <- rep(1, length(delta))
  nonzero <- abs(delta) > 1e-12
  sinc[nonzero] <- sin(delta[nonzero]) / delta[nonzero]
  -1i * 0.1 / pi * sinc * b0
}

test_that("KRM low mode matches the published single-sinc cylinder", {
  frequency <- c(100, 1000, 35000)
  grids <- list(
    c(-0.05, 0.05), seq(-0.05, 0.05, length.out = 11),
    seq(-0.05, 0.05, length.out = 101),
    c(-0.05, -0.048, -0.02, 0.003, 0.049, 0.05)
  )
  for (variant in c("lowcontrast", "mixed", "body_embedded")) {
    c_medium <- if (variant == "body_embedded") 1575 else 1500
    rho_medium <- if (variant == "body_embedded") 1070 else 1025
    for (theta in c(75, 90, 105) * pi / 180) {
      reference <- .krm_published_modal(frequency, theta, c_medium, rho_medium)
      for (x in grids) {
        actual <- .krm_reference_result(.krm_reference_fish(x, theta),
                                        frequency, variant)
        expect_equal(as.vector(actual$f_bladder), reference, tolerance = 1e-12)
      }
    }
  }
})

test_that("KRM high-frequency bladder phase follows the published coordinates", {
  frequency <- c(38000, 120000)
  theta <- 75 * pi / 180
  for (variant in c("lowcontrast", "mixed", "body_embedded")) {
    c_high <- if (variant == "lowcontrast") 1500 else 1575
    k <- 2 * pi * frequency / c_high
    original <- .krm_reference_result(.krm_reference_fish(theta = theta),
                                      frequency, variant)
    for (shift in list(c(0.02, 0), c(0, 0.003), c(0.02, 0.003))) {
      moved <- .krm_reference_result(
        .krm_reference_fish(theta = theta, dx = shift[1], dz = shift[2]),
        frequency, variant
      )
      phase <- exp(-2i * k * (shift[1] * cos(theta) + shift[2] * sin(theta)))
      expect_equal(moved$f_bladder, original$f_bladder * phase, tolerance = 1e-12)
      expect_equal(moved$sigma_bladder, original$sigma_bladder, tolerance = 1e-12)
      expect_equal(moved$f_body, original$f_body, tolerance = 1e-12)
    }
  }
})

test_that("KRM high-frequency coherent sum preserves axial translation phase", {
  frequency <- c(38000, 120000)
  theta <- 75 * pi / 180
  original <- .krm_reference_result(.krm_reference_fish(theta = theta), frequency)
  moved <- .krm_reference_result(
    .krm_reference_fish(theta = theta, dx = 0.02, body_dx = 0.02), frequency
  )
  phase <- exp(-2i * 2 * pi * frequency / 1500 * 0.02 * cos(theta))
  expect_equal(moved$f_bs, original$f_bs * phase, tolerance = 1e-12)
  expect_equal(moved$TS, original$TS, tolerance = 1e-10)
})

test_that("KRM low mode retains one equivalent cylinder for a bent outline", {
  frequency <- c(100, 1000, 35000)
  theta <- 75 * pi / 180
  for (x in list(c(-0.05, 0, 0.05), seq(-0.05, 0.05, length.out = 101))) {
    centerline <- 0.004 * abs(x) / 0.05
    # The equivalent-cylinder formula uses volume and axial length only.
    original <- .krm_reference_result(
      .krm_reference_fish(x, theta, centerline = centerline), frequency
    )
    reference <- .krm_published_modal(frequency, theta)
    expect_equal(as.vector(original$f_bladder), reference, tolerance = 1e-12)
  }
})

test_that("KRM low-frequency offsets retain the published centered-cylinder expression", {
  frequency <- c(100, 1000, 35000)
  for (variant in c("lowcontrast", "mixed", "body_embedded")) {
    c_medium <- if (variant == "body_embedded") 1575 else 1500
    rho_medium <- if (variant == "body_embedded") 1070 else 1025
    for (theta in c(75, 90, 105) * pi / 180) {
      expected <- .krm_published_modal(frequency, theta, c_medium, rho_medium)
      for (shift in list(c(0.02, 0), c(0, 0.003), c(0.02, 0.003))) {
        result <- .krm_reference_result(
          .krm_reference_fish(theta = theta, dx = shift[1], dz = shift[2]),
          frequency, variant
        )
        expect_equal(as.vector(result$f_bladder), expected, tolerance = 1e-12)
        expect_equal(result$f_bs, result$f_body + result$f_bladder, tolerance = 1e-12)
      }
    }
  }
})
