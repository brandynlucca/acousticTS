library(acousticTS)

test_that("canonicalize_shape fits all supported canonical target families", {
  shape_in <- arbitrary(
    x_body = c(0, 0.01, 0.02, 0.03, 0.04),
    radius_body = c(0, 0.004, 0.006, 0.004, 0)
  )

  expect_s4_class(
    canonicalize_shape(shape_in, to = "Sphere"),
    "Sphere"
  )
  expect_s4_class(
    canonicalize_shape(shape_in, to = "Cylinder"),
    "Cylinder"
  )
  expect_s4_class(
    canonicalize_shape(shape_in, to = "ProlateSpheroid"),
    "ProlateSpheroid"
  )
  expect_s4_class(
    canonicalize_shape(shape_in, to = "OblateSpheroid"),
    "OblateSpheroid"
  )
})

test_that("canonicalize_shape returns fit diagnostics when requested", {
  shape_in <- arbitrary(
    x_body = c(0, 0.01, 0.02, 0.03, 0.04),
    radius_body = c(0, 0.004, 0.006, 0.004, 0)
  )

  fit <- canonicalize_shape(
    shape_in,
    to = "Cylinder",
    diagnostics = TRUE
  )

  expect_named(fit, c("shape", "diagnostics"))
  expect_s4_class(fit$shape, "Cylinder")
  expect_true(is.list(fit$diagnostics$source))
  expect_true(is.list(fit$diagnostics$target))
  expect_true(is.list(fit$diagnostics$fit))
  expect_true(is.finite(fit$diagnostics$fit$radius_rmse))
  expect_true(is.finite(fit$diagnostics$fit$radius_nrmse))
})

test_that("length-volume canonicalization preserves cylinder length and volume", {
  shape_in <- arbitrary(
    x_body = c(0, 0.01, 0.02, 0.03, 0.04),
    radius_body = c(0.003, 0.004, 0.004, 0.004, 0.003)
  )

  fit <- canonicalize_shape(
    shape_in,
    to = "Cylinder",
    method = "length_volume",
    diagnostics = TRUE
  )

  expect_equal(fit$diagnostics$fit$length_ratio, 1, tolerance = 1e-12)
  expect_equal(fit$diagnostics$fit$volume_ratio, 1, tolerance = 5e-3)
})

test_that("length-volume spheroidal canonicalization errors when the target family is incompatible", {
  shape_in <- arbitrary(
    x_body = c(0, 0.01, 0.02, 0.03, 0.04),
    radius_body = c(0, 0.004, 0.006, 0.004, 0)
  )

  expect_error(
    canonicalize_shape(
      shape_in,
      to = "OblateSpheroid",
      method = "length_volume"
    ),
    "cannot preserve both length and volume while remaining an oblate spheroid"
  )
})

test_that("canonicalize_shape handles equal-area width-height arbitrary profiles", {
  shape_in <- arbitrary(
    x_body = c(0, 0.01, 0.02, 0.03),
    w_body = c(0.004, 0.008, 0.008, 0.004),
    zU_body = c(0.001, 0.002, 0.002, 0.001),
    zL_body = c(-0.001, -0.002, -0.002, -0.001)
  )

  fit <- canonicalize_shape(
    shape_in,
    to = "Sphere",
    diagnostics = TRUE
  )

  expect_s4_class(fit$shape, "Sphere")
  expect_true(is.finite(fit$diagnostics$fit$volume_ratio))
})
