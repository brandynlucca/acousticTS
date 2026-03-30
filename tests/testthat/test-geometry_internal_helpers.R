library(acousticTS)

test_that("canonicalize_shape internal helpers cover validation and fitting branches", {
  shape_in <- arbitrary(
    x_body = c(0, 0.01, 0.02, 0.03, 0.04),
    radius_body = c(0, 0.004, 0.006, 0.004, 0)
  )
  source_metrics <- acousticTS:::.canonicalize_shape_metrics(shape_in)

  expect_error(canonicalize_shape("not_a_shape"), "'shape' must be a Shape object.")
  expect_error(
    canonicalize_shape(shape_in, n_segments = 1),
    "'n_segments' must be one integer greater than or equal to 2."
  )

  local({
    testthat::local_mocked_bindings(
      .canonicalize_shape_metrics = function(shape) {
        list(length = 0, volume = 1, max_radius = 1, n_segments = 10)
      },
      .package = "acousticTS"
    )
    expect_error(canonicalize_shape(shape_in), "Source shape must have positive length.")
  })

  local({
    testthat::local_mocked_bindings(
      .canonicalize_shape_metrics = function(shape) {
        list(length = 1, volume = 0, max_radius = 1, n_segments = 10)
      },
      .package = "acousticTS"
    )
    expect_error(
      canonicalize_shape(shape_in),
      "Source shape must have positive enclosed volume."
    )
  })

  local({
    testthat::local_mocked_bindings(
      .canonicalize_shape_metrics = function(shape) {
        list(length = 1, volume = 1, max_radius = 0, n_segments = 10)
      },
      .package = "acousticTS"
    )
    expect_error(
      canonicalize_shape(shape_in),
      "Source shape must have positive radius somewhere along the profile."
    )
  })

  expect_error(
    acousticTS:::.canonicalize_shape_fit_volume(source_metrics, "Cylinder"),
    "only supported for 'to = \"Sphere\"'"
  )
  expect_equal(acousticTS:::.shape_equivalent_volume(c(0), c(1)), 0)
  expect_equal(
    acousticTS:::.canonical_radius_profile(c(0, 1), "Sphere", -1, 1),
    c(0, 0)
  )

  expect_error(
    acousticTS:::.shape_equivalent_radius_profile(1),
    "'position_matrix' must be a matrix."
  )
  expect_equal(
    acousticTS:::.shape_equivalent_radius_profile(cbind(w = c(2, 4))),
    c(1, 2)
  )

  local({
    testthat::local_mocked_bindings(
      .shape_radius_profile = function(position_matrix) rep(9, nrow(position_matrix)),
      .package = "acousticTS"
    )
    expect_equal(
      acousticTS:::.shape_equivalent_radius_profile(cbind(foo = c(1, 2))),
      c(9, 9)
    )
  })

  expect_equal(
    acousticTS:::.canonicalize_shape_fit_length_volume(
      list(length = 0.04, volume = 1e-5),
      "Sphere"
    )$length,
    2 * acousticTS:::.canonicalize_shape_fit_length_volume(
      list(length = 0.04, volume = 1e-5),
      "Sphere"
    )$radius,
    tolerance = 1e-12
  )
  expect_error(
    acousticTS:::.canonicalize_shape_fit_length_volume(
      list(length = 0.04, volume = 5e-5),
      "ProlateSpheroid"
    ),
    "cannot preserve both length and volume while remaining a prolate spheroid"
  )
  oblate_fit <- acousticTS:::.canonicalize_shape_fit_length_volume(
    list(length = 0.04, volume = 5e-5),
    "OblateSpheroid"
  )
  expect_true(oblate_fit$radius > oblate_fit$length / 2)
  expect_error(
    acousticTS:::.canonicalize_shape_fit_length_volume(
      list(length = 0.04, volume = 1e-5),
      "Weird"
    ),
    "Unsupported canonical target."
  )

  sphere_fit <- acousticTS:::.canonicalize_shape_fit_profile_l2(source_metrics, "Sphere")
  expect_true(is.finite(sphere_fit$radius))
  expect_true(is.finite(sphere_fit$objective))

  prolate_fit <- acousticTS:::.canonicalize_shape_fit_profile_l2(
    list(
      x = c(0, 0.5, 1),
      radius_eq = c(0.6, 0.6, 0.6),
      length = 1,
      volume = 1,
      max_radius = 0.6
    ),
    "ProlateSpheroid"
  )
  expect_true(is.finite(prolate_fit$length))
  expect_true(is.finite(prolate_fit$radius))
})

test_that("geometry profile helpers cover fallbacks, aliases, and geometry reconstruction", {
  row_major_body <- list(
    rpos = rbind(
      x = c(0, 1, 2),
      y = c(0, 0, 0),
      z = c(0, 0, 0)
    ),
    radius = 0.2
  )
  width_matrix <- cbind(x = c(0, 1, 2), w = c(2, 4, 6))
  widthless_matrix <- cbind(x = c(0, 1, 2), q = c(0, 0, 0))
  radius_body <- list(
    rpos = rbind(
      x = c(0, 1, 2),
      z = c(0, 0, 0)
    ),
    radius = c(1, 2, 3)
  )

  expect_error(acousticTS:::.shape_x(1), "'position_matrix' must be a matrix.")
  expect_error(acousticTS:::.canonicalize_x_order("bad"), "'x' must be numeric.")
  expect_equal(
    acousticTS:::.canonicalize_position_matrix(
      cbind(x = c(2, 0, 1), radius = c(3, 1, 2))
    )[, "x"],
    c(0, 1, 2)
  )
  expect_equal(acousticTS:::.shape_length(body = row_major_body), 2)
  expect_equal(
    acousticTS:::.shape_radius_profile(body = row_major_body, row_major = TRUE),
    rep(0.2, 3)
  )
  expect_equal(
    acousticTS:::.shape_radius_profile(position_matrix = width_matrix),
    c(1, 2, 3)
  )
  expect_error(
    acousticTS:::.shape_radius_profile(position_matrix = widthless_matrix),
    "Unable to resolve the nodewise radius profile"
  )

  expect_equal(
    acousticTS:::.shape_width_profile(
      body = radius_body,
      row_major = TRUE
    ),
    c(2, 4, 6)
  )
  expect_equal(
    acousticTS:::.shape_height_profile(
      body = radius_body,
      row_major = TRUE
    ),
    c(2, 4, 6)
  )

  expect_error(
    acousticTS:::brake_df(list(), 1),
    "Body shape information must be a list with a matrix element 'rpos'"
  )
  expect_error(
    acousticTS:::brake_df(list(rpos = matrix(1:4, nrow = 2), radius = c(0.1, 0.1)), -1),
    "Radius of curvature must be a positive-only, real number."
  )
  expect_error(
    acousticTS:::brake_df(list(rpos = matrix(1:4, nrow = 2), radius = c(0.1, 0.1)), 1, mode = "bad"),
    "Radius-of-curvature 'mode' must be either 'ratio' or 'measurement'."
  )
  expect_error(
    acousticTS:::brake_df(
      list(rpos = rbind(x = 0, y = 0, z = 0), radius = 0.1),
      1
    ),
    "Position matrix 'rpos' must have at least two columns."
  )

  bent <- acousticTS:::brake_df(
    list(
      rpos = rbind(
        x = c(0, 0.5, 1),
        y = c(0, 0, 0),
        z = c(0, 0, 0)
      ),
      radius = c(0.05, 0.05, 0.05)
    ),
    radius_curvature = 2,
    mode = "measurement"
  )
  expect_true(is.finite(bent$radius_curvature_ratio))

  sphere_obj <- fls_generate(
    shape = sphere(radius_body = 0.01, n_segments = 29),
    density_body = 1028.9,
    sound_speed_body = 1480.3
  )
  arbitrary_obj <- fls_generate(
    shape = arbitrary(
      x_body = c(0, 0.01, 0.02),
      radius_body = c(0, 0.005, 0)
    ),
    density_body = 1028.9,
    sound_speed_body = 1480.3
  )

  dwba_profile <- acousticTS:::.as_dwba_profile(sphere_obj)
  dwba_arbitrary <- acousticTS:::.as_dwba_profile(arbitrary_obj)
  expect_true(all(c("x", "y", "z", "zU", "zL") %in% rownames(dwba_profile@body$rpos)))
  expect_equal(dwba_arbitrary@body$radius, c(0, 0.005, 0))
  expect_equal(
    acousticTS:::.dwba_profile_object(sphere_obj)@body$radius,
    dwba_profile@body$radius
  )

  krm_arbitrary <- acousticTS:::.as_krm_profile(arbitrary_obj)
  expect_true(is.unsorted(krm_arbitrary@body$rpos["x", ], strictly = TRUE))
  expect_equal(
    acousticTS:::.krm_profile_object(arbitrary_obj)@body$radius,
    krm_arbitrary@body$radius
  )
})
