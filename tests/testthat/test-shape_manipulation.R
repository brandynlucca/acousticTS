library(acousticTS)

test_that("translate_shape() and reanchor_shape() move Shape geometry predictably", {
  shape_obj <- cylinder(length_body = 0.05, radius_body = 0.003, n_segments = 10)
  moved_shape <- translate_shape(shape_obj, x_offset = 0.01, z_offset = 0.002)

  expect_s4_class(moved_shape, "Cylinder")
  expect_equal(
    range(extract(moved_shape, c("position_matrix", "x"))),
    range(extract(shape_obj, c("position_matrix", "x"))) + 0.01
  )
  expect_equal(
    extract(moved_shape, c("position_matrix", "zU")),
    extract(shape_obj, c("position_matrix", "zU")) + 0.002
  )

  centered_shape <- reanchor_shape(shape_obj, anchor = "center", at = 0)
  expect_equal(
    mean(range(extract(centered_shape, c("position_matrix", "x")))),
    0,
    tolerance = 1e-12
  )
})

test_that("flip_shape() preserves the x grid and reverses or mirrors profiles", {
  shape_obj <- arbitrary(
    x_body = c(0, 0.01, 0.02, 0.03),
    radius_body = c(0, 0.004, 0.002, 0)
  )

  flipped_x <- flip_shape(shape_obj, axis = "x")
  expect_equal(
    extract(flipped_x, c("position_matrix", "x")),
    extract(shape_obj, c("position_matrix", "x"))
  )
  expect_equal(
    extract(flipped_x, c("position_matrix", "zU")),
    rev(extract(shape_obj, c("position_matrix", "zU")))
  )

  flipped_z <- flip_shape(shape_obj, axis = "z")
  expect_equal(
    extract(flipped_z, c("position_matrix", "zU")),
    -extract(shape_obj, c("position_matrix", "zL"))
  )
  expect_equal(
    extract(flipped_z, c("position_matrix", "zL")),
    -extract(shape_obj, c("position_matrix", "zU"))
  )
})

test_that("resample_shape() updates shape and scatterer segment counts", {
  shape_obj <- sphere(radius_body = 0.01, n_segments = 12)
  shape_fine <- resample_shape(shape_obj, n_segments = 40)
  expect_equal(extract(shape_fine, c("shape_parameters", "n_segments")), 40)

  obj <- fls_generate(
    shape = cylinder(length_body = 0.05, radius_body = 0.003, n_segments = 12),
    density_body = 1045,
    sound_speed_body = 1520
  )
  obj_fine <- resample_shape(obj, n_segments = 60)
  expect_equal(extract(obj_fine, c("shape_parameters", "n_segments")), 60)
  expect_equal(ncol(extract(obj_fine, "body")$rpos), 61)
})

test_that("inflate_shape() and smooth_shape() relabel profile-edited shapes as Arbitrary", {
  shape_obj <- arbitrary(
    x_body = c(0, 0.01, 0.02, 0.03, 0.04),
    radius_body = c(0, 0.002, 0.004, 0.002, 0)
  )
  pinched_shape <- inflate_shape(
    shape_obj,
    x_range = c(0.01, 0.03),
    scale = 0.5
  )

  expect_s4_class(pinched_shape, "Arbitrary")
  expect_lt(
    max(extract(pinched_shape, c("shape_parameters", "radius")), na.rm = TRUE),
    max(extract(shape_obj, c("shape_parameters", "radius")), na.rm = TRUE)
  )

  smoothed_shape <- smooth_shape(shape_obj, span = 3)
  expect_s4_class(smoothed_shape, "Arbitrary")
  expect_equal(
    extract(smoothed_shape, c("shape_parameters", "n_segments")),
    extract(shape_obj, c("shape_parameters", "n_segments"))
  )
})

test_that("offset_component() repositions internal components and enforces containment", {
  fish <- sbf_generate(
    x_body = c(0, 0.1),
    w_body = c(0.006, 0.008),
    zU_body = c(0.001, 0.002),
    zL_body = c(-0.001, -0.002),
    x_bladder = c(0.02, 0.08),
    w_bladder = c(0, 0),
    zU_bladder = c(0.0012, 0.0012),
    zL_bladder = c(-0.0012, -0.0012),
    density_body = 1040,
    density_bladder = 1.2,
    sound_speed_body = 1500,
    sound_speed_bladder = 340
  )

  shifted_fish <- offset_component(fish, component = "bladder", x_offset = 0.003)
  expect_equal(
    min(extract(shifted_fish, "bladder")$rpos["x_bladder", ]),
    min(extract(fish, "bladder")$rpos["x_bladder", ]) + 0.003
  )

  expect_error(
    offset_component(
      fish,
      component = "bladder",
      z_offset = 0.02,
      containment = "error"
    ),
    "Swimbladder exceeds body bounds"
  )
})

test_that("translate_shape() warns when y_offset is unsupported for row-major profiles", {
  obj <- fls_generate(
    shape = cylinder(length_body = 0.05, radius_body = 0.003, n_segments = 12),
    density_body = 1045,
    sound_speed_body = 1520
  )

  expect_warning(
    translate_shape(obj, y_offset = 0.01),
    "y_offset"
  )
})
