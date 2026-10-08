library(acousticTS)

test_that("profile operations agree across matrix orientations", {
  x <- seq(0, 0.1, length.out = 9)
  z <- 0.002 * sin(pi * x / 0.1)
  radius <- c(1, 2, 2, 3, 2, 3, 2, 2, 1) * 0.001
  columns <- cbind(x = x, w = radius, z = z, zU = z + radius, zL = z - radius)
  rows <- t(columns)
  expect_equal(acousticTS:::.shape_centerline_z(body = list(rpos = rows)), z)
  expect_equal(acousticTS:::.shape_centerline_z(columns), z)
  reversed <- rows[, rev(seq_len(ncol(rows)))]
  expect_equal(
    acousticTS:::.canonicalize_position_matrix(reversed, row_major = TRUE),
    rows
  )
  flip <- acousticTS:::.shape_flip_matrix
  expect_equal(
    flip(columns, "shape_column_major", "z"),
    t(flip(rows, "profile_row_major", "z"))
  )
  expect_equal(flip(columns, "shape_column_major", "z")[, "z"], -z)
  scale <- acousticTS:::.shape_local_scale_matrix
  expect_equal(
    scale(columns, "shape_column_major", scale = 1.2, axis = "height"),
    t(scale(rows, "profile_row_major", scale = 1.2, axis = "height"))
  )
  smooth <- acousticTS:::.shape_smooth_matrix
  expect_equal(
    smooth(columns, "shape_column_major", span = 3),
    t(smooth(rows, "profile_row_major", span = 3))
  )
  expect_warning(
    shifted <- acousticTS:::.translate_column_position_matrix(
      columns, x_offset = 1, y_offset = 2, z_offset = 3
    ),
    "does not store an explicit lateral centerline"
  )
  expect_equal(shifted[, "x"], x + 1)
  expect_equal(shifted[, "z"], z + 3)
  expect_equal(shifted[, "w"], radius)
  zero_fields <- list(x = x)
  expect_equal(acousticTS:::.reforge_profile_centerline(zero_fields), rep(0, 9))
  expect_equal(
    acousticTS:::.reforge_profile_half_height(zero_fields), rep(0, 9)
  )
})

test_that("backbone edits synchronize geometry and enforce containment", {
  fish <- bbf_generate(
    body_shape = cylinder(0.08, 0.006),
    backbone_shape = cylinder(
      0.04,
      8e-04
    ),
    g_body = 1.04,
    h_body = 1.02,
    density_backbone = 1900,
    sound_speed_longitudinal_backbone = 3500,
    sound_speed_transversal_backbone = 1700,
    x_offset_backbone = 0.02
  )
  moved <- offset_component(fish, "backbone", x_offset = 0.005)
  expect_equal(
    moved@backbone$rpos["x", ],
    fish@backbone$rpos["x", ] +
      0.005
  )
  expect_equal(moved@components$backbone, moved@backbone)
  expect_equal(moved@body, fish@body)
  expect_error(
    offset_component(fish, "backbone", z_offset = 0.01, containment = "error"),
    "Backbone exceeds body bounds"
  )

  shell <- fixture_sphere("shelled_liquid")
  inflated <- inflate_shape(shell, scale = 1.2, profile = "box")
  expect_equal(inflated@shell$radius, 1.2 * shell@shell$radius)
  expect_equal(inflated@fluid, shell@fluid)
})

test_that("profile metadata updates retain dimensions in metres", {
  summary <- list(
    length = 0.08,
    n_segments = 4L,
    radius_profile = c(
      0,
      0.01,
      0.02,
      0.01,
      0
    ),
    max_radius = 0.02
  )
  params <- list(
    shape = "ProlateSpheroid",
    length = 0.04,
    n_segments = 2L,
    radius_shape = numeric(3),
    diameter_shape = numeric(3),
    radius = 0.01,
    diameter = 0.02,
    mean_radius = 0,
    max_radius = 0,
    semimajor_length = 0.02,
    semiminor_length = 0.01,
    length_radius_ratio = 4
  )
  actual <- acousticTS:::.update_manipulated_shape_params(
    params,
    summary,
    force_arbitrary = TRUE
  )
  expect_equal(actual$shape, "Arbitrary")
  expect_equal(actual$radius_shape, summary$radius_profile)
  expect_equal(actual$diameter_shape, 2 * summary$radius_profile)
  expect_equal(actual$diameter, 0.04)
  expect_equal(actual$semimajor_length, 0.04)
  expect_equal(actual$semiminor_length, 0.02)
  expect_equal(actual$length_radius_ratio, 4)
  expect_equal(actual$mean_radius, mean(summary$radius_profile))

  profile <- cbind(x = 0:4, radius = c(0, 1, 2, 1, 0))
  scaled <- acousticTS:::.shape_local_scale_matrix(
    profile,
    "shape_column_major",
    scale = 2,
    profile = "box"
  )
  expect_equal(scaled[, "radius"], 2 * profile[, "radius"])
  smoothed <- acousticTS:::.shape_smooth_matrix(
    profile,
    "shape_column_major",
    span = 3
  )
  expect_equal(smoothed[, "radius"], c(0, 1, 4 / 3, 1, 0))
  expect_equal(smoothed[, "x"], profile[, "x"])
  bad <- fls_generate(
    shape = cylinder(0.04, 0.003),
    g_body = 1.03,
    h_body = 1.02
  )
  bad@body$rpos <- NULL
  expect_error(translate_shape(bad), "valid position matrix")
  expect_error(
    acousticTS:::.geometry_storage(matrix(1, 1, 1)),
    "Unable to determine"
  )
})

test_that("profile edits preserve coordinates and dimensions", {
  shape <- arbitrary(
    x_body = seq(0, 0.04, length.out = 5),
    w_body = c(
      0,
      0.002,
      0.006,
      0.004,
      0
    ),
    zU_body = c(0, 0.003, 0.005, 0.002, 0),
    zL_body = c(0, -0.001, -0.003, -0.002, 0)
  )
  object <- fls_generate(shape = shape, g_body = 1.03, h_body = 1.02)
  original <- object@body$rpos
  for (axis in c("x", "z")) {
    flipped <- flip_shape(object, axis = axis)
    expect_equal(flip_shape(flipped, axis = axis)@body$rpos, original)
    expect_equal(flipped@body$rpos[1, ], original[1, ])
  }
  for (axis in c("radius", "width", "height")) {
    inflated <- inflate_shape(object, scale = 2, axis = axis, profile = "box")
    expected <- inflate_shape(shape, scale = 2, axis = axis, profile = "box")
    actual_fields <- acousticTS:::.profile_fields(inflated@body$rpos)
    expected_fields <- acousticTS:::.profile_fields(
      fls_generate(
        shape = expected,
        g_body = 1.03,
        h_body = 1.02
      )@body$rpos
    )
    expect_equal(actual_fields, expected_fields, tolerance = 1e-12)
  }
  smooth <- smooth_shape(object, span = 4)
  smooth_odd <- smooth_shape(object, span = 5)
  expect_equal(smooth@body$rpos, smooth_odd@body$rpos)
  expect_equal(smooth@body$rpos[, c(1, 5)], original[, c(1, 5)])
  expect_equal(smooth@body$rpos[1, ], original[1, ])
  expect_error(translate_shape(1), "Shape or Scatterer")
  expect_error(
    translate_shape(object, component = "bladder"),
    "does not contain"
  )
  expect_error(offset_component(shape), "Scatterer object")
  expect_error(resample_shape(object, n_segments = 0), "positive integer")
  expect_error(inflate_shape(object, scale = 0), "positive number")
  expect_error(inflate_shape(object, x_range = 1), "length two")
  expect_error(inflate_shape(object, x_range = c(0, Inf)), "finite axial")
  expect_error(smooth_shape(object, span = 2), "integer >= 3")
  expect_equal(
    acousticTS:::.shape_window(0:4, c(2, 2)),
    c(
      0,
      0,
      1,
      0,
      0
    )
  )
})

test_that(
  "translate_shape() and reanchor_shape() move Shape geometry predictably",
  {
    shape_obj <- cylinder(
      length_body = 0.05, radius_body = 0.003, n_segments =
        10
    )
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
  }
)

test_that(
  "flip_shape() preserves the x grid and reverses or mirrors profiles",
  {
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
  }
)

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

test_that(
  paste0(
    "inflate_shape() and smooth_shape() relabel profile-edited shapes as ",
    "Arbitrary"
  ),
  {
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
      max(
        extract(pinched_shape, c("shape_parameters", "radius")),
        na.rm = TRUE
      ),
      max(extract(shape_obj, c("shape_parameters", "radius")), na.rm = TRUE)
    )

    smoothed_shape <- smooth_shape(shape_obj, span = 3)
    expect_s4_class(smoothed_shape, "Arbitrary")
    expect_equal(
      extract(smoothed_shape, c("shape_parameters", "n_segments")),
      extract(shape_obj, c("shape_parameters", "n_segments"))
    )
  }
)

test_that(
  paste0(
    "offset_component() repositions internal components and enforces ",
    "containment"
  ),
  {
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

    shifted_fish <- offset_component(fish,
      component = "bladder", x_offset =
        0.003
    )
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
  }
)

test_that(
  paste0(
    "translate_shape() warns when y_offset is unsupported for row-major ",
    "profiles"
  ),
  {
    obj <- fls_generate(
      shape = cylinder(
        length_body = 0.05,
        radius_body = 0.003,
        n_segments = 12
      ),
      density_body = 1045,
      sound_speed_body = 1520
    )

    expect_warning(
      translate_shape(obj, y_offset = 0.01),
      "y_offset"
    )
  }
)
