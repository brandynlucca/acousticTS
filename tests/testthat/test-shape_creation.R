library(acousticTS)

test_that("Shape creation functions work correctly", {
  # Test sphere function
  radius <- 0.05 # 5 cm radius
  n_segments <- 10
  sphere_obj <- sphere(radius_body = radius, n_segments = n_segments)

  # Check that it returns a Sphere object
  expect_s4_class(sphere_obj, "Sphere")

  # Check shape parameters
  expect_equal(sphere_obj@shape_parameters$radius, radius)
  expect_equal(sphere_obj@shape_parameters$n_segments, n_segments)
  expect_equal(sphere_obj@shape_parameters$diameter_shape, radius * 2)
  expect_equal(sphere_obj@shape_parameters$diameter_units, "m")

  # Check position matrix structure
  expect_true(is.matrix(sphere_obj@position_matrix))
  expect_equal(ncol(sphere_obj@position_matrix), 5) # x, y, z, zU, zL
  expect_true("x" %in% colnames(sphere_obj@position_matrix))
  expect_true("y" %in% colnames(sphere_obj@position_matrix))
  expect_true("z" %in% colnames(sphere_obj@position_matrix))
  expect_true("zU" %in% colnames(sphere_obj@position_matrix))
  expect_true("zL" %in% colnames(sphere_obj@position_matrix))

  # Check that zU and zL are symmetric around zero
  expect_equal(
    sphere_obj@position_matrix[, "zU"],
    -sphere_obj@position_matrix[, "zL"]
  )
})

test_that("Prolate spheroid creation works", {

  # Test prolate spheroid
  length_body <- 0.1 # 10 cm
  radius_body <- 0.02 # 2 cm
  n_segments <- 15

  prolate_obj <- prolate_spheroid(
    length_body = length_body,
    radius_body = radius_body,
    n_segments = n_segments
  )

  # Check that it returns a ProlateSpheroid object
  expect_s4_class(prolate_obj, "ProlateSpheroid")
  expect_s4_class(prolate_obj, "Shape") # ProlateSpheroid inherits from Shape

  # Check shape parameters
  expect_equal(prolate_obj@shape_parameters$length, length_body)
  expect_equal(prolate_obj@shape_parameters$radius,
    c(
      0.00000000, 0.00997775, 0.01359739, 0.01600000, 0.01768867,
      0.01885618, 0.01959592, 0.01995551, 0.01995551, 0.01959592,
      0.01885618, 0.01768867, 0.01600000, 0.01359739, 0.00997775,
      0.00000000
    ),
    tolerance = 1e-6
  )
  expect_equal(prolate_obj@shape_parameters$n_segments, n_segments)

  # Check position matrix structure
  expect_true(is.matrix(prolate_obj@position_matrix))
  expect_equal(ncol(prolate_obj@position_matrix), 5)

  # Test case where `length_radius_ratio` is used in lieu of radius_body
  prolate_obj <- prolate_spheroid(
    length_body = length_body,
    length_radius_ratio = 20,
    n_segments = n_segments
  )

  # Check that it returns a ProlateSpheroid object
  expect_s4_class(prolate_obj, "ProlateSpheroid")
  expect_s4_class(prolate_obj, "Shape") # ProlateSpheroid inherits from Shape
  expect_equal(prolate_obj@shape_parameters$length, length_body)
  expect_equal(prolate_obj@shape_parameters$radius,
    c(
      0.000000000, 0.002494438, 0.003399346, 0.004000000,
      0.004422166, 0.004714045, 0.004898979, 0.004988877,
      0.004988877, 0.004898979, 0.004714045, 0.004422166,
      0.004000000, 0.003399346, 0.002494438, 0.000000000
    ),
    tolerance = 1e-6
  )
  expect_equal(prolate_obj@shape_parameters$n_segments, n_segments)


  # Test error for under-specified shape info
  expect_error(
    prolate_spheroid(length_body = length_body),
    "Either 'radius' or 'length_radius_ratio' must be provided."
  )
})

test_that("Cylinder creation works", {

  # Test cylinder
  length_body <- 0.08
  radius_body <- 0.015
  n_segments <- 12

  cylinder_obj <- cylinder(
    length_body = length_body,
    radius_body = radius_body,
    n_segments = n_segments
  )

  # Check that it returns a Cylinder object
  expect_s4_class(cylinder_obj, "Cylinder")

  # Check shape parameters
  expect_equal(cylinder_obj@shape_parameters$length, length_body)
  expect_equal(
    cylinder_obj@shape_parameters$radius,
    rep(radius_body, n_segments + 1)
  )
  expect_equal(cylinder_obj@shape_parameters$n_segments, n_segments)

  # Check position matrix structure
  expect_true(is.matrix(cylinder_obj@position_matrix))
  expect_equal(ncol(cylinder_obj@position_matrix), 5)

  # Under-specified function
  expect_error(
    cylinder(length_body = 1),
    paste0(
      "Either 'radius' or 'length_radius_ratio' must be provided."
    )
  )
})

test_that("Arbitrary shape creation works", {
  # Test arbitrary shape
  x_body <- c(0, 0.01, 0.02, 0.03, 0.04)
  y_body <- rep(0, 5)
  z_body <- rep(0, 5)
  radius_body <- c(0.005, 0.008, 0.01, 0.008, 0.005)

  arbitrary_obj <- arbitrary(
    x_body = x_body,
    y_body = y_body,
    z_body = z_body,
    radius_body = radius_body
  )

  # Check that it returns an Arbitrary object
  expect_s4_class(arbitrary_obj, "Arbitrary")

  # Check shape parameters
  expect_equal(arbitrary_obj@shape_parameters$n_segments, length(x_body) - 1)
  expect_equal(arbitrary_obj@shape_parameters$diameter_units, "m")

  # Check position matrix structure
  expect_true(is.matrix(arbitrary_obj@position_matrix))
  expect_equal(ncol(arbitrary_obj@position_matrix), 6)
  expect_equal(nrow(arbitrary_obj@position_matrix), length(x_body))
})

test_that("Canonical KRM profiles match the common benchmark geometries after reversal", {
  krm_common_shape <- function(shape, n_points = 30L) {
    ang <- seq(-pi / 2, pi / 2, length.out = n_points)

    switch(
      shape,
      sphere = data.frame(
        x = 0.01 * (sin(ang) + 1),
        w = 0.02 * cos(ang),
        zU = 0.01 * cos(ang),
        zL = -0.01 * cos(ang)
      ),
      prolate = data.frame(
        x = 0.07 * (sin(ang) + 1),
        w = 0.02 * cos(ang),
        zU = 0.01 * cos(ang),
        zL = -0.01 * cos(ang)
      ),
      cylinder = data.frame(
        x = seq(0, 0.07, length.out = n_points),
        w = rep(0.02, n_points),
        zU = rep(0.01, n_points),
        zL = rep(-0.01, n_points)
      )
    )
  }

  for (shape_name in c("sphere", "prolate", "cylinder")) {
    obj <- switch(
      shape_name,
      sphere = fls_generate(
        sphere(radius_body = 0.01, n_segments = 29),
        density_body = 1028.9,
        sound_speed_body = 1480.3
      ),
      prolate = fls_generate(
        prolate_spheroid(
          length_body = 0.14,
          radius_body = 0.01,
          n_segments = 29
        ),
        density_body = 1028.9,
        sound_speed_body = 1480.3
      ),
      cylinder = fls_generate(
        cylinder(
          length_body = 0.07,
          radius_body = 0.01,
          n_segments = 29
        ),
        density_body = 1028.9,
        sound_speed_body = 1480.3
      )
    )
    profiled <- acousticTS:::.krm_profile_object(obj)
    body <- profiled@body$rpos
    acoustic <- data.frame(x = rev(body["x", ]),
                           w = rev(body["w", ]),
                           zU = rev(body["zU", ]),
                           zL = rev(body["zL", ]))

    expect_equal(acoustic, krm_common_shape(shape_name), tolerance = 1e-12)
  }
})

test_that("arbitrary FLS radius metadata uses semi-diameter for zU/zL inputs", {
  obj <- fls_generate(
    arbitrary(
      x_body = c(0, 0.01, 0.02),
      w_body = c(0, 0.01, 0),
      zU_body = c(0, 0.005, 0),
      zL_body = c(0, -0.005, 0)
    ),
    density_body = 1028.9,
    sound_speed_body = 1480.3
  )

  expect_equal(obj@shape_parameters$radius, 0.005)
})

test_that("gas_generate preserves absolute fluid properties when supplied", {
  obj <- gas_generate(
    shape = "sphere",
    radius_body = 0.01,
    density_fluid = 1.24,
    sound_speed_fluid = 345
  )

  expect_null(obj@body$g)
  expect_null(obj@body$h)
  expect_equal(obj@body$density, 1.24)
  expect_equal(obj@body$sound_speed, 345)
})

test_that("gas_generate preserves canonical prolate metadata for modal-series models", {
  obj <- gas_generate(
    shape = "prolate_spheroid",
    length_body = 0.14,
    radius_body = 0.01,
    density_fluid = 1.24,
    sound_speed_fluid = 345
  )

  expect_equal(obj@shape_parameters$length, 0.14)
  expect_equal(obj@shape_parameters$radius, 0.01)
  expect_equal(obj@shape_parameters$semimajor_length, 0.07)
  expect_equal(obj@shape_parameters$semiminor_length, 0.01)
  expect_length(obj@body$radius, 101)
})

test_that("'polyonmial_cylinder' function works", {
  # Test polynomial-based cylinder shape generation
  poly_vec <- c(0.83, 0.36, -2.10, -1.20, 0.63, 0.82, 0.64)
  length_body <- 1.0
  radius_body <- 0.2
  pcyl_obj <- polynomial_cylinder(
    length_body = length_body,
    radius_body = radius_body,
    polynomial = poly_vec
  )

  # Check that it returns a PolynomialCylinder object
  expect_s4_class(pcyl_obj, "PolynomialCylinder")
  expect_s4_class(pcyl_obj, "Shape")
})

test_that("create_shape wrapper function works", {
  # Test create_shape with sphere
  sphere_via_wrapper <- create_shape("sphere",
    radius_body = 0.03,
    n_segments = 8
  )
  expect_s4_class(sphere_via_wrapper, "Sphere")
  expect_equal(sphere_via_wrapper@shape_parameters$radius, 0.03)

  # Test create_shape with cylinder
  cylinder_via_wrapper <- create_shape("cylinder",
    length_body = 0.06,
    radius_body = 0.01,
    n_segments = 10
  )
  expect_s4_class(cylinder_via_wrapper, "Cylinder")
  expect_equal(cylinder_via_wrapper@shape_parameters$length, 0.06)

  # Test error for invalid shape type
  expect_error(create_shape("invalid_shape", radius_body = 0.01),
    regexp = "'invalid_shape' of mode 'function' was not found"
  )
})
