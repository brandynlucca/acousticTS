test_that("Shape creation functions work correctly", {
  library(acousticTS)
  
  # Test sphere function
  radius <- 0.05  # 5 cm radius
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
  expect_equal(ncol(sphere_obj@position_matrix), 5)  # x, y, z, zU, zL
  expect_true("x" %in% colnames(sphere_obj@position_matrix))
  expect_true("y" %in% colnames(sphere_obj@position_matrix))
  expect_true("z" %in% colnames(sphere_obj@position_matrix))
  expect_true("zU" %in% colnames(sphere_obj@position_matrix))
  expect_true("zL" %in% colnames(sphere_obj@position_matrix))
  
  # Check that zU and zL are symmetric around zero
  expect_equal(sphere_obj@position_matrix[, "zU"], -sphere_obj@position_matrix[, "zL"])
})

test_that("Prolate spheroid creation works", {
  library(acousticTS)
  
  # Test prolate spheroid
  length_body <- 0.1  # 10 cm
  radius_body <- 0.02  # 2 cm
  n_segments <- 15
  
  prolate_obj <- prolate_spheroid(length_body = length_body, 
                                  radius_body = radius_body,
                                  n_segments = n_segments)
  
  # Check that it returns a ProlateSpheroid object
  expect_s4_class(prolate_obj, "ProlateSpheroid")
  expect_s4_class(prolate_obj, "Shape") # ProlateSpheroid inherits from Shape
  
  # Check shape parameters
  expect_equal(prolate_obj@shape_parameters$length, length_body)
  expect_equal(prolate_obj@shape_parameters$radius, 
               c(0.00000000, 0.00718022, 0.01200000, 0.01490712, 0.01691810, 
                 0.01833030, 0.01927578, 0.01982142, 0.02000000, 0.01982142, 
                 0.01927578, 0.01833030, 0.01691810, 0.01490712, 0.01200000, 
                 0.00718022),
               tolerance=1e-6)
  expect_equal(prolate_obj@shape_parameters$n_segments, n_segments)
  
  # Check position matrix structure
  expect_true(is.matrix(prolate_obj@position_matrix))
  expect_equal(ncol(prolate_obj@position_matrix), 5)
})

test_that("Cylinder creation works", {
  library(acousticTS)
  
  # Test cylinder
  length_body <- 0.08
  radius_body <- 0.015
  n_segments <- 12
  
  cylinder_obj <- cylinder(length_body = length_body,
                          radius_body = radius_body,
                          n_segments = n_segments)
  
  # Check that it returns a Cylinder object
  expect_s4_class(cylinder_obj, "Cylinder")
  
  # Check shape parameters
  expect_equal(cylinder_obj@shape_parameters$length, length_body)
  expect_equal(cylinder_obj@shape_parameters$radius, 
               c(0.0, rep(radius_body, n_segments)))
  expect_equal(cylinder_obj@shape_parameters$n_segments, n_segments)
  
  # Check position matrix structure
  expect_true(is.matrix(cylinder_obj@position_matrix))
  expect_equal(ncol(cylinder_obj@position_matrix), 5)
})

test_that("Arbitrary shape creation works", {
  library(acousticTS)
  
  # Test arbitrary shape
  x_body <- c(0, 0.01, 0.02, 0.03, 0.04)
  y_body <- rep(0, 5)
  z_body <- rep(0, 5)
  radius_body <- c(0.005, 0.008, 0.01, 0.008, 0.005)
  
  arbitrary_obj <- arbitrary(x_body = x_body,
                            y_body = y_body,
                            z_body = z_body,
                            radius_body = radius_body)
  
  # Check that it returns an Arbitrary object
  expect_s4_class(arbitrary_obj, "Arbitrary")
  
  # Check shape parameters
  expect_equal(arbitrary_obj@shape_parameters$n_segments, length(x_body) - 1)
  expect_equal(arbitrary_obj@shape_parameters$diameter_units, "m")
  
  # Check position matrix structure
  expect_true(is.matrix(arbitrary_obj@position_matrix))
  expect_equal(ncol(arbitrary_obj@position_matrix), 5)
  expect_equal(nrow(arbitrary_obj@position_matrix), length(x_body))
})

test_that("create_shape wrapper function works", {
  library(acousticTS)
  
  # Test create_shape with sphere
  sphere_via_wrapper <- create_shape("sphere", radius_body = 0.03, 
                                     n_segments = 8)
  expect_s4_class(sphere_via_wrapper, "Sphere")
  expect_equal(sphere_via_wrapper@shape_parameters$radius, 0.03)
  
  # Test create_shape with cylinder  
  cylinder_via_wrapper <- create_shape("cylinder", 
                                      length_body = 0.06, 
                                      radius_body = 0.01,
                                      n_segments = 10)
  expect_s4_class(cylinder_via_wrapper, "Cylinder")
  expect_equal(cylinder_via_wrapper@shape_parameters$length, 0.06)
  
  # Test error for invalid shape type
  expect_error(create_shape("invalid_shape", radius_body = 0.01), 
               regexp = "'invalid_shape' of mode 'function' was not found")
})
