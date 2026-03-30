library(acousticTS)

test_that("shape resolution and validator helpers cover canonical constructor branches", {
  sphere_shape <- sphere(radius_body = 0.01, n_segments = 12)
  arbitrary_shape <- arbitrary(
    x_body = c(0, 0.01, 0.02),
    radius_body = c(0, 0.01, 0)
  )

  expect_identical(acousticTS:::.resolve_shape(sphere_shape, list()), sphere_shape)
  expect_s4_class(
    acousticTS:::.resolve_shape(
      "sphere",
      list(radius_body = 0.01, n_segments = 12, ignored = TRUE)
    ),
    "Sphere"
  )
  expect_s4_class(
    acousticTS:::.resolve_shape(
      "arbitrary",
      list(
        x_body = c(0, 0.01, 0.02),
        radius_body = c(0, 0.01, 0),
        ignored = TRUE
      )
    ),
    "Arbitrary"
  )

  expect_error(
    acousticTS:::.validate_shape_requirements("sphere", list(radius_body = NA_real_)),
    "Sphere requires 'radius_body'"
  )
  expect_error(
    acousticTS:::.validate_shape_requirements("cylinder", list(radius_body = 0.01)),
    "Cylinder requires 'length_body'"
  )
  expect_error(
    acousticTS:::.validate_shape_requirements("cylinder", list(length_body = 0.02)),
    "Cylinder requires either 'radius_body' or 'length_radius_ratio'"
  )
  expect_error(
    acousticTS:::.validate_shape_requirements("prolate_spheroid", list(radius_body = 0.01)),
    "Prolate spheroid requires 'length_body' or 'semimajor_length'"
  )
  expect_error(
    acousticTS:::.validate_shape_requirements("prolate_spheroid", list(length_body = 0.02)),
    "Prolate spheroid requires 'radius_body', 'semiminor_length', or 'length_radius_ratio'"
  )
  expect_error(
    acousticTS:::.validate_shape_requirements("oblate_spheroid", list(radius_body = 0.01)),
    "Oblate spheroid requires 'length_body' or 'semiminor_length'"
  )
  expect_error(
    acousticTS:::.validate_shape_requirements("oblate_spheroid", list(length_body = 0.02)),
    "Oblate spheroid requires 'radius_body', 'semimajor_length', or 'length_radius_ratio'"
  )
  expect_error(
    acousticTS:::.validate_shape_requirements(
      "polynomial_cylinder",
      list(radius_body = 0.01, polynomial = 1)
    ),
    "Polynomial cylinder requires 'length_body'"
  )
  expect_error(
    acousticTS:::.validate_shape_requirements(
      "polynomial_cylinder",
      list(length_body = 0.02, polynomial = 1)
    ),
    "Polynomial cylinder requires 'radius_body'"
  )
  expect_error(
    acousticTS:::.validate_shape_requirements(
      "polynomial_cylinder",
      list(length_body = 0.02, radius_body = 0.01)
    ),
    "Polynomial cylinder requires 'polynomial' coefficients"
  )
  expect_null(acousticTS:::.validate_shape_requirements("arbitrary", list()))

  expect_error(
    acousticTS:::.validate_material_pair(
      contrast = NULL,
      absolute = NULL,
      contrast_name = "g_body",
      absolute_name = "density_body",
      component_label = "body"
    ),
    "Supply either g_body or density_body for the body"
  )
  expect_invisible(
    acousticTS:::.validate_material_pair(
      contrast = 1,
      absolute = NULL,
      contrast_name = "g_body",
      absolute_name = "density_body",
      component_label = "body"
    )
  )

  expect_identical(
    acousticTS:::.resolve_scatterer_shape_input(sphere_shape, list()),
    sphere_shape
  )
  expect_s4_class(
    acousticTS:::.resolve_scatterer_shape_input(
      NULL,
      list(x_body = c(0, 0.01, 0.02), radius_body = c(0, 0.01, 0))
    ),
    "Arbitrary"
  )
  expect_error(
    acousticTS:::.resolve_scatterer_shape_input(NULL, list()),
    "Supply 'shape' as a pre-built Shape object"
  )
  expect_warning(
    legacy_shape <- acousticTS:::.resolve_scatterer_shape_input(
      "sphere",
      list(radius_body = 0.01, n_segments = 12)
    ),
    "Character-based 'shape' dispatch"
  )
  expect_s4_class(legacy_shape, "Sphere")
  expect_s4_class(arbitrary_shape, "Arbitrary")
})

test_that("geometry and component helpers cover fallbacks and validation branches", {
  shape_matrix <- cbind(
    x = c(0, 1, 2),
    w = c(1, 2, 3),
    zU = c(0.5, 0.6, 0.7),
    zL = c(-0.5, -0.6, -0.7)
  )
  row_major <- rbind(
    x = c(0, 1, 2),
    w = c(1, 2, 3),
    zU = c(0.5, 0.6, 0.7),
    zL = c(-0.5, -0.6, -0.7)
  )

  expect_error(
    acousticTS:::.geometry_axis_values(1, axis = "x"),
    "Geometry must be stored as a matrix"
  )
  expect_equal(
    acousticTS:::.geometry_axis_values(shape_matrix, axis = "w", context = "Shape"),
    c(1, 2, 3)
  )
  expect_equal(
    acousticTS:::.geometry_axis_values(row_major, axis = "zU", row_major = TRUE),
    c(0.5, 0.6, 0.7)
  )
  expect_equal(
    acousticTS:::.geometry_axis_values(
      matrix(c(1, 2, 3), ncol = 1),
      axis = "x"
    ),
    c(1, 2, 3)
  )
  expect_equal(
    acousticTS:::.geometry_axis_values(
      matrix(c(1, 2, 3), ncol = 1),
      axis = "w",
      default = c(4, 5, 6)
    ),
    c(4, 5, 6)
  )

  expect_error(
    acousticTS:::.validate_geometry_contract(1),
    "Geometry must be stored as a matrix"
  )
  expect_error(
    acousticTS:::.validate_geometry_contract(matrix(1:4, nrow = 2)),
    "missing column names required by the internal geometry contract"
  )
  expect_error(
    acousticTS:::.validate_geometry_contract(cbind(y = c(0, 1), w = c(1, 1))),
    "must define an x-axis column"
  )
  expect_error(
    acousticTS:::.validate_geometry_contract(cbind(x = c(0, 1), zU = c(1, 1))),
    "must define both upper and lower height coordinates"
  )

  expect_error(acousticTS:::.shape_segment_count(1), "'position_matrix' must be a matrix")
  expect_equal(acousticTS:::.shape_segment_count(shape_matrix), 2L)
  expect_equal(acousticTS:::.shape_segment_count(row_major, row_major = TRUE), 2L)

  expect_equal(acousticTS:::.scatterer_shape_name(sphere(radius_body = 0.01), "arbitrary"), "Arbitrary")
  expect_match(
    acousticTS:::.scatterer_shape_name(sphere(radius_body = 0.01)),
    "Sphere"
  )

  expect_equal(
    acousticTS:::.extract_shape_component_row(row_major, candidates = c("radius"), default = c(9, 9, 9)),
    c(9, 9, 9)
  )
  expect_equal(
    acousticTS:::.extract_shape_component_row(matrix(1:6, nrow = 2), candidates = "x"),
    c(1, 3, 5)
  )

  translated <- acousticTS:::.translate_shape_position_matrix(
    shape_matrix,
    x_offset = 1,
    y_offset = 2,
    z_offset = -1
  )
  expect_equal(translated[, "x"], c(1, 2, 3))
  expect_equal(translated[, "w"], c(3, 4, 5))
  expect_equal(translated[, "zU"], c(-0.5, -0.4, -0.3))
  expect_equal(translated[, "zL"], c(-1.5, -1.6, -1.7))

  expect_error(
    acousticTS:::.validate_component_rpos(matrix(1, nrow = 4, ncol = 1), "Body"),
    "Body shape must have at least two points"
  )
  expect_invisible(acousticTS:::.validate_component_rpos(matrix(1:8, nrow = 4, ncol = 2), "Body"))
})

test_that("scatterer component builders and ESS helpers cover remaining branches", {
  sphere_shape <- sphere(radius_body = 0.01, n_segments = 12)
  cylinder_shape <- cylinder(length_body = 0.04, radius_body = 0.01, n_segments = 12)
  prolate_shape <- prolate_spheroid(length_body = 0.08, radius_body = 0.01, n_segments = 12)
  oblate_shape <- oblate_spheroid(length_body = 0.04, radius_body = 0.03, n_segments = 12)
  arbitrary_shape <- arbitrary(
    x_body = c(0, 0.02, 0.04),
    w_body = c(0, 0.02, 0),
    zU_body = c(0, 0.01, 0),
    zL_body = c(0, -0.01, 0),
    radius_body = c(0, 0.01, 0)
  )

  fluid_component <- acousticTS:::.build_row_major_fluid_component(
    sphere_shape,
    theta = pi / 3,
    density = 1025,
    sound_speed = 1500
  )
  elastic_component <- acousticTS:::.build_row_major_elastic_component(
    sphere_shape,
    density = 1100,
    sound_speed_longitudinal = 2500,
    sound_speed_transversal = 1200
  )
  elastic_override <- acousticTS:::.build_row_major_elastic_component(
    sphere_shape,
    position_matrix_override = extract(sphere_shape, "position_matrix"),
    density = 1100,
    sound_speed_longitudinal = 2500,
    sound_speed_transversal = 1200
  )

  expect_true(is.matrix(fluid_component$rpos))
  expect_true(is.matrix(elastic_component$rpos))
  expect_true(is.matrix(elastic_override$rpos))
  expect_equal(fluid_component$theta, pi / 3)

  diameter_radius <- acousticTS:::.normalize_shape_diameter(list(radius = c(1, 2, 3)))
  diameter_both <- acousticTS:::.normalize_shape_diameter(list(
    radius = c(1, 2, 3),
    diameter = c(2, 4, 6)
  ))
  diameter_shape <- acousticTS:::.normalize_shape_diameter(list(
    diameter_shape = c(2, 4, 6)
  ))
  expect_equal(diameter_radius$diameter, 6)
  expect_equal(diameter_both$radius_shape, c(1, 2, 3))
  expect_equal(diameter_both$diameter_shape, c(2, 4, 6))
  expect_equal(diameter_shape$diameter, 6)

  expect_null(acousticTS:::.inner_shape_from_shell(sphere_shape, NULL))
  expect_error(
    acousticTS:::.inner_shape_from_shell(sphere_shape, -0.001),
    "'shell_thickness' must be a single positive number"
  )
  expect_error(
    acousticTS:::.inner_shape_from_shell(sphere_shape, 0.02),
    "'shell_thickness' must be smaller than the shell radius"
  )

  inner_sphere <- acousticTS:::.inner_shape_from_shell(sphere_shape, 0.002)
  inner_cylinder <- acousticTS:::.inner_shape_from_shell(cylinder_shape, 0.002)
  inner_prolate <- acousticTS:::.inner_shape_from_shell(prolate_shape, 0.002)
  inner_oblate <- acousticTS:::.inner_shape_from_shell(oblate_shape, 0.002)
  inner_arbitrary <- acousticTS:::.inner_shape_from_shell(arbitrary_shape, 0.002)
  expect_s4_class(inner_sphere, "Sphere")
  expect_s4_class(inner_cylinder, "Cylinder")
  expect_s4_class(inner_prolate, "ProlateSpheroid")
  expect_s4_class(inner_oblate, "OblateSpheroid")
  expect_s4_class(inner_arbitrary, "Arbitrary")
  expect_equal(extract(inner_sphere, "shape_parameters")$radius, 0.008)

  shape_object_components <- acousticTS:::.resolve_ess_shape_components(
    shape = sphere_shape,
    arguments = list(shape = sphere_shape),
    shell_thickness = 0.002
  )
  null_shape_components <- acousticTS:::.resolve_ess_shape_components(
    shape = NULL,
    arguments = list(),
    radius_shell = c(0, 0.01, 0),
    shell_thickness = 0.002,
    x_body = c(0, 0.01, 0.02),
    y_body = c(0, 0, 0),
    z_body = c(0, 0, 0)
  )
  arbitrary_components <- acousticTS:::.resolve_ess_shape_components(
    shape = "arbitrary",
    arguments = list(),
    radius_shell = c(0, 0.01, 0),
    shell_thickness = 0.002,
    x_body = c(0, 0.01, 0.02),
    y_body = c(0, 0, 0),
    z_body = c(0, 0, 0)
  )
  expect_warning(
    legacy_components <- acousticTS:::.resolve_ess_shape_components(
      shape = "sphere",
      arguments = list(radius_body = 0.01, n_segments = 12),
      radius_shell = 0.01,
      shell_thickness = 0.002
    ),
    "Character-based 'shape' dispatch"
  )

  expect_s4_class(shape_object_components$shell, "Sphere")
  expect_s4_class(null_shape_components$shell, "Arbitrary")
  expect_s4_class(arbitrary_components$shell, "Arbitrary")
  expect_s4_class(legacy_components$shell, "Sphere")
  expect_error(
    acousticTS:::.resolve_ess_shape_components(
      shape = NULL,
      arguments = list(),
      shell_thickness = 0.002
    ),
    "Supply 'shape' as a pre-built Shape object"
  )
})
