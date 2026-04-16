library(acousticTS)

test_that("internal TMM setup helpers cover boundary, truncation, and geometry branches", {
  sphere_shape <- sphere(radius_body = 0.01, n_segments = 20)
  sphere_params <- fls_generate(shape = sphere_shape, g_body = 1, h_body = 1)@shape_parameters
  prolate_params <- fls_generate(
    shape = prolate_spheroid(length_body = 0.2, radius_body = 0.01, n_segments = 20),
    g_body = 1,
    h_body = 1
  )@shape_parameters
  prolate_short <- fls_generate(
    shape = prolate_spheroid(length_body = 0.04, radius_body = 0.01, n_segments = 20),
    g_body = 1,
    h_body = 1
  )@shape_parameters
  oblate_params <- fls_generate(
    shape = oblate_spheroid(length_body = 0.04, radius_body = 0.05, n_segments = 20),
    g_body = 1,
    h_body = 1
  )@shape_parameters
  cylinder_params <- fls_generate(
    shape = cylinder(length_body = 0.2, radius_body = 0.01, n_segments = 20),
    g_body = 1,
    h_body = 1
  )@shape_parameters

  fls_obj <- fls_generate(
    shape = sphere_shape,
    density_body = 1025,
    sound_speed_body = 1500,
    theta_body = pi / 3
  )
  gas_obj <- gas_generate(
    shape = sphere_shape,
    density_fluid = 1.2,
    sound_speed_fluid = 343,
    theta_body = pi / 4
  )
  shell_obj <- ess_generate(
    shape = sphere_shape,
    radius_shell = 0.01,
    shell_thickness = 0.001,
    density_shell = 1028.9,
    sound_speed_shell = 1480.3,
    density_fluid = 1031,
    sound_speed_fluid = 1483.3,
    theta_shell = pi / 2
  )
  elastic_shell_obj <- fixture_sphere("elastic_shelled")
  elastic_prolate_obj <- fixture_ps("elastic_shelled")

  expect_equal(acousticTS:::.tmm_boundary_default(fls_obj, NULL), "liquid_filled")
  expect_equal(acousticTS:::.tmm_boundary_default(gas_obj, NULL), "gas_filled")
  expect_equal(acousticTS:::.tmm_boundary_default(shell_obj, NULL), "shelled_liquid")
  expect_equal(
    acousticTS:::.tmm_boundary_default(elastic_shell_obj, NULL),
    "elastic_shelled"
  )
  expect_error(
    acousticTS:::.tmm_boundary_default(elastic_prolate_obj, NULL),
    "Specify 'boundary' explicitly"
  )
  expect_equal(
    acousticTS:::.tmm_boundary_default(fls_obj, "fixed_rigid"),
    "fixed_rigid"
  )
  expect_error(
    acousticTS:::.tmm_boundary_default(cal_generate(), NULL),
    "Specify 'boundary' explicitly"
  )

  expect_invisible(acousticTS:::.tmm_validate_object_scope(fls_obj))
  expect_invisible(acousticTS:::.tmm_validate_object_scope(shell_obj))
  expect_error(
    acousticTS:::.tmm_validate_object_scope(cal_generate()),
    "requires the scatterer to be either 'FLS', 'GAS', or a supported 'ESS'"
  )

  expect_invisible(
    acousticTS:::.tmm_require_homogeneous_body(
      list(density = c(1, 2), sound_speed = 1500),
      "fixed_rigid"
    )
  )
  expect_error(
    acousticTS:::.tmm_require_homogeneous_body(
      list(density = c(1, 2), sound_speed = 1500),
      "liquid_filled"
    ),
    "property 'density'"
  )
  expect_error(
    acousticTS:::.tmm_resolve_boundary(fls_obj, "elastic"),
    "Only the following values for 'boundary' are available"
  )
  expect_equal(
    acousticTS:::.tmm_resolve_boundary(shell_obj, "shelled_gas"),
    "shelled_gas"
  )
  expect_equal(
    acousticTS:::.tmm_resolve_boundary(elastic_shell_obj, "elastic_shelled"),
    "elastic_shelled"
  )
  expect_equal(
    acousticTS:::.tmm_resolve_boundary(elastic_prolate_obj, "elastic_shelled"),
    "elastic_shelled"
  )
  expect_error(
    acousticTS:::.tmm_validate_store_t_matrix(c(TRUE, FALSE)),
    "'store_t_matrix' must be either TRUE or FALSE"
  )
  expect_error(
    acousticTS:::.tmm_validate_penetrable_body(list(density = NULL), "liquid_filled"),
    "Penetrable TMM boundaries require scalar body density and sound speed"
  )
  expect_invisible(
    acousticTS:::.tmm_validate_penetrable_body(
      list(density = 1025, sound_speed = 1500),
      "liquid_filled"
    )
  )

  expect_equal(acousticTS:::.tmm_bounding_radius(sphere_params), 0.01)
  expect_equal(acousticTS:::.tmm_bounding_radius(prolate_params), 0.1)
  expect_equal(acousticTS:::.tmm_bounding_radius(oblate_params), 0.05)
  expect_equal(
    acousticTS:::.tmm_bounding_radius(cylinder_params),
    sqrt(0.1^2 + 0.01^2),
    tolerance = 1e-12
  )

  expect_true(is.na(acousticTS:::.tmm_prolate_nmax_override(sphere_params, "fixed_rigid")))
  expect_true(is.na(acousticTS:::.tmm_prolate_nmax_override(prolate_short, "fixed_rigid")))
  expect_equal(acousticTS:::.tmm_prolate_nmax_override(prolate_params, "fixed_rigid"), 30L)
  expect_equal(acousticTS:::.tmm_prolate_nmax_override(prolate_params, "pressure_release"), 24L)
  expect_equal(acousticTS:::.tmm_prolate_nmax_override(prolate_params, "liquid_filled"), 36L)
  expect_equal(acousticTS:::.tmm_prolate_nmax_override(prolate_params, "gas_filled"), 36L)

  expect_true(is.na(acousticTS:::.tmm_cylinder_nmax_floor(sphere_params, "fixed_rigid")))
  expect_equal(acousticTS:::.tmm_cylinder_nmax_floor(cylinder_params, "fixed_rigid"), 56L)
  expect_equal(acousticTS:::.tmm_cylinder_nmax_floor(cylinder_params, "pressure_release"), 52L)
  expect_equal(acousticTS:::.tmm_cylinder_nmax_floor(cylinder_params, "liquid_filled"), 64L)
  expect_equal(acousticTS:::.tmm_cylinder_nmax_floor(cylinder_params, "gas_filled"), 56L)
  expect_true(is.na(acousticTS:::.tmm_oblate_nmax_floor(sphere_params, "fixed_rigid")))
  expect_equal(acousticTS:::.tmm_oblate_nmax_floor(oblate_params, "fixed_rigid"), 24L)
  expect_equal(acousticTS:::.tmm_oblate_nmax_floor(oblate_params, "liquid_filled"), 24L)
  expect_true(is.na(acousticTS:::.tmm_oblate_nmax_floor(oblate_params, "pressure_release")))

  expect_equal(
    acousticTS:::.tmm_prepare_n_max(
      NULL,
      frequency = c(1, 2),
      k_sw = c(1, 2),
      shape_parameters = prolate_params,
      boundary = "liquid_filled"
    ),
    c(36L, 36L)
  )
  expect_true(all(
    acousticTS:::.tmm_prepare_n_max(
      NULL,
      frequency = c(1, 2),
      k_sw = c(0.5, 1),
      shape_parameters = cylinder_params,
      boundary = "fixed_rigid"
    ) >= 56L
  ))
  expect_equal(
    acousticTS:::.tmm_prepare_n_max(
      NULL,
      frequency = c(1, 2),
      k_sw = c(0.5, 1),
      shape_parameters = oblate_params,
      boundary = "fixed_rigid"
    ),
    c(24L, 24L)
  )
  expect_equal(
    acousticTS:::.tmm_prepare_n_max(
      NULL,
      frequency = c(1, 2),
      k_sw = c(0.5, 1),
      shape_parameters = oblate_params,
      boundary = "liquid_filled"
    ),
    c(24L, 24L)
  )
  expect_equal(
    acousticTS:::.tmm_prepare_n_max(
      7,
      frequency = c(1, 2),
      k_sw = c(1, 2),
      shape_parameters = sphere_params,
      boundary = "fixed_rigid"
    ),
    c(7L, 7L)
  )
  expect_equal(
    acousticTS:::.tmm_prepare_n_max(
      c(5, 6),
      frequency = c(1, 2),
      k_sw = c(1, 2),
      shape_parameters = sphere_params,
      boundary = "fixed_rigid"
    ),
    c(5L, 6L)
  )
  expect_error(
    acousticTS:::.tmm_prepare_n_max(
      0,
      frequency = c(1, 2),
      k_sw = c(1, 2),
      shape_parameters = sphere_params,
      boundary = "fixed_rigid"
    ),
    "'n_max' must be NULL or a positive finite integer vector"
  )
  expect_error(
    acousticTS:::.tmm_prepare_n_max(
      c(5, 6, 7),
      frequency = c(1, 2),
      k_sw = c(1, 2),
      shape_parameters = sphere_params,
      boundary = "fixed_rigid"
    ),
    "must be length 1 or match the length of 'frequency'"
  )

  expect_equal(
    acousticTS:::.tmm_prepare_cylinder_n_max(
      NULL,
      frequency = c(1, 2),
      k_sw = c(100, 200),
      shape_parameters = cylinder_params
    ),
    c(11L, 12L)
  )
  expect_equal(
    acousticTS:::.tmm_prepare_cylinder_n_max(
      4,
      frequency = c(1, 2),
      k_sw = c(1, 2),
      shape_parameters = cylinder_params
    ),
    c(4L, 4L)
  )
  expect_equal(
    acousticTS:::.tmm_prepare_cylinder_n_max(
      c(4, 5),
      frequency = c(1, 2),
      k_sw = c(1, 2),
      shape_parameters = cylinder_params
    ),
    c(4L, 5L)
  )
  expect_error(
    acousticTS:::.tmm_prepare_cylinder_n_max(
      -1,
      frequency = c(1, 2),
      k_sw = c(1, 2),
      shape_parameters = cylinder_params
    ),
    "'n_max' must be NULL or a positive finite integer vector"
  )
  expect_error(
    acousticTS:::.tmm_prepare_cylinder_n_max(
      c(4, 5, 6),
      frequency = c(1, 2),
      k_sw = c(1, 2),
      shape_parameters = cylinder_params
    ),
    "must be length 1 or match the length of 'frequency'"
  )

  old_opt <- options(acousticTS.tmm_cylinder_backend = "spherical")
  on.exit(options(old_opt), add = TRUE)
  expect_false(acousticTS:::.tmm_is_cylindrical_branch(cylinder_params))
  expect_equal(
    acousticTS:::.tmm_prepare_cylinder_spherical_n_max(
      NULL,
      frequency = c(1, 2),
      k_sw = c(100, 200),
      shape_parameters = cylinder_params,
      boundary = "fixed_rigid"
    ),
    c(47L, 48L)
  )
  expect_equal(
    acousticTS:::.tmm_prepare_cylinder_spherical_n_max(
      NULL,
      frequency = c(1, 2),
      k_sw = c(100, 200),
      shape_parameters = cylinder_params,
      boundary = "pressure_release"
    ),
    c(37L, 38L)
  )
  expect_equal(
    acousticTS:::.tmm_prepare_cylinder_spherical_n_max(
      NULL,
      frequency = c(1, 2),
      k_sw = c(100, 200),
      shape_parameters = cylinder_params,
      boundary = "gas_filled"
    ),
    c(16L, 17L)
  )
  expect_equal(
    acousticTS:::.tmm_prepare_cylinder_spherical_n_max(
      NULL,
      frequency = c(1, 2),
      k_sw = c(100, 200),
      shape_parameters = cylinder_params,
      boundary = "liquid_filled"
    ),
    c(5L, 6L)
  )

  expect_equal(acousticTS:::.tmm_collocation_nodes(cylinder_params, "fixed_rigid", 3), 128L)
  expect_equal(acousticTS:::.tmm_collocation_nodes(sphere_params, "fixed_rigid", 3), 64L)

  expect_equal(acousticTS:::.tmm_shape_values(sphere_params), 0.01)
  expect_equal(acousticTS:::.tmm_shape_values(prolate_params), c(0.1, 0.01))
  expect_equal(acousticTS:::.tmm_shape_values(oblate_params), c(0.02, 0.05))
  expect_equal(acousticTS:::.tmm_shape_values(cylinder_params), c(0.1, 0.01))
  expect_error(
    acousticTS:::.tmm_shape_values(list(shape = "Weird")),
    "Unsupported TMM shape geometry"
  )

  acoustics <- data.frame(
    frequency = 38000,
    k_sw = wavenumber(38000, 1500),
    n_max = 8L
  )
  body <- list(theta_body = pi / 3, g_body = 1.05, h_body = 0.97)
  rigid_out <- acousticTS:::.tmm_run_cylindrical_branch(
    shape_parameters = cylinder_params,
    acoustics = acoustics,
    body = body,
    boundary = "fixed_rigid"
  )
  fluid_out <- acousticTS:::.tmm_run_cylindrical_branch(
    shape_parameters = cylinder_params,
    acoustics = acoustics,
    body = body,
    boundary = "liquid_filled"
  )
  expect_equal(nrow(rigid_out$model), 1L)
  expect_equal(nrow(fluid_out$model), 1L)
  expect_true(all(is.finite(rigid_out$model$TS)))
  expect_true(all(is.finite(fluid_out$model$TS)))
})

test_that("internal TMM spherical helpers cover surface and solver branches", {
  sphere_params <- fls_generate(
    shape = sphere(radius_body = 0.01, n_segments = 20),
    g_body = 1,
    h_body = 1
  )@shape_parameters
  prolate_params <- fls_generate(
    shape = prolate_spheroid(length_body = 0.2, radius_body = 0.01, n_segments = 20),
    g_body = 1,
    h_body = 1
  )@shape_parameters
  oblate_params <- fls_generate(
    shape = oblate_spheroid(length_body = 0.04, radius_body = 0.05, n_segments = 20),
    g_body = 1,
    h_body = 1
  )@shape_parameters
  cylinder_params <- fls_generate(
    shape = cylinder(length_body = 0.2, radius_body = 0.01, n_segments = 20),
    g_body = 1,
    h_body = 1
  )@shape_parameters
  tapered_cylinder_params <- fls_generate(
    shape = cylinder(length_body = 0.2, radius_body = 0.01, taper = 8, n_segments = 20),
    g_body = 1,
    h_body = 1
  )@shape_parameters
  theta <- c(0, pi / 4, pi / 2)

  expect_equal(acousticTS:::.tmm_double_factorial_odd(0), 1)
  expect_equal(acousticTS:::.tmm_double_factorial_odd(3), 15)

  sphere_radius <- acousticTS:::.tmm_surface_radius(sphere_params, theta)
  prolate_radius <- acousticTS:::.tmm_surface_radius(prolate_params, theta)
  oblate_radius <- acousticTS:::.tmm_surface_radius(oblate_params, theta)
  cylinder_radius <- acousticTS:::.tmm_surface_radius(cylinder_params, theta)
  tapered_radius <- acousticTS:::.tmm_surface_radius(tapered_cylinder_params, theta)

  expect_equal(sphere_radius$radius, rep(0.01, length(theta)))
  expect_equal(sphere_radius$radius_derivative, rep(0, length(theta)))
  expect_true(all(prolate_radius$radius > 0))
  expect_true(all(oblate_radius$radius > 0))
  expect_true(all(cylinder_radius$radius > 0))
  expect_true(all(tapered_radius$radius > 0))
  expect_equal(length(cylinder_radius$radius_derivative), length(theta))
  expect_equal(acousticTS:::.tmm_cylinder_endcap_fraction(tapered_cylinder_params), 0)
  expect_equal(acousticTS:::.tmm_cylinder_endcap_fraction(cylinder_params), 0)
  expect_equal(
    acousticTS:::.tmm_cylinder_endcap_fraction(
      cylinder_params,
      cylinder_endcap_fraction = 0.05
    ),
    0.05
  )
  expect_equal(acousticTS:::.tmm_resolve_cylinder_endcap_fraction(0.6), 0.45)
  expect_error(
    acousticTS:::.tmm_resolve_cylinder_endcap_fraction(-0.01),
    "'cylinder_endcap_fraction' must be non-negative"
  )
  piecewise_surface <- acousticTS:::.tmm_cylinder_piecewise_surface(
    cylinder_params,
    n_terms = 6
  )
  expect_true(all(piecewise_surface$radius > 0))
  expect_true(all(piecewise_surface$area_weight > 0))
  expect_equal(length(piecewise_surface$mu), length(piecewise_surface$normal_r))
  expect_equal(length(piecewise_surface$mu), length(piecewise_surface$normal_theta))
  expect_error(
    acousticTS:::.tmm_surface_radius(list(shape = "Weird"), theta),
    "Unsupported TMM shape geometry"
  )

  singular_square <- matrix(c(1, 1, 2, 2), nrow = 2)
  rhs_square <- c(1, 2)
  expect_error(
    acousticTS:::.tmm_solve_linear_system(singular_square, rhs_square),
    "singular matrix"
  )

  rectangular <- matrix(c(1, 0, 0, 1, 1, 1), nrow = 3, byrow = TRUE)
  rhs_rect <- c(1, 1, 2)
  expect_equal(
    acousticTS:::.tmm_solve_linear_system(rectangular, rhs_rect),
    qr.solve(rectangular, rhs_rect)
  )

  mu <- cos(c(pi / 6, pi / 3))
  p_mat <- acousticTS:::.tmm_assoc_legendre_table(0, 1, mu)
  dp_dtheta <- acousticTS:::.tmm_assoc_legendre_theta_derivative(0, 0:1, mu, p_mat)
  radial <- matrix(c(1 + 0i, 2 + 0i, 3 + 0i, 4 + 0i), nrow = 2, byrow = TRUE)
  radial_deriv <- matrix(c(0.5 + 0i, 0.25 + 0i, 0.75 + 0i, 0.1 + 0i), nrow = 2, byrow = TRUE)
  radius <- c(1.2, 0.9)
  radius_derivative <- c(0.3, -0.25)

  normal_scale <- sqrt(1 + (radius_derivative / radius)^2)
  normal_r <- 1 / normal_scale
  normal_theta <- (radius_derivative / radius) / normal_scale

  expect_equal(
    acousticTS:::.tmm_normal_derivative_matrix(
      radial = radial,
      radial_deriv = radial_deriv,
      angular = p_mat,
      angular_theta_deriv = dp_dtheta,
      k = 2.3,
      radius = radius,
      radius_derivative = radius_derivative
    ),
    acousticTS:::.tmm_normal_derivative_explicit(
      radial = radial,
      radial_deriv = radial_deriv,
      angular = p_mat,
      angular_theta_deriv = dp_dtheta,
      k = 2.3,
      radius = radius,
      normal_r = normal_r,
      normal_theta = normal_theta
    ),
    tolerance = 1e-12
  )
})

test_that("internal TMM diagnostics helpers cover validation and continuation branches", {
  density_sw <- 1026.8
  sound_speed_sw <- 1477.3

  sphere_obj <- target_strength(
    object = fls_generate(
      shape = sphere(radius_body = 0.01, n_segments = 40),
      g_body = 1,
      h_body = 1,
      theta_body = pi / 2
    ),
    frequency = c(38e3, 70e3),
    model = "tmm",
    boundary = "pressure_release",
    density_sw = density_sw,
    sound_speed_sw = sound_speed_sw,
    store_t_matrix = TRUE
  )
  oblate_obj <- target_strength(
    object = fls_generate(
      shape = oblate_spheroid(length_body = 0.03, radius_body = 0.05, n_segments = 40),
      g_body = 1,
      h_body = 1,
      theta_body = pi / 2
    ),
    frequency = 38e3,
    model = "tmm",
    boundary = "pressure_release",
    density_sw = density_sw,
    sound_speed_sw = sound_speed_sw,
    store_t_matrix = TRUE
  )

  expect_equal(acousticTS:::.tmm_frequency_indices(NULL, c(38e3, 70e3)), 1:2)
  expect_equal(
    acousticTS:::.tmm_frequency_indices(c(38e3, 70e3), c(38e3, 70e3)),
    1:2
  )
  expect_error(
    acousticTS:::.tmm_frequency_indices("bad", c(38e3, 70e3)),
    "'frequency' must be NULL or a numeric vector of finite values in Hz"
  )

  default_pairs <- acousticTS:::.tmm_validate_reciprocity_pairs(NULL)
  expect_true(all(
    c("theta_body", "phi_body", "theta_scatter", "phi_scatter") %in%
      names(default_pairs)
  ))
  expect_error(
    acousticTS:::.tmm_validate_reciprocity_pairs(data.frame(theta_body = 1)),
    "must contain the columns"
  )
  expect_error(
    acousticTS:::.tmm_validate_reciprocity_pairs(data.frame(
      theta_body = "bad",
      phi_body = 0,
      theta_scatter = 0,
      phi_scatter = 0
    )),
    "must be finite numeric angles"
  )

  sphere_params <- fls_generate(
    shape = sphere(radius_body = 0.01, n_segments = 20),
    g_body = 1,
    h_body = 1
  )@shape_parameters
  oblate_params <- fls_generate(
    shape = oblate_spheroid(length_body = 0.04, radius_body = 0.05, n_segments = 20),
    g_body = 1,
    h_body = 1
  )@shape_parameters
  expect_equal(acousticTS:::.tmm_equivalent_volume_radius(sphere_params), 0.01)
  expect_equal(
    acousticTS:::.tmm_equivalent_volume_radius(oblate_params),
    (0.05^2 * 0.02)^(1 / 3),
    tolerance = 1e-12
  )
  expect_equal(acousticTS:::.tmm_spheroidal_aspect_ratio(oblate_params), 2.5)

  axes_oblate <- acousticTS:::.tmm_spheroid_axes_from_equal_volume(
    "OblateSpheroid",
    radius_eq = 1,
    aspect_ratio = 2
  )
  expect_equal(axes_oblate$length_body, 2 / (2^(2 / 3)), tolerance = 1e-12)
  expect_equal(axes_oblate$radius_body, 2^(1 / 3), tolerance = 1e-12)
  expect_error(
    acousticTS:::.tmm_spheroid_axes_from_equal_volume("Cylinder", 1, 2),
    "Unsupported shape for spheroidal continuation"
  )

  gas_obj <- gas_generate(
    shape = sphere(radius_body = 0.01, n_segments = 20),
    density_fluid = 1.2,
    sound_speed_fluid = 343,
    theta_body = pi / 4
  )
  gas_rebuilt <- acousticTS:::.tmm_rebuild_shape_like(
    gas_obj,
    shape = sphere(radius_body = 0.02, n_segments = 20),
    body = list(theta = pi / 5, density = 1.3, sound_speed = 345)
  )
  fls_rebuilt <- acousticTS:::.tmm_rebuild_shape_like(
    fls_generate(
      shape = sphere(radius_body = 0.01, n_segments = 20),
      g_body = 0.9,
      h_body = 1.1,
      theta_body = pi / 6
    ),
    shape = sphere(radius_body = 0.02, n_segments = 20),
    body = list(theta = pi / 7, g = 0.8, h = 1.2)
  )
  expect_s4_class(gas_rebuilt, "GAS")
  expect_s4_class(fls_rebuilt, "FLS")
  expect_equal(extract(gas_rebuilt, "body")$density, 1.3)
  expect_equal(extract(fls_rebuilt, "body")$g, 0.8)

  sphere_model_params <- sphere_obj@model_parameters$TMM
  oblate_model_params <- oblate_obj@model_parameters$TMM
  expect_null(acousticTS:::.tmm_sphere_to_spheroid_path(
    sphere_obj,
    sphere_model_params,
    continuation_steps = 3L
  ))
  expect_null(acousticTS:::.tmm_sphere_to_spheroid_path(
    oblate_obj,
    oblate_model_params,
    continuation_steps = 1L
  ))
  expect_error(
    acousticTS:::.tmm_sphere_to_spheroid_path(
      oblate_obj,
      oblate_model_params,
      continuation_steps = 1.5
    ),
    "'continuation_steps' must be a single integer >= 0"
  )
  oblate_path <- acousticTS:::.tmm_sphere_to_spheroid_path(
    oblate_obj,
    oblate_model_params,
    continuation_steps = 3L
  )
  expect_true(is.data.frame(oblate_path))
  expect_true(any(oblate_path$shape == "Sphere"))
  expect_true(any(oblate_path$shape == "OblateSpheroid"))

  one_step_summary <- acousticTS:::.tmm_sphere_to_spheroid_summary(data.frame(
    step = 1,
    shape = "Sphere",
    aspect_ratio = 1,
    target_aspect_ratio = 1,
    frequency = 38e3,
    sigma_bs = 1,
    TS = -50
  ))
  expect_equal(one_step_summary$continuation_max_abs_step_TS, 0)
  expect_equal(one_step_summary$continuation_max_abs_second_diff_TS, 0)

  expect_equal(nrow(acousticTS:::.tmm_block_metrics(NULL)), 0L)
  zero_block <- acousticTS:::.tmm_block_metrics(list(list(
    T = matrix(0 + 0i, nrow = 2, ncol = 2),
    m = 0,
    n_seq = 0:1
  )))
  expect_equal(zero_block$transpose_residual, 0)

  expect_error(acousticTS:::.tmm_validate_diagnostic_grid(4, 5), "single integer >= 5")
  expect_error(acousticTS:::.tmm_validate_diagnostic_grid(5, 4), "single integer >= 5")

  expect_true(all(is.na(unlist(acousticTS:::.tmm_continuation_metrics(NULL, 38e3)))))
  expect_true(all(is.na(unlist(acousticTS:::.tmm_continuation_metrics(
    data.frame(
      frequency = 70e3,
      continuation_target_aspect_ratio = 2,
      continuation_max_abs_step_TS = 1,
      continuation_max_abs_second_diff_TS = 0.1,
      continuation_any_nonfinite = FALSE
    ),
    38e3
  )))))
})

test_that("internal TMM scattering and orientation helpers cover remaining validation branches", {
  density_sw <- 1026.8
  sound_speed_sw <- 1477.3

  raw_obj <- fls_generate(
    shape = sphere(radius_body = 0.01, n_segments = 40),
    g_body = 1,
    h_body = 1,
    theta_body = pi / 2
  )
  stored_obj <- target_strength(
    object = raw_obj,
    frequency = c(38e3, 70e3),
    model = "tmm",
    boundary = "pressure_release",
    density_sw = density_sw,
    sound_speed_sw = sound_speed_sw,
    store_t_matrix = TRUE
  )
  storeless_obj <- target_strength(
    object = raw_obj,
    frequency = 38e3,
    model = "tmm",
    boundary = "pressure_release",
    density_sw = density_sw,
    sound_speed_sw = sound_speed_sw,
    store_t_matrix = FALSE
  )
  sphere_params <- extract(sphere(radius_body = 0.01, n_segments = 20), "shape_parameters")
  cylinder_params <- extract(
    cylinder(length_body = 0.2, radius_body = 0.01, n_segments = 20),
    "shape_parameters"
  )

  expect_error(
    acousticTS:::.tmm_require_stored_blocks(raw_obj),
    "No stored TMM model was found"
  )
  expect_error(
    acousticTS:::.tmm_require_stored_blocks(storeless_obj),
    "Stored T-matrix blocks are required for this helper"
  )
  expect_equal(acousticTS:::.tmm_scalar_angle(NULL, pi / 3, "theta_body"), pi / 3)
  expect_error(
    acousticTS:::.tmm_scalar_angle(c(1, 2), pi / 3, "theta_body"),
    "'theta_body' must be a single finite angle"
  )

  expect_equal(
    acousticTS:::.tmm_angle_vector(
      NULL,
      default = c(0, pi),
      n_default = 3,
      lower = 0,
      upper = pi,
      name = "theta_scatter"
    ),
    c(0, pi / 2, pi)
  )
  expect_error(
    acousticTS:::.tmm_angle_vector(NULL, name = "theta_scatter"),
    "'theta_scatter' could not be resolved"
  )
  expect_error(
    acousticTS:::.tmm_angle_vector("bad", name = "theta_scatter"),
    "must be a non-empty numeric vector of finite angles"
  )
  expect_error(
    acousticTS:::.tmm_angle_vector(c(-0.1, 0.2), lower = 0, name = "theta_scatter"),
    "must be >= 0 radians"
  )
  expect_error(
    acousticTS:::.tmm_angle_vector(c(0.1, pi + 0.1), upper = pi, name = "theta_scatter"),
    "must be <="
  )

  expect_equal(acousticTS:::.tmm_interval_weights(0.5), 1)
  expect_equal(
    acousticTS:::.tmm_interval_weights(c(0, pi / 2, pi), lower = 0, upper = pi),
    c(pi / 4, pi / 2, pi / 4)
  )
  expect_error(
    acousticTS:::.tmm_interval_weights(c(1, 0.5)),
    "must be strictly increasing"
  )
  expect_equal(acousticTS:::.tmm_grid_edges(0.5, 0, 1), c(0, 1))
  expect_equal(
    acousticTS:::.tmm_grid_edges(c(1, 2, 4), 0, 5),
    c(0, 1.5, 3, 5)
  )

  expect_error(
    acousticTS:::.tmm_validate_scattering_grid_dims(1, 2),
    "'n_theta' must be a single integer >= 2"
  )
  expect_error(
    acousticTS:::.tmm_validate_scattering_grid_dims(2, 1),
    "'n_phi' must be a single integer >= 2"
  )

  expect_error(
    acousticTS:::.tmm_cylindrical_monostatic_f_bs(
      acoustics_row = data.frame(frequency = 38000, k_sw = wavenumber(38000, 1500), n_max = 8L),
      body_defaults = list(g_body = 1, h_body = 1, theta_body = pi / 2),
      shape_parameters = cylinder_params,
      boundary = "elastic",
      theta_body = pi / 2
    ),
    "Unsupported boundary for cylindrical TMM branch"
  )

  fake_params <- list(
    parameters = list(
      coordinate_system = "weird",
      acoustics = data.frame(frequency = 38000, k_sw = 1),
      t_matrix = list(list())
    )
  )
  expect_error(
    acousticTS:::.tmm_scattering_points(
      model_params = list(
        parameters = list(
          coordinate_system = "spherical",
          acoustics = data.frame(frequency = 38000, k_sw = 1),
          t_matrix = list(list())
        )
      ),
      frequency_idx = 1L,
      shape_parameters = sphere_params,
      theta_body = c(pi / 2, pi / 3),
      phi_body = pi,
      theta_scatter = c(pi / 2, pi / 3),
      phi_scatter = c(pi, 0)
    ),
    "require equal-length angle vectors"
  )
  expect_error(
    acousticTS:::.tmm_scattering_grid_matrix(
      model_params = fake_params,
      frequency_idx = 1L,
      shape_parameters = sphere_params,
      theta_body = pi / 2,
      phi_body = pi,
      theta_scatter = c(0, pi),
      phi_scatter = c(0, pi)
    ),
    "Unsupported stored TMM coordinate system"
  )

  expect_equal(acousticTS:::.tmm_plot_frequency_index(NULL, 38000), 1L)
  expect_error(
    acousticTS:::.tmm_plot_frequency_index(Inf, c(38000, 70000)),
    "'frequency' must be a single finite value in Hz"
  )
  expect_warning(
    idx <- acousticTS:::.tmm_plot_frequency_index(38100, c(38000, 70000)),
    "nearest stored frequency"
  )
  expect_equal(idx, 1L)

  valid_dist <- structure(
    data.frame(
      theta_body = c(0, pi / 2),
      phi_body = c(pi, pi),
      weights = c(2, 1),
      distribution = "quadrature"
    ),
    class = c("TMMOrientationDistribution", "data.frame")
  )
  normalized <- acousticTS:::.tmm_validate_orientation_distribution(valid_dist)
  expect_equal(sum(normalized$weights), 1)
  expect_error(
    acousticTS:::.tmm_validate_orientation_distribution(data.frame(theta_body = 0)),
    "must be created by 'tmm_orientation_distribution"
  )
  expect_error(
    acousticTS:::.tmm_validate_orientation_distribution(structure(
      data.frame(theta_body = 0, phi_body = 0, weights = 1),
      class = c("TMMOrientationDistribution", "data.frame")
    )),
    "missing required column"
  )
  expect_error(
    acousticTS:::.tmm_validate_orientation_distribution(structure(
      data.frame(theta_body = Inf, phi_body = 0, weights = 1, distribution = "quadrature"),
      class = c("TMMOrientationDistribution", "data.frame")
    )),
    "theta_body"
  )
  expect_error(
    acousticTS:::.tmm_validate_orientation_distribution(structure(
      data.frame(theta_body = 0, phi_body = Inf, weights = 1, distribution = "quadrature"),
      class = c("TMMOrientationDistribution", "data.frame")
    )),
    "phi_body"
  )
  expect_error(
    acousticTS:::.tmm_validate_orientation_distribution(structure(
      data.frame(theta_body = 0, phi_body = 0, weights = -1, distribution = "quadrature"),
      class = c("TMMOrientationDistribution", "data.frame")
    )),
    "weights"
  )

  expect_equal(
    acousticTS:::.tmm_distribution_weights(c(0, pi / 2), c(1, 1), lower = 0, upper = pi / 2),
    c(0.5, 0.5)
  )
  expect_error(
    acousticTS:::.tmm_distribution_weights(c(0, pi / 2), c(1, -1)),
    "'density_values' must be a non-negative numeric vector"
  )
  expect_error(
    acousticTS:::.tmm_distribution_weights(c(0, pi / 2), c(0, 0)),
    "Orientation weights must sum to a positive value"
  )

  expect_error(
    tmm_orientation_distribution(distribution = "quadrature"),
    "'theta_body' must be supplied for 'distribution = \"quadrature\"'"
  )
  expect_error(
    tmm_orientation_distribution(distribution = "pdf", pdf = function(x) x),
    "'theta_body' must be supplied for 'distribution = \"pdf\"'"
  )

  expect_error(
    acousticTS:::.tmm_average_orientation_direct_inputs("bad", NULL, pi),
    "'theta_body' must be a non-empty numeric vector of angles in radians"
  )
  expect_error(
    acousticTS:::.tmm_average_orientation_direct_inputs(c(pi / 2), NULL, NA_real_),
    "'phi_body' must be finite"
  )
  expect_error(
    acousticTS:::.tmm_average_orientation_weights(c(1, -1), 2, c(pi / 4, pi / 3)),
    "'weights' must be a non-negative numeric vector"
  )

  default_scatter <- acousticTS:::.tmm_average_scatter_angles(
    theta_body = c(pi / 4, pi / 3),
    phi_body = c(pi, pi / 2),
    theta_scatter = NULL,
    phi_scatter = NULL
  )
  explicit_scatter <- acousticTS:::.tmm_average_scatter_angles(
    theta_body = c(pi / 4, pi / 3),
    phi_body = c(pi, pi / 2),
    theta_scatter = 0.1,
    phi_scatter = 0.2
  )
  expect_equal(default_scatter$theta_scatter, c(3 * pi / 4, 2 * pi / 3))
  expect_equal(explicit_scatter$theta_scatter, c(0.1, 0.1))
  expect_equal(explicit_scatter$phi_scatter, c(0.2, 0.2))

  expect_equal(acousticTS:::.tmm_resolve_orientation_phi(pi / 4, 3), rep(pi / 4, 3))
  expect_error(
    acousticTS:::.tmm_resolve_orientation_phi(c(pi / 4, Inf), 2),
    "'phi_body' must be a finite numeric scalar or vector"
  )

  avg <- tmm_average_orientation(
    object = stored_obj,
    theta_body = pi / 2,
    weights = 1,
    phi_body = pi
  )
  mono <- tmm_scattering(stored_obj, theta_body = pi / 2, phi_body = pi)
  expect_equal(avg$sigma_bs, mono$sigma_scat, tolerance = 1e-10)
})
