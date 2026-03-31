library(acousticTS)

test_that("general utilities resolve scalar and complex helpers correctly", {
  expect_equal(acousticTS:::.calculate_max_radius(0.02, 0.1, NULL), 0.02)
  expect_equal(acousticTS:::.calculate_max_radius(NULL, 0.1, 5), 0.02)
  expect_error(
    acousticTS:::.calculate_max_radius(NULL, 0.1, NULL),
    "Either 'radius' or 'length_radius_ratio' must be provided."
  )

  expect_equal(acousticTS:::`%||%`(3, 9), 3)
  expect_equal(acousticTS:::`%||%`(NA_real_, 9), 9)
  expect_equal(acousticTS:::`%||%`(NULL, 9), 9)
  expect_equal(acousticTS:::`%||%`(list(alpha = 1), 9), list(alpha = 1))

  expect_equal(acousticTS:::`%R%`(1 + 0i), 1)
  expect_equal(
    acousticTS:::`%R%`(c(1 + 0i, 2 + 0i)),
    c(1, 2)
  )
  expect_equal(
    acousticTS:::`%R%`(c(1 + 1i, 2 + 0i)),
    c(1 + 1i, 2 + 0i)
  )
  expect_equal(acousticTS:::`%R%`(5), 5)
})

test_that("simulation parameter resolver handles batch, function, scalar, and vector inputs", {
  simulation_grid <- data.frame(alpha_idx = c(1, 2, 1))
  batch_values <- list(alpha = c(10, 20))

  expect_equal(
    acousticTS:::.resolve_param_value(
      param_name = "alpha",
      param_value = 999,
      batch_by = "alpha",
      batch_values = batch_values,
      grid_size = 3,
      simulation_grid = simulation_grid
    ),
    c(10, 20, 10)
  )

  expect_equal(
    acousticTS:::.resolve_param_value(
      param_name = "beta",
      param_value = function() 7,
      batch_by = NULL,
      batch_values = list(),
      grid_size = 4,
      simulation_grid = data.frame()
    ),
    rep(7, 4)
  )

  structured_scalar <- acousticTS:::.resolve_param_value(
    param_name = "beta_target",
    param_value = c(length = 0.02),
    batch_by = NULL,
    batch_values = list(),
    grid_size = 3,
    simulation_grid = data.frame()
  )
  expect_true(is.list(structured_scalar))
  expect_equal(structured_scalar[[1]], c(length = 0.02))
  expect_equal(structured_scalar[[3]], c(length = 0.02))

  structured_draws <- acousticTS:::.resolve_param_value(
    param_name = "gamma_target",
    param_value = function() c(length = 0.03),
    batch_by = NULL,
    batch_values = list(),
    grid_size = 2,
    simulation_grid = data.frame()
  )
  expect_true(is.list(structured_draws))
  expect_equal(structured_draws[[1]], c(length = 0.03))
  expect_equal(structured_draws[[2]], c(length = 0.03))

  expect_equal(
    acousticTS:::.resolve_param_value(
      param_name = "gamma",
      param_value = 3,
      batch_by = NULL,
      batch_values = list(),
      grid_size = 4,
      simulation_grid = data.frame()
    ),
    rep(3, 4)
  )

  expect_equal(
    acousticTS:::.resolve_param_value(
      param_name = "delta",
      param_value = c(1, 2, 3),
      batch_by = NULL,
      batch_values = list(),
      grid_size = 3,
      simulation_grid = data.frame()
    ),
    c(1, 2, 3)
  )

  expect_error(
    acousticTS:::.resolve_param_value(
      param_name = "epsilon",
      param_value = c(1, 2),
      batch_by = NULL,
      batch_values = list(),
      grid_size = 3,
      simulation_grid = data.frame()
    ),
    "Length of parameter 'epsilon' \\[2\\] does not match number of realizations \\[3\\]."
  )
})

test_that("simulation batch-value normalizer preserves structured candidates", {
  expect_equal(
    acousticTS:::.normalize_simulation_batch_values(
      param_name = "body_target",
      param_value = c(length = 0.02)
    ),
    list(c(length = 0.02))
  )

  expect_equal(
    acousticTS:::.normalize_simulation_batch_values(
      param_name = "body_target",
      param_value = list(c(length = 0.02), c(length = 0.03))
    ),
    list(c(length = 0.02), c(length = 0.03))
  )
})

test_that("extract walks nested list, matrix, and vector paths", {
  obj <- cal_generate()
  obj@model <- list(
    helper = list(
      mat = matrix(
        c(1, 2, 3, 4),
        nrow = 2,
        dimnames = list(c("r1", "r2"), c("c1", "c2"))
      ),
      vec = c(alpha = 10, beta = 20)
    )
  )

  expect_equal(acousticTS::extract(obj, c("model", "helper", "mat", "r1")), c(c1 = 1, c2 = 3))
  expect_equal(acousticTS::extract(obj, c("model", "helper", "mat", "c2")), c(r1 = 3, r2 = 4))
  expect_equal(acousticTS::extract(obj, c("model", "helper", "vec", "beta")), c(beta = 20))

  expect_error(
    acousticTS::extract(obj, c("model", "helper", "missing")),
    "No feature 'missing'"
  )
  expect_error(
    acousticTS::extract(obj, c("model", "helper", "vec", "gamma")),
    "No feature 'gamma'"
  )
})

test_that("validation helpers enforce expected argument contracts", {
  valid_dims <- c("length", "width", "height")

  expect_null(
    acousticTS:::.validate_dimensions_target(
      target = NULL,
      target_name = "body_target",
      valid_dims = valid_dims
    )
  )
  expect_equal(
    acousticTS:::.validate_dimensions_target(
      target = c(length = 2),
      target_name = "body_target",
      valid_dims = valid_dims
    ),
    c(length = 2)
  )
  expect_error(
    acousticTS:::.validate_dimensions_target("x", "body_target", valid_dims),
    "'body_target' must be numeric."
  )
  expect_error(
    acousticTS:::.validate_dimensions_target(c(1, 2), "body_target", valid_dims),
    "'body_target' must either be a scalar or a named vector"
  )
  expect_error(
    acousticTS:::.validate_dimensions_target(c(depth = 2), "body_target", valid_dims),
    "'body_target' has one or more invalid dimensions: 'depth'."
  )

  expect_equal(
    acousticTS:::.validate_dimension_scaling(
      dims = 2,
      dims_name = "body_scale",
      valid_dims = valid_dims,
      isometry = TRUE,
      iso_name = "isometric_body"
    ),
    stats::setNames(rep(2, length(valid_dims)), valid_dims)
  )
  expect_equal(
    acousticTS:::.validate_dimension_scaling(
      dims = c(width = 1.5),
      dims_name = "body_scale",
      valid_dims = valid_dims,
      isometry = FALSE,
      iso_name = "isometric_body"
    ),
    c(length = 1, width = 1.5, height = 1)
  )
  expect_error(
    acousticTS:::.validate_dimension_scaling(
      dims = "x",
      dims_name = "body_scale",
      valid_dims = valid_dims,
      isometry = FALSE,
      iso_name = "isometric_body"
    ),
    "'body_scale' must be numeric."
  )
  expect_error(
    acousticTS:::.validate_dimension_scaling(
      dims = c(length = 2, width = 3),
      dims_name = "body_scale",
      valid_dims = valid_dims,
      isometry = TRUE,
      iso_name = "isometric_body"
    ),
    "'body_scale' contains more than 1 dimension while 'isometric_body' is TRUE."
  )
  expect_error(
    acousticTS:::.validate_dimension_scaling(
      dims = c(depth = 2),
      dims_name = "body_scale",
      valid_dims = valid_dims,
      isometry = FALSE,
      iso_name = "isometric_body"
    ),
    "'body_scale' has one or more invalid dimensions: 'depth'."
  )

  expect_equal(
    acousticTS:::.validate_elastic_inputs(K = 1, E = 2, G = NULL, nu = NULL),
    c("K", "E")
  )
  expect_error(
    acousticTS:::.validate_elastic_inputs(K = 1, E = NULL, G = NULL, nu = NULL),
    "At least two elasticity moduli values are required."
  )
  expect_error(
    acousticTS:::.validate_elastic_inputs(
      K = 1,
      E = NULL,
      G = NULL,
      nu = NULL,
      param_name = "Young's modulus"
    ),
    paste0(
      "At least two elasticity moduli values are required to calculate ",
      "Young's modulus."
    )
  )
})

test_that("plotting and brake validators reject malformed inputs cleanly", {
  expect_invisible(
    acousticTS:::.validate_brake_params(
      body_df = list(rpos = matrix(0, nrow = 4, ncol = 4)),
      radius_curvature = 2,
      mode = "ratio"
    )
  )
  expect_error(
    acousticTS:::.validate_brake_params(
      body_df = list(),
      radius_curvature = 2,
      mode = "ratio"
    ),
    "Body shape information must be a list with a matrix element 'rpos'."
  )
  expect_error(
    acousticTS:::.validate_brake_params(
      body_df = list(rpos = matrix(0, nrow = 4, ncol = 4)),
      radius_curvature = 0,
      mode = "ratio"
    ),
    "Radius of curvature must be a positive-only, real number."
  )
  expect_error(
    acousticTS:::.validate_brake_params(
      body_df = list(rpos = matrix(0, nrow = 4, ncol = 4)),
      radius_curvature = 2,
      mode = "bad"
    ),
    "Radius-of-curvature 'mode' must be either 'ratio' or 'measurement'."
  )

  empty_obj <- cal_generate()
  expect_error(
    acousticTS:::.validate_and_extract_model(empty_obj, "calibration"),
    "ERROR: no model results detected in object."
  )

  empty_obj@model <- list(other = list())
  empty_obj@model_parameters <- list(other = list())
  expect_error(
    acousticTS:::.validate_and_extract_model(empty_obj, "calibration"),
    "Model 'calibration' not found. Available: other"
  )

  good_obj <- cal_generate()
  good_obj@model <- list(calibration = list(TS = 1))
  good_obj@model_parameters <- list(calibration = list(quantity = "TS"))
  info <- acousticTS:::.validate_and_extract_model(good_obj, "calibration")

  expect_equal(info$model$TS, 1)
  expect_equal(info$parameters$quantity, "TS")
  expect_true(is.list(info$shape))
})
