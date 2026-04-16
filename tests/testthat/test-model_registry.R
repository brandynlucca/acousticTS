library(acousticTS)

capture_model_registry_state <- function() {
  state <- get(".model_registry_state", envir = asNamespace("acousticTS"))
  old_user <- state$user
  old_loaded <- state$loaded

  function() {
    state$user <- old_user
    state$loaded <- old_loaded
  }
}

tsl_initialize <- function(object,
                           frequency,
                           intercept = -70,
                           slope = 20) {
  shape <- acousticTS::extract(object, "shape_parameters")

  if (is.null(shape$length) || is.na(shape$length)) {
    stop("TSL requires the target shape to have a defined length.")
  }

  methods::slot(object, "model_parameters")$TSL <- list(
    parameters = data.frame(frequency = frequency),
    body = data.frame(length_m = shape$length),
    coefficients = data.frame(intercept = intercept, slope = slope)
  )
  methods::slot(object, "model")$TSL <- data.frame(
    frequency = frequency,
    f_bs = rep(NA_real_, length(frequency)),
    sigma_bs = rep(NA_real_, length(frequency)),
    TS = rep(NA_real_, length(frequency))
  )

  object
}

TSL <- function(object) {
  model <- acousticTS::extract(object, "model_parameters")$TSL
  length_mm <- model$body$length_m * 1e3
  intercept <- model$coefficients$intercept
  slope <- model$coefficients$slope

  TS <- intercept + slope * log10(length_mm)
  sigma_bs <- acousticTS::linear(TS)

  methods::slot(object, "model")$TSL <- data.frame(
    frequency = model$parameters$frequency,
    f_bs = rep(sqrt(sigma_bs), nrow(model$parameters)),
    sigma_bs = rep(sigma_bs, nrow(model$parameters)),
    TS = rep(TS, nrow(model$parameters))
  )

  object
}

test_that("available_models lists built-ins and target_strength resolves aliases", {
  models <- acousticTS::available_models()

  expect_true("calibration" %in% models$model)
  expect_true(any(models$model == "calibration" & grepl("soems", models$aliases)))
  expect_true("espsms" %in% models$model)
  expect_true(any(models$model == "espsms" & grepl("epsms", models$aliases)))

  cal_obj <- target_strength(cal_generate(), frequency = 38e3, model = "soems")

  expect_true("calibration" %in% names(cal_obj@model))
  expect_true(all(is.finite(cal_obj@model$calibration$TS)))
})

test_that("user-registered models work in target_strength and simulate_ts", {
  restore_registry <- capture_model_registry_state()
  on.exit(restore_registry(), add = TRUE)

  acousticTS::register_model(
    name = "tsl",
    initialize = tsl_initialize,
    solver = TSL,
    slot = "TSL",
    aliases = "toy_tsl"
  )

  models <- acousticTS::available_models()
  tsl_row <- models[models$model == "tsl", , drop = FALSE]

  expect_equal(nrow(tsl_row), 1)
  expect_equal(tsl_row$source, "user")
  expect_match(tsl_row$aliases, "toy_tsl")

  obj <- fls_generate(
    shape = cylinder(length_body = 0.07, radius_body = 0.01, n_segments = 80),
    density_body = 1028.9,
    sound_speed_body = 1480.3
  )

  out <- target_strength(
    object = obj,
    frequency = c(38e3, 70e3),
    model = "toy_tsl",
    model_args = list(toy_tsl = list(intercept = -68, slope = 19.5))
  )

  expect_true("TSL" %in% names(out@model))
  expect_true(all(is.finite(out@model$TSL$TS)))
  expect_equal(length(unique(out@model$TSL$TS)), 1)

  sim <- simulate_ts(
    object = obj,
    frequency = c(38e3, 70e3),
    model = "tsl",
    n_realizations = 2,
    parameters = list(intercept = -66),
    parallel = FALSE,
    verbose = FALSE
  )

  expect_true("TSL" %in% names(sim))
  expect_equal(nrow(sim$TSL), 4)
  expect_true(all(is.finite(sim$TSL$TS)))
})

test_that("model registry guards collisions and unregisters user models", {
  restore_registry <- capture_model_registry_state()
  on.exit(restore_registry(), add = TRUE)

  acousticTS::register_model(
    name = "tsl",
    initialize = tsl_initialize,
    solver = TSL,
    slot = "TSL",
    aliases = "toy_tsl"
  )

  expect_error(
    acousticTS::register_model(
      name = "dwba",
      initialize = tsl_initialize,
      solver = TSL
    ),
    "Built-in model"
  )
  expect_error(
    acousticTS::register_model(
      name = "another_tsl",
      initialize = tsl_initialize,
      solver = TSL,
      aliases = "toy_tsl"
    ),
    "already in use"
  )

  expect_invisible(acousticTS::unregister_model("toy_tsl"))
  expect_false("tsl" %in% acousticTS::available_models()$model)
  expect_error(
    target_strength(
      object = fls_generate(
        shape = cylinder(length_body = 0.07, radius_body = 0.01, n_segments = 80),
        density_body = 1028.9,
        sound_speed_body = 1480.3
      ),
      frequency = 38e3,
      model = "tsl"
    ),
    "Unknown target strength model"
  )
})
