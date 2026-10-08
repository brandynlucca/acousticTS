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

test_that("persistent registrations round trip in an isolated registry", {
  restore <- capture_model_registry_state()
  on.exit(restore(), add = TRUE)
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  local_mocked_bindings(
    .model_registry_user_path = function() path,
    .package = "acousticTS"
  )
  reset_model_registry(remove_persisted = TRUE)
  register_model(
    "persisted_dwba",
    initialize = "acousticTS:::dwba_initialize",
    solver = "acousticTS:::DWBA",
    slot = "DWBA",
    aliases = "persisted_alias",
    persist = TRUE
  )
  expect_true(file.exists(path))
  stored <- readRDS(path)
  expect_equal(stored$persisted_dwba$canonical, "persisted_dwba")
  state <- get(".model_registry_state", asNamespace("acousticTS"))
  state$user <- list()
  state$loaded <- FALSE
  expect_true("persisted_dwba" %in% available_models()$model)
  object <- fls_generate(
    shape = cylinder(0.03, 0.003),
    g_body = 1.03,
    h_body = 1.02
  )
  actual <- target_strength(object, 38000, model = "persisted_alias")
  expected <- target_strength(object, 38000, model = "dwba")
  expect_equal(actual@model$DWBA, expected@model$DWBA)
  unregister_model("persisted_dwba", remove_persisted = TRUE)
  expect_false(file.exists(path))
  register_model(
    "persisted_dwba",
    initialize = "acousticTS:::dwba_initialize",
    solver = "acousticTS:::DWBA",
    persist = TRUE
  )
  reset_model_registry(remove_persisted = TRUE)
  expect_false(file.exists(path))
  expect_false("persisted_dwba" %in% available_models()$model)

  writeLines("invalid registry", path)
  state$loaded <- FALSE
  expect_warning(available_models(), "Could not read persisted")
  saveRDS(
    list(list(
      canonical = "broken",
      slot = "BROKEN",
      initialize_ref = "missing_function",
      solver_ref = "base::identity"
    )),
    path
  )
  state$loaded <- FALSE
  expect_warning(available_models(), "Could not resolve")
})

test_that("registry validates identifiers and resolves qualified callables", {
  restore <- capture_model_registry_state()
  on.exit(restore(), add = TRUE)
  expect_error(register_model(1, identity, identity), "must be character")
  expect_error(register_model("", identity, identity), "non-empty model names")
  expect_error(
    register_model("test_model", identity, identity, slot = ""),
    "non-empty string"
  )
  expect_error(register_model("test_model", 1, identity), "function reference")
  expect_error(
    register_model("test_model", "missing_function", identity),
    "Could not resolve"
  )
  expect_error(
    register_model("test_model", identity, identity, persist = TRUE),
    "package-qualified"
  )
  register_model("test_model", identity, identity)
  expect_error(
    register_model("test_model", identity, identity),
    "already registered"
  )
  register_model(
    "test_model",
    "base::identity",
    "base::identity",
    overwrite = TRUE
  )
  expect_identical(
    acousticTS:::.resolve_model_function_reference("base::identity"),
    identity
  )
  expect_error(
    acousticTS:::.resolve_model_function_reference(NA_character_),
    "non-empty function reference"
  )
  expect_error(unregister_model("dwba"), "Only user-registered")
  expect_equal(acousticTS:::.default_model_slot("CALIBRATION"), "calibration")
  solvers <- acousticTS:::.get_models()
  expect_identical(solvers$DWBA, acousticTS:::DWBA)
  expect_identical(
    solvers$calibration,
    acousticTS:::.resolve_model_function_reference(
      acousticTS:::.resolve_model_registry_entry("calibration")$solver
    )
  )
})

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

test_that(
  "available_models lists built-ins and target_strength resolves aliases",
  {
    models <- acousticTS::available_models()

    expect_true("calibration" %in% models$model)
    expect_true(any(models$model == "calibration" & grepl(
      "soems",
      models$aliases
    )))
    expect_false("espsms" %in% models$model)
    expect_false("epsms" %in% models$model)

    cal_obj <- target_strength(
      cal_generate(),
      frequency = 38e3,
      model = "soems"
    )

    expect_true("calibration" %in% names(cal_obj@model))
    expect_true(all(is.finite(cal_obj@model$calibration$TS)))
    expect_error(
      target_strength(cal_generate(), frequency = 38e3, model = "epsms"),
      "Unknown target strength model 'epsms'"
    )
  }
)

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
        shape = cylinder(
          length_body = 0.07, radius_body = 0.01, n_segments =
            80
        ),
        density_body = 1028.9,
        sound_speed_body = 1480.3
      ),
      frequency = 38e3,
      model = "tsl"
    ),
    "Unknown target strength model"
  )
})
