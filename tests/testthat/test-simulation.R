library(acousticTS)

test_that("simulate_ts function works with empty parameters", {
  # Test with a simple CAL object
  cal_obj <- cal_generate()
  frequency <- c(38e3, 120e3)

  # Basic simulation parameters -- nothing varies
  parameters <- list()

  # Run simulation
  result <- simulate_ts(
    object = cal_obj,
    frequency = frequency,
    model = "calibration",
    n_realizations = 5,
    parameters = parameters,
    parallel = FALSE,
    verbose = FALSE
  )

  # Check that result is a data frame
  expect_type(result, "list")

  # Get the results
  result_df <- result$calibration

  # Check that we get the expected number of rows
  expect_equal(nrow(result_df), 5 * length(frequency))

  # Check that required columns exist
  expect_true("frequency" %in% colnames(result_df))
  expect_true("TS" %in% colnames(result_df))

  # Check values
  cal_reference <- cal_generate()
  cal_reference <- target_strength(
    cal_reference,
    frequency = frequency,
    model = "calibration"
  )
  reference_df <- cal_reference@model$calibration

  expect_equal(unique(result_df$TS), unique(reference_df$TS))
})

test_that("simulate_ts is the stable simulation entry point", {
  cal_obj <- cal_generate()

  expect_silent(
    simulate_ts(
      object = cal_obj,
      frequency = 38e3,
      model = "calibration",
      n_realizations = 1,
      parameters = list(),
      parallel = FALSE,
      verbose = FALSE
    )
  )
})

test_that("simulate_ts function works with single sets of parameters", {
  # Test with a simple CAL object
  data(krill)
  frequency <- c(38e3, 120e3)

  # Basic simulation parameters -- nothing varies
  parameters <- list(
    length = 20e-3,
    radius = 4e-3
  )

  # Run simulation
  result <- simulate_ts(
    object = krill,
    frequency = frequency,
    model = "DWBA",
    n_realizations = 1,
    parameters = parameters,
    parallel = FALSE,
    verbose = FALSE
  )

  # Check that result is a data frame
  expect_type(result, "list")

  # Get the results
  result_df <- result$DWBA

  # Check that we get the expected number of rows
  expect_equal(nrow(result_df), length(frequency))

  # Check that required columns exist
  expect_true("frequency" %in% colnames(result_df))
  expect_true("TS" %in% colnames(result_df))

  # Check values
  krill_reference <- target_strength(
    krill,
    frequency = frequency,
    model = "DWBA"
  )
  reference_df <- krill_reference@model$DWBA

  expect_false(all(unique(result_df$TS) == unique(reference_df$TS)))

  # Pass undefined variable (i.e. not already within the object)
  expect_warning(
    result_curved <- simulate_ts(
      object = krill,
      frequency = frequency,
      model = "DWBA_curved",
      n_realizations = 1,
      parameters = list(radius_curvature_ratio = 3),
      parallel = FALSE,
      verbose = FALSE
    ),
    "deprecated"
  )

  # Check values
  expect_false(all(result_curved$DWBA_curved$TS == reference_df$TS))
})

test_that("simulate_ts accepts model names case-insensitively", {
  data(krill)
  frequency <- c(38e3, 120e3)
  parameters <- list(length = 20e-3)

  result_upper <- simulate_ts(
    object = krill,
    frequency = frequency,
    model = "DWBA",
    n_realizations = 1,
    parameters = parameters,
    parallel = FALSE,
    verbose = FALSE
  )

  result_lower <- simulate_ts(
    object = krill,
    frequency = frequency,
    model = "dwba",
    n_realizations = 1,
    parameters = parameters,
    parallel = FALSE,
    verbose = FALSE
  )

  expect_true("DWBA" %in% names(result_upper))
  expect_true("DWBA" %in% names(result_lower))
  expect_equal(result_upper$DWBA$TS, result_lower$DWBA$TS)
})

test_that("simulate_ts works with batch_by parameter", {
  # Test batching with different parameter values
  data(krill)
  frequency <- c(38e3, 120e3)

  # Parameters with batch_by
  parameters <- list(
    length = c(10e-3, 20e-3, 30e-3) # Will batch over these values
  )

  # Run simulation with batching
  result <- simulate_ts(
    object = krill,
    frequency = frequency,
    model = "DWBA",
    n_realizations = 1,
    parameters = parameters,
    batch_by = "length",
    parallel = FALSE,
    verbose = FALSE
  )

  # Check that result is a data frame
  expect_type(result, "list")
  result_df <- result$DWBA
  expect_s3_class(result_df, "data.frame")

  # Should have more rows due to batching
  expect_true(nrow(result_df) == 6)

  # Check values against direct reforge() + target_strength() references
  reference_df <- do.call(
    rbind,
    lapply(parameters$length, function(length_value) {
      obj <- reforge(krill, length = length_value)
      obj <- target_strength(
        object = obj,
        frequency = frequency,
        model = "DWBA"
      )
      df <- extract(obj, "model")$DWBA
      data.frame(
        length = length_value,
        frequency = df$frequency,
        TS = df$TS
      )
    })
  )

  result_df <- result_df[order(result_df$length, result_df$frequency), ]
  reference_df <- reference_df[order(reference_df$length, reference_df$frequency), ]
  rownames(result_df) <- NULL
  rownames(reference_df) <- NULL

  expect_true(length(unique(result_df$TS)) == 6)
  expect_equal(result_df$length, reference_df$length)
  expect_equal(result_df$frequency, reference_df$frequency)
  expect_equal(result_df$TS, reference_df$TS, tolerance = 1e-6)
})

test_that("simulate_ts supports structured reforge parameters", {
  data(krill)

  result <- simulate_ts(
    object = krill,
    frequency = 38e3,
    model = "DWBA",
    n_realizations = 1,
    parameters = list(
      body_target = c(length = 0.02),
      n_segments_body = 120
    ),
    parallel = FALSE,
    verbose = FALSE
  )

  expect_type(result, "list")
  expect_s3_class(result$DWBA, "data.frame")
  expect_true(is.list(result$DWBA$body_target))
  expect_equal(result$DWBA$body_target[[1]], c(length = 0.02))
  expect_equal(result$DWBA$n_segments_body, 120)
  expect_true(all(is.finite(result$DWBA$TS)))
})

test_that("simulate_ts supports convenience FLS reforge aliases", {
  data(krill)

  result <- simulate_ts(
    object = krill,
    frequency = 38e3,
    model = "DWBA",
    n_realizations = 1,
    parameters = list(
      length_body = 0.02,
      n_segments_body = 120
    ),
    parallel = FALSE,
    verbose = FALSE
  )

  expect_type(result, "list")
  expect_s3_class(result$DWBA, "data.frame")
  expect_equal(result$DWBA$length_body, 0.02)
  expect_equal(result$DWBA$n_segments_body, 120)
  expect_true(all(is.finite(result$DWBA$TS)))
})

test_that("simulate_ts supports batched structured reforge parameters", {
  data(krill)

  result <- simulate_ts(
    object = krill,
    frequency = 38e3,
    model = "DWBA",
    n_realizations = 1,
    parameters = list(
      body_target = list(
        c(length = 0.02),
        c(length = 0.03)
      )
    ),
    batch_by = "body_target",
    parallel = FALSE,
    verbose = FALSE
  )

  expect_type(result, "list")
  expect_s3_class(result$DWBA, "data.frame")
  expect_equal(nrow(result$DWBA), 2)
  expect_true(is.list(result$DWBA$body_target))
  expect_equal(result$DWBA$body_target[[1]], c(length = 0.02))
  expect_equal(result$DWBA$body_target[[2]], c(length = 0.03))
  expect_true(length(unique(result$DWBA$TS)) == 2)
})

test_that("simulate_ts supports batched convenience FLS reforge aliases", {
  data(krill)
  frequency <- c(38e3, 120e3)

  result <- simulate_ts(
    object = krill,
    frequency = frequency,
    model = "DWBA",
    n_realizations = 1,
    parameters = list(
      length_body = c(0.02, 0.03)
    ),
    batch_by = "length_body",
    parallel = FALSE,
    verbose = FALSE
  )

  reference_df <- do.call(
    rbind,
    lapply(c(0.02, 0.03), function(length_value) {
      obj <- reforge(krill, body_target = c(length = length_value))
      obj <- target_strength(
        object = obj,
        frequency = frequency,
        model = "DWBA"
      )
      df <- extract(obj, "model")$DWBA
      data.frame(
        length_body = length_value,
        frequency = df$frequency,
        TS = df$TS
      )
    })
  )

  result_df <- result$DWBA[order(result$DWBA$length_body, result$DWBA$frequency), ]
  reference_df <- reference_df[order(reference_df$length_body, reference_df$frequency), ]
  rownames(result_df) <- NULL
  rownames(reference_df) <- NULL

  expect_equal(result_df$length_body, reference_df$length_body)
  expect_equal(result_df$frequency, reference_df$frequency)
  expect_equal(result_df$TS, reference_df$TS, tolerance = 1e-6)
})

test_that("simulate_ts propagates invalid structured reforge inputs", {
  data(krill)

  expect_error(
    simulate_ts(
      object = krill,
      frequency = 38e3,
      model = "DWBA",
      n_realizations = 1,
      parameters = list(body_target = c(depth = 0.02)),
      parallel = FALSE,
      verbose = FALSE
    ),
    "invalid dimensions"
  )
})

test_that("simulate_ts rejects conflicting convenience and explicit targets", {
  data(krill)

  expect_error(
    simulate_ts(
      object = krill,
      frequency = 38e3,
      model = "DWBA",
      n_realizations = 1,
      parameters = list(
        body_target = c(length = 0.02),
        length_body = 0.03
      ),
      parallel = FALSE,
      verbose = FALSE
    ),
    "Specify either 'body_target' or its convenience aliases"
  )
})

test_that("simulate_ts works with generating functions", {
  # Test with a generating function
  data(krill)
  frequency <- c(38e3)

  # Parameters with a generating function
  parameters <- list(
    length = function() stats::rnorm(1, mean = 40e-3, sd = 5e-3)
  )

  # Run simulation
  result <- simulate_ts(
    object = krill,
    frequency = frequency,
    model = "DWBA",
    n_realizations = 10,
    parameters = parameters,
    parallel = FALSE,
    verbose = FALSE
  )

  # Check that result is a data frame
  expect_type(result, "list")
  result_df <- result$DWBA
  expect_s3_class(result_df, "data.frame")

  # Get the results
  result_df <- result$DWBA

  # Check that we get the expected number of rows
  expect_equal(nrow(result_df), length(frequency) * 10)

  # Check that required columns exist
  expect_true("frequency" %in% colnames(result_df))
  expect_true("TS" %in% colnames(result_df))

  # Ensure that values are not duplicated
  expect_true(length(unique(result_df$TS)) == length(frequency) * 10)
})

test_that("simulate_ts supports convenience reforge aliases from generators", {
  data(krill)

  result <- simulate_ts(
    object = krill,
    frequency = 38e3,
    model = "DWBA",
    n_realizations = 10,
    parameters = list(
      length_body = function() stats::rnorm(1, mean = 40e-3, sd = 5e-3)
    ),
    parallel = FALSE,
    verbose = FALSE
  )

  expect_type(result, "list")
  expect_s3_class(result$DWBA, "data.frame")
  expect_equal(nrow(result$DWBA), 10)
  expect_true("length_body" %in% colnames(result$DWBA))
  expect_true(length(unique(result$DWBA$length_body)) > 1)
  expect_true(length(unique(result$DWBA$TS)) > 1)
})

test_that("simulate_ts preserves legacy curved krill workflows via length_body", {
  data(krill)

  frequency <- 120e3
  lengths <- c(0.015, 0.02)
  params_legacy <- list(
    length = lengths,
    sound_speed_sw = 1500,
    density_sw = 1026,
    g = 1.02,
    h = 1.03,
    theta = function() pi / 2,
    radius_curvature_ratio = function() 5,
    n_iterations = 5,
    n_segments_init = 15,
    length_init = 17.9e-3,
    frequency_init = 120e3,
    phase_sd_init = 0.31
  )
  params_alias <- params_legacy
  params_alias$length <- NULL
  params_alias$length_body <- lengths

  set.seed(20260331)
  legacy <- suppressWarnings(
    simulate_ts(
      object = krill,
      batch_by = "length",
      n_realizations = 2,
      frequency = frequency,
      model = "SDWBA_curved",
      parameters = params_legacy,
      parallel = FALSE,
      verbose = FALSE
    )
  )

  set.seed(20260331)
  alias <- suppressWarnings(
    simulate_ts(
      object = krill,
      batch_by = "length_body",
      n_realizations = 2,
      frequency = frequency,
      model = "SDWBA_curved",
      parameters = params_alias,
      parallel = FALSE,
      verbose = FALSE
    )
  )

  legacy_df <- legacy$SDWBA_curved[
    order(legacy$SDWBA_curved$length, legacy$SDWBA_curved$realization),
    c("length", "TS")
  ]
  alias_df <- alias$SDWBA_curved[
    order(alias$SDWBA_curved$length_body, alias$SDWBA_curved$realization),
    c("length_body", "TS")
  ]

  expect_equal(legacy_df$length, alias_df$length_body)
  expect_equal(legacy_df$TS, alias_df$TS)
})

test_that("simulate_ts works with multiple generating functions", {
  # Test with distribution parameters
  data(krill)
  frequency <- c(120e3)

  # Parameters with distribution
  parameters <- list(
    length = function() runif(1, min = 0.01, max = 0.05),
    theta = function() runif(1, min = 0, max = pi)
  )

  # Run simulation
  result <- simulate_ts(
    object = krill,
    frequency = frequency,
    model = "DWBA",
    n_realizations = 10,
    parameters = parameters,
    parallel = FALSE,
    verbose = FALSE
  )

  # Check that result is a data frame
  expect_type(result, "list")
  result_df <- result$DWBA
  expect_s3_class(result_df, "data.frame")

  # Check number of rows
  expect_equal(nrow(result_df), 10 * length(frequency))

  # Ensure that values are not duplicated
  expect_true(length(unique(result_df$TS)) == length(frequency) * 10)
})

test_that("simulate_ts handles mixed parameter types", {
  # Test mixing single values, vectors, and generating functions
  data(krill)
  frequency <- c(38e3, 70e3)

  parameters <- list(
    length = 20e-3, # Single value
    radius_body = seq(1e-3, 3e-3, length.out = 10),
    theta = function() runif(1, min = 0, max = pi)
  )

  # Run simulation with batching on vector_values
  result <- simulate_ts(
    object = krill,
    frequency = frequency,
    model = "DWBA",
    parameters = parameters,
    n_realizations = 10,
    parallel = FALSE,
    verbose = FALSE
  )

  # Check that result is a data frame
  expect_type(result, "list")
  result_df <- result$DWBA
  expect_s3_class(result_df, "data.frame")

  # Get the results
  result_df <- result$DWBA

  # Check that we get the expected number of rows
  expect_equal(nrow(result_df), length(frequency) * 10)

  # Ensure that values are not duplicated
  expect_true(length(unique(result_df$TS)) == length(frequency) * 10)
})

test_that("simulate_ts validates inputs correctly", {
  # Test mixing single values, vectors, and generating functions
  data(krill)
  frequency <- c(38e3)
  parameters <- list(
    length = 20e-3, # Single value
    radius_body = seq(1e-3, 3e-3, length.out = 10),
    theta = function() runif(1, min = 0, max = pi)
  )

  # Test invalid model
  expect_error(
    simulate_ts(krill, frequency, "invalid_model", 5, parameters),
    "Unknown target strength model 'invalid_model'"
  )

  # Test missing batch_by parameter
  expect_error(
    simulate_ts(krill, frequency, "calibration", 5, parameters,
      batch_by = "nonexistent_param"
    ),
    "missing from 'parameters'"
  )

  # Test invalid object class
  expect_error(
    simulate_ts("not_a_scatterer", frequency, "calibration", 5, parameters),
    "must be a 'scatterer'-based class"
  )
})

test_that("simulate_ts works with parallel processing", {
  # Test mixing single values, vectors, and generating functions
  data(krill)
  frequency <- c(38e3, 70e3)

  parameters <- list(
    length = 20e-3, # Single value
    theta = function() {
      set.seed(999)
      runif(1, min = 0, max = pi)
    }
  )

  # Test with parallel = TRUE (default)
  result_parallel <- simulate_ts(
    object = krill,
    frequency = frequency,
    model = "DWBA",
    n_realizations = 10,
    parameters = parameters,
    parallel = TRUE,
    n_cores = 1,
    verbose = FALSE
  )

  # Test with parallel = FALSE
  result_sequential <- simulate_ts(
    object = krill,
    frequency = frequency,
    model = "DWBA",
    n_realizations = 10,
    parameters = parameters,
    parallel = FALSE,
    verbose = FALSE
  )

  # Test batch_by parameter with generating function
  result_gen <- simulate_ts(
    object = krill,
    frequency = frequency,
    model = "DWBA",
    n_realizations = 10,
    parameters = list(length = function(x) {
      set.seed(999)
      rnorm(1, 20e-3, 1e-6)
    }),
    parallel = FALSE,
    verbose = FALSE
  )

  # Check that result is a data frame
  expect_type(result_parallel, "list")
  result_df <- result_parallel$DWBA
  expect_s3_class(result_df, "data.frame")

  # Extract sequential
  sequential_df <- result_sequential$DWBA

  # Check that we get the expected number of rows
  expect_equal(nrow(result_df), length(frequency) * 10)

  # Ensure that values are not duplicated
  expect_true(length(unique(result_df$TS)) == 2)
  expect_true(all(result_df$TS == sequential_df$TS))

  # Extract generative
  gen_df <- result_gen$DWBA

  # Check that we get the expected number of rows
  expect_equal(nrow(gen_df), length(frequency) * 10)

  # Ensure that values are not duplicated
  expect_true(length(unique(gen_df$TS)) == 2)
  expect_true(all(gen_df$TS != sequential_df$TS))
})

test_that("simulate_ts works with PSOCK clusters when n_cores > 1", {
  skip_on_cran()

  data(krill)
  frequency <- c(38e3, 70e3)

  parameters <- list(
    length = 20e-3,
    theta = function() {
      set.seed(999)
      runif(1, min = 0, max = pi)
    }
  )

  result_parallel <- simulate_ts(
    object = krill,
    frequency = frequency,
    model = "dwba",
    n_realizations = 4,
    parameters = parameters,
    parallel = TRUE,
    n_cores = 2,
    verbose = FALSE
  )

  result_sequential <- simulate_ts(
    object = krill,
    frequency = frequency,
    model = "DWBA",
    n_realizations = 4,
    parameters = parameters,
    parallel = FALSE,
    verbose = FALSE
  )

  expect_true("DWBA" %in% names(result_parallel))
  expect_true("DWBA" %in% names(result_sequential))
  expect_equal(result_parallel$DWBA$theta, result_sequential$DWBA$theta)
  expect_equal(result_parallel$DWBA$length, result_sequential$DWBA$length)
  expect_equal(result_parallel$DWBA$frequency, result_sequential$DWBA$frequency)
  if (dir.exists(file.path(getNamespaceInfo(asNamespace("acousticTS"), "path"), "Meta"))) {
    expect_equal(result_parallel$DWBA$TS, result_sequential$DWBA$TS)
  } else {
    expect_true(all(is.finite(result_parallel$DWBA$TS)))
    expect_true(all(is.finite(result_sequential$DWBA$TS)))
  }
})

test_that("simulate_ts supports theta_body for FLS objects in PSOCK mode", {
  skip_on_cran()

  obj <- fls_generate(
    shape = cylinder(
      length_body = 0.05,
      radius_body = 0.003,
      n_segments = 80
    ),
    density_body = 1045,
    sound_speed_body = 1520
  )

  frequency <- seq(38e3, 50e3, by = 2e3)
  parameters <- list(
    theta_body = function() {
      set.seed(999)
      runif(1, min = 0.5 * pi, max = pi)
    }
  )

  result_parallel <- simulate_ts(
    object = obj,
    frequency = frequency,
    model = "DWBA",
    n_realizations = 4,
    parameters = parameters,
    parallel = TRUE,
    n_cores = 2,
    verbose = FALSE
  )

  result_sequential <- simulate_ts(
    object = obj,
    frequency = frequency,
    model = "dwba",
    n_realizations = 4,
    parameters = parameters,
    parallel = FALSE,
    verbose = FALSE
  )

  expect_true("DWBA" %in% names(result_parallel))
  expect_true("DWBA" %in% names(result_sequential))
  expect_equal(result_parallel$DWBA$theta_body, result_sequential$DWBA$theta_body)
  expect_equal(result_parallel$DWBA$frequency, result_sequential$DWBA$frequency)
  if (dir.exists(file.path(getNamespaceInfo(asNamespace("acousticTS"), "path"), "Meta"))) {
    expect_equal(result_parallel$DWBA$TS, result_sequential$DWBA$TS)
  } else {
    expect_true(all(is.finite(result_parallel$DWBA$TS)))
    expect_true(all(is.finite(result_sequential$DWBA$TS)))
  }
})

test_that("Simulation errors are raised as expected", {
  # Test mixing single values, vectors, and generating functions
  data(krill)
  frequency <- c(38e3, 70e3)
  parameters <- list(
    length = function(x) {
      c()
    } # NULL value
  )

  # Test case with invalid output [empty] from generating function
  expect_error(
    simulate_ts(
      object = krill,
      frequency = frequency,
      model = "DWBA",
      batch_by = "length",
      n_realizations = 10,
      parameters = parameters,
      parallel = FALSE,
      verbose = FALSE
    ),
    "Batch parameter 'length' function must return at least 1 valid value."
  )
})

test_that("simulation helper utilities print headers and prepare optional clusters", {
  skip_on_cran()

  cal_obj <- cal_generate()
  simulation_grid <- data.frame(realization = 1:2)

  header <- capture.output(
    acousticTS:::.print_simulation_header(
      object = cal_obj,
      model = "calibration",
      batch_by = NULL,
      parameters = list(length = 0.02),
      parallel = FALSE,
      simulation_grid = simulation_grid
    )
  )
  expect_true(any(grepl("Scatterer-class: CAL", header, fixed = TRUE)))
  expect_true(any(grepl("Total simulation realizations: 2", header, fixed = TRUE)))

  sequential <- capture.output(
    cluster <- acousticTS:::.prepare_simulation_cluster(
      parallel = FALSE,
      n_cores = 2,
      object = cal_obj,
      frequency = 38000,
      normalized_model = "calibration",
      simulation_grid = simulation_grid,
      verbose = TRUE
    )
  )
  expect_null(cluster)
  expect_true(any(grepl("Preparing sequential simulations", sequential, fixed = TRUE)))

  parallel_out <- capture.output(
    cluster <- acousticTS:::.prepare_simulation_cluster(
      parallel = TRUE,
      n_cores = 2,
      object = cal_obj,
      frequency = 38000,
      normalized_model = "calibration",
      simulation_grid = simulation_grid,
      verbose = TRUE
    )
  )
  on.exit(parallel::stopCluster(cluster), add = TRUE)

  expect_s3_class(cluster, "cluster")
  expect_true(any(grepl("Preparing parallelized simulations", parallel_out, fixed = TRUE)))
})

test_that("simulation helper utilities are exercised under coverage runs", {
  cal_obj <- cal_generate()
  simulation_grid <- data.frame(realization = 1:2)

  header <- capture.output(
    acousticTS:::.print_simulation_header(
      object = cal_obj,
      model = "calibration",
      batch_by = "length",
      parameters = list(length = 0.02),
      parallel = FALSE,
      simulation_grid = simulation_grid
    )
  )
  expect_true(any(grepl("Batching parameter(s): length", header, fixed = TRUE)))

  sequential <- capture.output(
    cluster <- acousticTS:::.prepare_simulation_cluster(
      parallel = FALSE,
      n_cores = 2,
      object = cal_obj,
      frequency = 38000,
      normalized_model = "calibration",
      simulation_grid = simulation_grid,
      verbose = TRUE
    )
  )
  expect_null(cluster)
  expect_true(any(grepl("Preparing sequential simulations", sequential, fixed = TRUE)))

  parallel_out <- capture.output(
    cluster <- acousticTS:::.prepare_simulation_cluster(
      parallel = TRUE,
      n_cores = 2,
      object = cal_obj,
      frequency = 38000,
      normalized_model = "calibration",
      simulation_grid = simulation_grid,
      verbose = TRUE
    )
  )
  on.exit(parallel::stopCluster(cluster), add = TRUE)

  expect_s3_class(cluster, "cluster")
  expect_true(any(grepl("Preparing parallelized simulations", parallel_out, fixed = TRUE)))
})

test_that("simulation helpers cover scalar batch values, empty combines, and verbose sequential execution", {
  expect_equal(
    acousticTS:::.prepare_simulation_batch_values(
      batch_by = c("length", "density_body"),
      parameters = list(length = c(0.02, 0.03), density_body = 1028.9)
    ),
    list(length = c(0.02, 0.03), density_body = 1028.9)
  )
  expect_null(acousticTS:::.combine_simulation_results(list()))

  cal_obj <- cal_generate()
  verbose_run <- capture.output(
    result <- simulate_ts(
      object = cal_obj,
      frequency = 38e3,
      model = "calibration",
      n_realizations = 1,
      parameters = list(),
      parallel = FALSE,
      verbose = TRUE
    )
  )
  expect_true(is.list(result))
  expect_s3_class(result$calibration, "data.frame")
  expect_true(any(grepl("Scatterer-class: CAL", verbose_run, fixed = TRUE)))
  expect_true(any(grepl("Simulations complete!", verbose_run, fixed = TRUE)))
})
