test_that("simulate_ts function works with empty parameters", {
  library(acousticTS)

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

test_that("simulate_ts function works with single sets of parameters", {
  library(acousticTS)

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
  result_curved <- simulate_ts(
    object = krill,
    frequency = frequency,
    model = "DWBA_curved",
    n_realizations = 1,
    parameters = list(radius_curvature_ratio = 3),
    parallel = FALSE,
    verbose = FALSE
  )

  # Check values
  expect_false(all(result_curved$DWBA_curved$TS == reference_df$TS))
})

test_that("simulate_ts works with batch_by parameter", {
  library(acousticTS)

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

  # Check values
  expect_true(length(unique(result_df$TS)) == 6)
  expect_equal(result_df$TS[result_df$frequency == 38e3],
    c(-120.848794, -102.844556, -92.375204),
    tolerance = 1e-6
  )
  expect_equal(result_df$TS[result_df$frequency == 120e3],
    c(-101.045293, -83.563902, -73.985651),
    tolerance = 1e-6
  )
})

test_that("simulate_ts works with generating functions", {
  library(acousticTS)

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

test_that("simulate_ts works with multiple generating functions", {
  library(acousticTS)

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
  library(acousticTS)

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
  library(acousticTS)

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
    "not supported"
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
  library(acousticTS)

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
