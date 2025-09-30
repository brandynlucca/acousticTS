#' Simulate target strength (TS) with flexible parameterization and batching
#' @inheritParams target_strength
#' @param model Model name. If multiple models are specified, the output will
#' be a list of data frames, one for each model.
#' @param n_realizations Number of realizations and output TS values.
#' @param parameters List containing the values, distributions, or generating
#' functions of parameter values that inform the TS model.
#' @param batch_by Optional. Specifies which parameters in \code{parameters} to
#' batch over. Simulations will be run for all combinations of these parameter
#' values. Default is \code{NULL}.
#' @param parallel Logical; whether to parallelize the simulations. Default is
#' \code{TRUE}.
#' @param n_cores Optional. Number of CPU cores to use for parallelization.
#' Default is \code{parallel::detectCores() - 1}.
#' @param verbose Logical; whether to print progress and status messages to the
#' console. Default is \code{TRUE}.
#'
#' @return A data frame (or list of data frames) with simulation results.
#'
#' @details
#' For example, if \code{batch_by = "length"} and \code{parameters["length"]}
#' is a vector of values, simulations will be run for each value of length
#' \code{n_realizations} times. If multiple parameters are specified in
#' \code{batch_by}, batching will occur over all combinations of their values.
#'
#' @section Parallelization:
#' This function uses \code{pbapply::pblapply()} for parallelized simulation
#' with progress bars. On Windows, parallelization uses PSOCK clusters, which
#' require all necessary objects and packages to be exported to worker
#' processes. On Unix-like systems, forking is used, which is generally simpler.
#'
#' @section Performance Issues:
#' Including too many parameters from \code{parameters}
#' within \code{batch_by} may cause significant performance issues or cause
#' \code{R} to crash. If intensive simulations are required, consider breaking
#' them into more manageable chunks
#'
#' @section DEPRECATION WARNING:
#' The `simulate_ts` function will be deprecated in future versions and will be
#' replaced by the `anneal` function in future versions.
#'
#' @importFrom pbapply pblapply
#' @export
simulate_ts <- function(object,
                        frequency,
                        model,
                        n_realizations,
                        parameters,
                        batch_by = NULL,
                        parallel = TRUE,
                        n_cores = parallel::detectCores() - 1,
                        verbose = TRUE) {
  # Future deprecation warning =================================================
  message(
    "DeprecationWarning: The `simulate_ts` function will be replaced by the ",
    "`anneal` function in future releases."
  )
  # Validate that object is of the correct class ===============================
  stopifnot(
    "'object' must be a 'scatterer'-based class" = inherits(object, "Scatterer")
  )
  # Validate model =============================================================
  unexpected_model <- model[!(model %in% names(model_registry))]
  if (length(unexpected_model) > 0) {
    stop(
      "The following user-defined models are not supported by `acousticTS`: ",
      paste(unexpected_model, collapse = ", ")
    )
  }
  # Handle batching parameter argument =========================================
  if (!is.null(batch_by)) {
    # ---- Validate ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    missing_parameters <- setdiff(batch_by, names(parameters))
    if (length(missing_parameters) > 0) {
      stop(
        "The following 'batch_by' are missing from 'parameters': ",
        paste(missing_parameters, collapse = ", ")
      )
    }
    # ---- Set up batching values ++++++++++++++++++++++++++++++++++++++++++++++
    batch_values <- lapply(batch_by, function(param) {
      value <- parameters[[param]]
      # ---- Case where parameter is a generating function +++++++++++++++++++++
      if (is.function(value)) {
        result <- value()
        if (length(result) == 0) {
          stop(
            paste0(
              "Batch parameter '", param, "' function must return at least 1 ",
              "valid value."
            )
          )
        }
        result
        # ---- Case where parameter is a defined value +++++++++++++++++++++++++
      } else {
        value
      }
    })
    # ---- Apply names +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    names(batch_values) <- batch_by
    # ---- Prepare parameter grid ++++++++++++++++++++++++++++++++++++++++++++++
    parameter_grid <- expand.grid(
      lapply(batch_values, function(x) seq_along(x)),
      stringsAsFactors = FALSE
    )
    names(parameter_grid) <- paste0(names(batch_values), "_idx")
    # ---- Create indexed dataframe ++++++++++++++++++++++++++++++++++++++++++++
    simulation_grid <- data.frame(
      realization = rep(seq_len(n_realizations), times = nrow(parameter_grid)),
      parameter_grid[rep(seq_len(nrow(parameter_grid)),
                         each = n_realizations
                         ), , drop = FALSE],
      row.names = NULL
      )
  } else {
    # ---- Non-batched case for dataframe ++++++++++++++++++++++++++++++++++++++
    simulation_grid <- data.frame(
      realization = seq_len(n_realizations),
      row.names = NULL
    )
  }
  # Simulate/map parameter values ==============================================
  parameter_names <- names(parameters)
  parameter_matrix <- as.data.frame(
    stats::setNames(
      lapply(parameter_names, function(param) {
        value <- parameters[[param]]
        if (!is.null(batch_by) && param %in% batch_by) {
          batch_values[[param]][
            simulation_grid[[paste0(param, "_idx")]]
          ]
        } else if (is.function(value)) {
          replicate(nrow(simulation_grid), value())
        } else if (length(value) == 1) {
          rep(value, nrow(simulation_grid))
        } else if (length(value) == nrow(simulation_grid)) {
          value
        } else {
          sim_type <- ifelse(is.null(batch_by),
            "realizations",
            "batched realizations"
          )
          stop(
            paste0(
              "Length of parameter '", param, "' does not match number of ",
              sim_type, " [n=", nrow(simulation_grid), "]."
            )
          )
        }
      }),
      parameter_names
    )
  )
  # ---- Bind to simulation grid +++++++++++++++++++++++++++++++++++++++++++++++
  if (length(parameters) > 0) {
    simulation_grid <- cbind(simulation_grid, parameter_matrix)
  }
  # Run simulations ============================================================
  if (verbose) {
    cat("====================================\n")
    cat("Scatterer-class:", class(object)[[1]], "\n")
    cat("Model(s):", paste(model, collapse = ", "), "\n")
    if (!is.null(batch_by)) {
      cat(
        "Batching parameter(s):",
        paste(batch_by, collapse = ", "), "\n"
      )
    } else {
      cat("")
    }
    cat(
      "Simulated parameters:",
      paste(names(parameters), collapse = ", "), "\n"
    )
    cat("Total simulation realizations:", nrow(simulation_grid), "\n")

    cat("Parallelize TS calculations:", parallel, "\n")
    cat("====================================\n")
  }
  # ---- Parallelized approach +++++++++++++++++++++++++++++++++++++++++++++++++
  if (parallel) {
    if (verbose) {
      cat("====================================\n")
      cat("Preparing parallelized simulations\n")
      cat("Number of cores:", paste0(n_cores), "\n")
    }
    # ---- Set up cluster ++++++++++++++++++++++++++++++++++++++++++++++++++++++
    cluster <- parallel::makeCluster(n_cores)
    on.exit(parallel::stopCluster(cluster))
    if (verbose) print(cluster)
    # ---- Make sure appropriate libraries are loaded ++++++++++++++++++++++++++
    parallel::clusterEvalQ(
      cluster,
      {
        loadNamespace("acousticTS")
        loadNamespace("methods")
        NULL
      }
    )
    # ---- Export the required objects/functions to each core ++++++++++++++++++
    parallel::clusterExport(
      cluster,
      c(
        # Arguments
        "object", "frequency", "model",
        # Intermediate variables
        "simulation_grid",
        # Functions
        ".discover_reforge_params", ".get_TS", "reforge", "target_strength"
      ),
      envir = environment()
    )
    # ---- Sequential approach +++++++++++++++++++++++++++++++++++++++++++++++++
  } else {
    if (verbose) {
      cat("====================================\n")
      cat("Preparing sequential simulations\n")
    }
    # ---- Set up cluster ++++++++++++++++++++++++++++++++++++++++++++++++++++++
    cluster <- NULL
  }
  # Run TS simulations =========================================================
  pbapply::pboptions(type = if (verbose) "txt" else "none")
  results_list <- pbapply::pblapply(
    seq_len(nrow(simulation_grid)),
    function(grid_index) {
      .get_TS(
        grid_index, object, parameters, simulation_grid, frequency,
        model
      )
    },
    cl = cluster
  )
  # Prepare output =============================================================
  if (length(results_list) == 0) {
    return(NULL)
  }
  # ---- Get model names +++++++++++++++++++++++++++++++++++++++++++++++++++++++
  model_names <- names(results_list[[1]])
  # ---- Concatenate into a single dataframe +++++++++++++++++++++++++++++++++++
  final_result <- stats::setNames(
    lapply(
      model_names,
      function(mod_name) {
        model_data <- lapply(results_list, function(x) x[[mod_name]])
        df <- do.call(rbind, model_data)
        rownames(df) <- NULL
        df
      }
    ),
    model_names
  )
  # Return output dataframe ====================================================
  if (verbose) {
    cat("====================================\n")
    cat("Simulations complete!\n")
    cat("====================================\n")
  }
  return(final_result)
}

#' Run a single simulation for a given parameter grid index
#'
#' This helper function extracts parameter values for a given simulation grid
#' index, updates the working scatterer object accordingly (including reforge
#' if needed), runs the target strength calculation, and formats the results
#' for output.
#'
#' @inheritParams simulate_ts
#' @param grid_index Integer index of the row in \code{simulation_grid} to
#' simulate.
#'
#' @return
#' A named list of data frames, one per model, each containing the simulation
#' results along with the model name and parameter values for this simulation.
#'
#' @details
#' This function is intended for internal use in both parallel and sequential
#' simulation workflows. It handles object re-forging, parameter assignment,
#' and result extraction for a single simulation instance.
#'
#' @seealso \code{\link{target_strength}}, \code{\link{reforge}}
#'
#' @keywords internal
#' @noRd
.get_TS <- function(grid_index,
                    object,
                    parameters,
                    simulation_grid,
                    frequency,
                    model) {
  # Extract parameter values for this grid index ===============================
  parameter_values <- as.list(simulation_grid[grid_index, , drop = FALSE])
  # Create working copy of scattering object ===================================
  working_object <- object
  # Reforge object, if parameters require ======================================
  tryCatch(
    {
      # ---- Get valid reforge parameters for this object type +++++++++++++++++
      valid_reforge_params <- .discover_reforge_params(
        class(working_object)
      )
      # ---- Filter parameters for reforge (only allowed parameters pass)
      reforge_parameters <- parameter_values[
        names(parameter_values) %in% valid_reforge_params
      ]
      # ---- Set up arguments ++++++++++++++++++++++++++++++++++++++++++++++++++
      reforge_args <- c(list(object = working_object), reforge_parameters)
      # ---- Reforge working object ++++++++++++++++++++++++++++++++++++++++++++
      working_object <- do.call(reforge, reforge_args)
    },
    error = function(e) NULL
  )
  # Override parameter definitions where appropriate ===========================
  for (param in names(parameter_values)) {
    tryCatch(
      {
        if (param %in% names(working_object@body)) {
          methods::slot(
            working_object,
            "body"
          )[[param]] <- parameter_values[[param]]
        }
      },
      error = function(e) {
        tryCatch(
          {
            if (param %in% names(working_object@bladder)) {
              methods::slot(
                working_object,
                "bladder"
              )[[param]] <- parameter_values[[param]]
            }
          },
          error = function(e) NULL
        )
      }
    )
  }
  # Calculate acousticTS =======================================================
  # [pun intended :)]
  # ---- Set up TS function arguments ++++++++++++++++++++++++++++++++++++++++++
  ts_args <- c(
    list(
      object = working_object,
      frequency = frequency,
      model = model,
      verbose = FALSE
    ),
    parameter_values
  )
  # ---- Compute TS ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  working_object <- do.call(target_strength, ts_args)
  # Extract the results ========================================================
  model_results <- extract(working_object, "model")
  # Process the results for output =============================================
  model_dfs <- lapply(
    names(model_results),
    function(mod_name) {
      # ---- Get model result and ensure it's a data frame +++++++++++++++++++++
      df <- if (is.data.frame(model_results[[mod_name]])) {
        model_results[[mod_name]]
      } else {
        data.frame(value = model_results[[mod_name]])
      }
      # ---- Add model name, parameter values, and realization number ++++++++++
      cbind(
        model = mod_name,
        as.data.frame(parameter_values, optional = TRUE),
        df
      )
    }
  )
  names(model_dfs) <- names(model_results)
  # Return output ==============================================================
  model_dfs
}
