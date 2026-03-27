################################################################################
# SIMULATION HELPERS
################################################################################
# Resolve the default core count used by `simulate_ts()`.
#' @noRd
.default_simulation_cores <- function() {
  # Cap the default at two cores to stay within CRAN's shared-check guidance ===
  detected <- parallel::detectCores()
  detected <- if (length(detected) && is.finite(detected)) detected else 2L

  min(2L, max(1L, as.integer(detected) - 1L))
}

# Resolve whether `simulate_ts()` should actually use a PSOCK cluster.
#' @noRd
.resolve_simulation_parallel <- function(parallel, n_cores) {
  # Require a valid multi-core request before enabling PSOCK execution =========
  isTRUE(parallel) &&
    is.numeric(n_cores) &&
    length(n_cores) == 1 &&
    !is.na(n_cores) &&
    n_cores > 1
}

# Validate the requested target-strength models for one simulation run.
#' @noRd
.validate_simulation_models <- function(model) {
  # Normalize model aliases and confirm each one is supported ==================
  normalized_model <- .normalize_simulation_models(model)
  unexpected_model <- model[!(normalized_model %in%
    .normalize_simulation_models(names(.get_models())))]
  if (length(unexpected_model) > 0) {
    stop(
      "The following user-defined models are not supported by acousticTS: ",
      paste(unexpected_model, collapse = ", ")
    )
  }

  normalized_model
}

# Build the simulation grid, including any batched parameter expansion.
#' @noRd
.prepare_simulation_grid <- function(n_realizations, parameters, batch_by) {
  # Fall back to the simple realization index when batching is disabled ========
  if (is.null(batch_by)) {
    return(list(
      simulation_grid = data.frame(
        realization = seq_len(n_realizations),
        row.names = NULL
      ),
      batch_values = NULL
    ))
  }

  # Validate that each batching parameter is present in `parameters` ===========
  missing_parameters <- setdiff(batch_by, names(parameters))
  if (length(missing_parameters) > 0) {
    stop(
      "The following 'batch_by' are missing from 'parameters': ",
      paste(missing_parameters, collapse = ", ")
    )
  }

  # Resolve the concrete batch-value vectors before expanding the grid =========
  batch_values <- .prepare_simulation_batch_values(batch_by, parameters)
  parameter_grid <- expand.grid(
    lapply(batch_values, function(x) seq_along(x)),
    stringsAsFactors = FALSE
  )
  names(parameter_grid) <- paste0(names(batch_values), "_idx")

  list(
    simulation_grid = data.frame(
      realization = rep(seq_len(n_realizations), times = nrow(parameter_grid)),
      parameter_grid[rep(seq_len(nrow(parameter_grid)), each = n_realizations), ,
        drop = FALSE
      ],
      row.names = NULL
    ),
    batch_values = batch_values
  )
}

# Resolve the concrete vectors used for each batched simulation parameter.
#' @noRd
.prepare_simulation_batch_values <- function(batch_by, parameters) {
  # Evaluate generating functions once per batch dimension =====================
  batch_values <- lapply(batch_by, function(param) {
    value <- parameters[[param]]
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
      return(result)
    }

    value
  })
  names(batch_values) <- batch_by

  batch_values
}

# Expand the user-supplied parameter definitions across the simulation grid.
#' @noRd
.simulation_parameter_matrix <- function(parameters,
                                         batch_by,
                                         batch_values,
                                         simulation_grid) {
  # Resolve one parameter vector per simulation input ==========================
  parameter_names <- names(parameters)
  as.data.frame(
    stats::setNames(
      lapply(parameter_names, function(param) {
        .resolve_param_value(
          param_name = param,
          param_value = parameters[[param]],
          batch_by = batch_by,
          batch_values = batch_values,
          grid_size = nrow(simulation_grid),
          simulation_grid = simulation_grid
        )
      }),
      parameter_names
    )
  )
}

# Print the high-level simulation summary shown before the TS runs begin.
#' @noRd
.print_simulation_header <- function(object,
                                     model,
                                     batch_by,
                                     parameters,
                                     parallel,
                                     simulation_grid) {
  # Print the standardized simulation summary =================================
  cat("====================================\n")
  cat("Scatterer-class:", class(object)[[1]], "\n")
  cat("Model(s):", paste(model, collapse = ", "), "\n")
  if (!is.null(batch_by)) {
    cat("Batching parameter(s):", paste(batch_by, collapse = ", "), "\n")
  } else {
    cat("")
  }
  cat("Simulated parameters:", paste(names(parameters), collapse = ", "), "\n")
  cat("Total simulation realizations:", nrow(simulation_grid), "\n")
  cat("Parallelize TS calculations:", parallel, "\n")
  cat("====================================\n")
}

# Prepare the optional PSOCK cluster used by `simulate_ts()`.
#' @noRd
.prepare_simulation_cluster <- function(parallel,
                                        n_cores,
                                        object,
                                        frequency,
                                        normalized_model,
                                        simulation_grid,
                                        verbose) {
  # Fall back to sequential execution when no cluster is requested =============
  if (!parallel) {
    if (verbose) {
      cat("====================================\n")
      cat("Preparing sequential simulations\n")
    }
    return(NULL)
  }

  # Build and prime the PSOCK cluster for worker-side TS evaluation ============
  if (verbose) {
    cat("====================================\n")
    cat("Preparing parallelized simulations\n")
    cat("Number of cores:", paste0(n_cores), "\n")
  }

  cluster <- parallel::makeCluster(n_cores)
  if (verbose) print(cluster)
  parallel::clusterCall(
    cluster,
    function() {
      loadNamespace("acousticTS")
      loadNamespace("methods")
      NULL
    }
  )
  parallel::clusterExport(
    cluster,
    c("object", "frequency", "normalized_model", "simulation_grid"),
    envir = environment()
  )
  parallel::clusterExport(
    cluster,
    c(
      ".discover_reforge_params", ".get_TS", "reforge", "target_strength",
      "extract"
    ),
    envir = asNamespace("acousticTS")
  )

  cluster
}

# Bind one list of per-model simulation results into the final return object.
#' @noRd
.combine_simulation_results <- function(results_list) {
  # Return NULL for the degenerate empty-simulation case =======================
  if (length(results_list) == 0) {
    return(NULL)
  }

  # Bind each model-specific result stack into one data frame ==================
  model_names <- names(results_list[[1]])
  stats::setNames(
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
}

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
#' Default is the smaller of 2 cores and \code{parallel::detectCores() - 1}.
#' @param verbose Logical; whether to print progress and status messages to the
#' console. Default is \code{TRUE}.
#'
#' @return
#' A data frame when a single model is requested, or a named list of data
#' frames when multiple models are requested. Each returned data frame contains
#' the realized parameter values together with the modeled acoustic output for
#' each simulated run.
#'
#' @details
#' `simulate_ts()` is a workflow wrapper around repeated `target_strength()`
#' calls. It supports three broad parameter modes inside `parameters`:
#'
#' \itemize{
#'   \item scalars that are recycled across every realization,
#'   \item explicit vectors that are either aligned with the full simulation
#'   grid or with one or more batched dimensions, and
#'   \item generating functions that are re-evaluated for each realization.
#' }
#'
#' If \code{batch_by = "length"} and \code{parameters[["length"]]} is a vector
#' of candidate values, then simulations are run for each length value,
#' repeated \code{n_realizations} times. When multiple parameters are supplied
#' through \code{batch_by}, the function builds the full Cartesian grid of
#' those parameter values and runs the requested number of realizations inside
#' each batch cell.
#'
#' Parameter names are interpreted in the same way they would be if supplied
#' directly to \code{target_strength()} or to the relevant object constructor /
#' \code{reforge()} path. This means \code{simulate_ts()} can be used for:
#'
#' \itemize{
#'   \item orientation perturbations,
#'   \item material-property perturbations,
#'   \item morphology studies that trigger shape rebuilding or reforging, and
#'   \item side-by-side comparisons across one or more model families.
#' }
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
#' @examples
#' shape_obj <- cylinder(
#'   length_body = 0.05,
#'   radius_body = 0.003,
#'   n_segments = 40
#' )
#'
#' obj <- fls_generate(
#'   shape = shape_obj,
#'   density_body = 1045,
#'   sound_speed_body = 1520
#' )
#'
#' res <- simulate_ts(
#'   object = obj,
#'   frequency = seq(38e3, 50e3, by = 6e3),
#'   model = "dwba",
#'   n_realizations = 2,
#'   parameters = list(
#'     theta_body = function() runif(1, min = 0.5 * pi, max = pi),
#'     density_body = 1045
#'   ),
#'   parallel = FALSE,
#'   verbose = FALSE
#' )
#'
#' head(res)
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
                        n_cores = .default_simulation_cores(),
                        verbose = TRUE) {
  # Validate that object is of the correct class ===============================
  stopifnot(
    "'object' must be a 'scatterer'-based class" = inherits(object, "Scatterer")
  )
  parallel <- .resolve_simulation_parallel(parallel, n_cores)
  normalized_model <- .validate_simulation_models(model)
  simulation_setup <- .prepare_simulation_grid(
    n_realizations = n_realizations,
    parameters = parameters,
    batch_by = batch_by
  )
  simulation_grid <- simulation_setup$simulation_grid
  batch_values <- simulation_setup$batch_values
  # Simulate/map parameter values ==============================================
  parameter_matrix <- .simulation_parameter_matrix(
    parameters = parameters,
    batch_by = batch_by,
    batch_values = batch_values,
    simulation_grid = simulation_grid
  )
  # ---- Bind to simulation grid +++++++++++++++++++++++++++++++++++++++++++++++
  if (length(parameters) > 0) {
    simulation_grid <- cbind(simulation_grid, parameter_matrix)
  }
  # Run simulations ============================================================
  if (verbose) {
    .print_simulation_header(
      object = object,
      model = model,
      batch_by = batch_by,
      parameters = parameters,
      parallel = parallel,
      simulation_grid = simulation_grid
    )
  }
  # Prepare the optional PSOCK cluster =========================================
  cluster <- .prepare_simulation_cluster(
    parallel = parallel,
    n_cores = n_cores,
    object = object,
    frequency = frequency,
    normalized_model = normalized_model,
    simulation_grid = simulation_grid,
    verbose = verbose
  )
  if (!is.null(cluster)) {
    on.exit(parallel::stopCluster(cluster))
  }
  # Run TS simulations =========================================================
  pbapply::pboptions(type = if (verbose) "txt" else "none")
  results_list <- pbapply::pblapply(
    seq_len(nrow(simulation_grid)),
    function(grid_index) {
      .get_TS(
        grid_index, object, parameters, simulation_grid, frequency,
        normalized_model
      )
    },
    cl = cluster
  )
  # Prepare output =============================================================
  final_result <- .combine_simulation_results(results_list)
  # Return output dataframe ====================================================
  if (verbose) {
    cat("====================================\n")
    cat("Simulations complete!\n")
    cat("====================================\n")
  }
  final_result
}

#' Normalize simulation model names to the strings accepted by target_strength()
#' @param model Character vector of model names.
#' @keywords internal
#' @noRd
.normalize_simulation_models <- function(model) {
  normalized_model <- tolower(model)
  normalized_model[normalized_model == "soems"] <- "calibration"
  normalized_model
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
