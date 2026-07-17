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

  normalized_model
}

# Identify a deterministic multi-level parameter that defines a grid axis.
#' @noRd
.is_auto_batch_parameter <- function(value) {
  # Generating functions are stochastic draws, not deterministic grid axes =====
  if (is.function(value)) {
    return(FALSE)
  }
  # Structured values are a single unit unless supplied as a list of levels ====
  if (.is_structured_simulation_value(value)) {
    return(is.list(value) && length(value) > 1)
  }
  # Plain atomic inputs define a grid axis only with more than one level =======
  length(value) > 1
}

# Resolve the deterministic grid axes and the per-cell repeat count.
#' @noRd
.resolve_simulation_realizations <- function(parameters, batch_by, n_realizations) {
  # Resolve the per-cell repeat count (defaults to a single evaluation) =========
  # `n_realizations` is purely a repeat/redraw multiplier: each deterministic
  # grid cell is evaluated this many times, redrawing any generating functions.
  if (is.null(n_realizations)) {
    n_realizations <- 1L
  } else if (!is.numeric(n_realizations) || length(n_realizations) != 1 ||
    is.na(n_realizations) || n_realizations < 1) {
    stop("'n_realizations' must be a single positive integer.", call. = FALSE)
  } else {
    n_realizations <- as.integer(n_realizations)
  }

  # Every deterministic multi-level parameter defines a grid axis ==============
  # This holds regardless of `n_realizations`, so a bare multi-value vector
  # always means "sweep these values" rather than "align one per realization".
  auto_batch <- names(parameters)[
    vapply(parameters, .is_auto_batch_parameter, logical(1))
  ]
  batch_by <- union(batch_by %||% character(0), auto_batch)
  if (length(batch_by) == 0) {
    batch_by <- NULL
  }
  list(n_realizations = n_realizations, batch_by = batch_by)
}

# Build the simulation grid, including any batched parameter expansion.
#' @noRd
.prepare_simulation_grid <- function(n_realizations, parameters, batch_by,
                                     permute = TRUE) {
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
  parameter_grid <- .simulation_index_grid(batch_values, permute)
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

# Build the per-axis index grid, crossing or pairing the varied parameters.
#' @noRd
.simulation_index_grid <- function(batch_values, permute) {
  # Cross every axis into the full Cartesian grid ==============================
  if (isTRUE(permute)) {
    return(expand.grid(
      lapply(batch_values, function(x) seq_along(x)),
      stringsAsFactors = FALSE
    ))
  }

  # Otherwise zip the axes element-wise into aligned combinations =============
  axis_lengths <- vapply(batch_values, length, integer(1))
  n_paired <- max(axis_lengths)
  mismatched <- axis_lengths != n_paired & axis_lengths != 1L
  if (any(mismatched)) {
    stop(
      "Paired simulation ('permute = FALSE') requires every varied parameter ",
      "to share the same number of values. Mismatched: ",
      paste(
        sprintf(
          "'%s' [%d]",
          names(batch_values)[mismatched],
          axis_lengths[mismatched]
        ),
        collapse = ", "
      ),
      ".",
      call. = FALSE
    )
  }

  # Length-one axes stay constant; the rest advance together through the pairs =
  as.data.frame(
    lapply(axis_lengths, function(n) {
      if (n == 1L) rep(1L, n_paired) else seq_len(n_paired)
    }),
    stringsAsFactors = FALSE
  )
}

# Resolve the concrete vectors used for each batched simulation parameter.
#' @noRd
.prepare_simulation_batch_values <- function(batch_by, parameters) {
  # Resolve one normalized candidate set per batched parameter ================
  batch_values <- lapply(batch_by, function(param) {
    .normalize_simulation_batch_values(
      param_name = param,
      param_value = parameters[[param]]
    )
  })
  names(batch_values) <- batch_by

  # Return the normalized batch-value definitions ==============================
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
  resolved_values <- stats::setNames(
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

  # Preserve structured parameters as list-columns in the simulation grid =====
  out <- data.frame(row.names = seq_len(nrow(simulation_grid)))
  for (param in parameter_names) {
    out[[param]] <- resolved_values[[param]]
  }

  # Return the resolved simulation parameter matrix ===========================
  out
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

# Resolve the package-library path used by PSOCK workers.
#' @noRd
.resolve_simulation_worker_library <- function(package = "acousticTS") {
  # Resolve the active namespace path used by the current R session ============
  package_dir <- normalizePath(
    getNamespaceInfo(asNamespace(package), "path"),
    winslash = "/",
    mustWork = FALSE
  )
  if (dir.exists(package_dir) &&
    file.exists(file.path(package_dir, "DESCRIPTION"))) {
    # Reuse the active library path when the current namespace is installed ====
    if (dir.exists(file.path(package_dir, "Meta"))) {
      return(dirname(package_dir))
    }
  } else {
    package_dir <- NULL
  }

  # Reuse another installed package when no source checkout is active =========
  if (is.null(package_dir)) {
    installed_paths <- file.path(.libPaths(), package)
    installed_paths <- installed_paths[
      file.exists(file.path(installed_paths, "DESCRIPTION"))
    ]
    if (length(installed_paths) > 0) {
      return(dirname(normalizePath(installed_paths[[1]], winslash = "/")))
    }
    return(NULL)
  }

  # Reuse a session-cached worker library when it matches this source tree ====
  cache_key <- paste0(package, ".simulation_worker_library")
  cache <- getOption(cache_key, NULL)
  cache_lib <- if (is.list(cache)) cache$library else NULL
  cache_src <- if (is.list(cache)) cache$source else NULL
  cache_pkg <- if (!is.null(cache_lib)) file.path(cache_lib, package) else NULL
  if (!is.null(cache_lib) &&
    identical(cache_src, package_dir) &&
    dir.exists(cache_pkg)) {
    return(cache_lib)
  }

  # Install the current source checkout into a temporary worker library =======
  worker_lib <- file.path(tempdir(), paste0(package, "-simulation-lib"))
  if (dir.exists(file.path(worker_lib, package))) {
    unlink(file.path(worker_lib, package), recursive = TRUE, force = TRUE)
  }
  dir.create(worker_lib, recursive = TRUE, showWarnings = FALSE)

  r_executable <- file.path(
    R.home("bin"),
    if (.Platform$OS.type == "windows") "Rcmd.exe" else "R"
  )
  install_args <- if (.Platform$OS.type == "windows") {
    c(
      "INSTALL",
      "-l", worker_lib,
      "--no-help",
      "--no-html",
      "--no-demo",
      package_dir
    )
  } else {
    c(
      "CMD", "INSTALL",
      "-l", worker_lib,
      "--no-help",
      "--no-html",
      "--no-demo",
      package_dir
    )
  }
  install_output <- system2(
    r_executable,
    args = install_args,
    stdout = TRUE,
    stderr = TRUE
  )
  install_status <- attr(install_output, "status") %||% 0L
  if (!identical(as.integer(install_status), 0L)) {
    stop(
      "Unable to prepare the temporary worker installation for parallel ",
      "simulate_ts() runs.\n",
      paste(install_output, collapse = "\n"),
      call. = FALSE
    )
  }

  # Cache and return the worker library used for PSOCK execution ==============
  options(
    stats::setNames(
      list(list(library = worker_lib, source = package_dir)),
      cache_key
    )
  )
  worker_lib
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

  # Resolve the library location that workers should load acousticTS from =====
  worker_lib <- .resolve_simulation_worker_library("acousticTS")
  cluster <- parallel::makeCluster(n_cores)
  if (verbose) print(cluster)
  parallel::clusterCall(
    cluster,
    function(worker_lib) {
      if (!is.null(worker_lib)) {
        .libPaths(c(worker_lib, .libPaths()))
      }
      loadNamespace("acousticTS")
      loadNamespace("methods")
      NULL
    },
    worker_lib = worker_lib
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
#' @param n_realizations Optional number of repeated evaluations of each
#' deterministic grid cell (default \code{1}). It is a pure repeat/redraw
#' multiplier: every deterministic multi-value parameter is always treated as a
#' grid axis, and \code{n_realizations} controls how many times each resulting
#' cell is evaluated, redrawing any generating functions each time. The total
#' number of realizations is therefore the size of the Cartesian grid multiplied
#' by \code{n_realizations} (for example, one length and two radii yield two
#' cells; with \code{n_realizations = 5} that is ten realizations). Use it to
#' draw repeated samples from generating functions (distributions).
#' @param parameters List containing the values, distributions, or generating
#' functions of parameter values that inform the TS model. Defaults to an empty
#' list, which runs a single realization of the unmodified object.
#' @param batch_by Optional. Specifies which parameters in \code{parameters} to
#' batch over. Simulations will be run for all combinations of these parameter
#' values. Default is \code{NULL}.
#' @param permute Logical; controls how varied parameters combine. When
#' \code{TRUE} (default) every deterministic multi-value parameter is crossed
#' into the full Cartesian grid. When \code{FALSE} the varied parameters are
#' instead paired (zipped) element-wise, so they advance together rather than
#' combinatorially; in that case every varied parameter must supply the same
#' number of values (length-one parameters are held constant).
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
#'   \item generating functions that are re-evaluated for each realization, and
#'   \item structured values such as named target-dimension vectors used by
#'   \code{reforge()} (for example \code{body_target = c(length = 0.03)}).
#' }
#'
#' If \code{batch_by = "length"} and \code{parameters[["length"]]} is a vector
#' of candidate values, then simulations are run for each length value,
#' repeated \code{n_realizations} times. When multiple parameters are supplied
#' through \code{batch_by}, the function builds the full Cartesian grid of
#' those parameter values and runs the requested number of realizations inside
#' each batch cell.
#'
#' Batching is automatic and follows one consistent rule: every deterministic
#' multi-value parameter is a grid axis, whether or not \code{n_realizations} or
#' \code{batch_by} are supplied. A bare vector therefore always means "sweep
#' these values" rather than "align one value per realization". \code{batch_by}
#' remains available to document intent, but naming an axis is no longer
#' required for it to be swept.
#'
#' To vary parameters \emph{together} (paired rather than crossed) set
#' \code{permute = FALSE}. The varied parameters are then zipped element-wise
#' and must share the same number of values. For example,
#' \code{parameters = list(theta_body = c(1, 2), density_body = c(1040, 1050))}
#' yields four runs by default but two paired runs -
#' \code{(1, 1040)} and \code{(2, 1050)} - under \code{permute = FALSE}.
#' Alternatively, values that belong to one \code{reforge()} target can be
#' paired within a single structured parameter, for example
#' \code{parameters = list(body_target = list(c(length = 0.02, radius = 0.002),
#' c(length = 0.03, radius = 0.003)))}.
#'
#' Structured batch values should be wrapped in a list so that each candidate is
#' preserved as one unit. For example, use
#' \code{parameters = list(body_target = list(c(length = 0.02), c(length = 0.03))))}
#' when batching across multiple explicit `reforge()` targets.
#'
#' Convenience dimension aliases are also supported for compatible
#' \code{reforge()} methods. For example, \code{length_body = 0.03} is treated
#' the same as \code{body_target = c(length = 0.03)} for fluid-like scatterers,
#' while retaining the original \code{length_body} column in the returned
#' simulation output.
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
#' with progress bars. The current implementation uses PSOCK clusters for
#' worker execution across platforms, including Windows, macOS, and Linux.
#' That means worker processes need access to the package namespace and any
#' required exported objects, and it also means startup overhead is more
#' noticeable for very small simulation jobs than for larger batched runs.
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
                        n_realizations = NULL,
                        parameters = list(),
                        batch_by = NULL,
                        permute = TRUE,
                        parallel = TRUE,
                        n_cores = .default_simulation_cores(),
                        verbose = TRUE) {
  # Validate that object is of the correct class ===============================
  stopifnot(
    "'object' must be a 'scatterer'-based class" = inherits(object, "Scatterer")
  )
  parallel <- .resolve_simulation_parallel(parallel, n_cores)
  normalized_model <- .validate_simulation_models(model)
  # Infer the realization count and grid axes when they are not supplied =======
  realization_setup <- .resolve_simulation_realizations(
    parameters = parameters,
    batch_by = batch_by,
    n_realizations = n_realizations
  )
  n_realizations <- realization_setup$n_realizations
  batch_by <- realization_setup$batch_by
  simulation_setup <- .prepare_simulation_grid(
    n_realizations = n_realizations,
    parameters = parameters,
    batch_by = batch_by,
    permute = permute
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
  vapply(
    .resolve_model_registry_entries(model),
    function(spec) spec$canonical,
    character(1)
  )
}

# Extract one simulation-grid row while preserving structured list-columns.
#' @noRd
.simulation_grid_row_values <- function(simulation_grid, grid_index) {
  # Resolve one value per simulation-grid column ==============================
  values <- lapply(names(simulation_grid), function(param) {
    column <- simulation_grid[[param]]
    if (is.list(column)) {
      return(column[[grid_index]])
    }

    column[[grid_index]]
  })
  names(values) <- names(simulation_grid)

  # Return the resolved parameter bundle ======================================
  values
}

# Build a one-row data frame from one simulation parameter bundle.
#' @noRd
.simulation_parameter_row_df <- function(parameter_values) {
  # Preserve structured values as list-columns in the returned results ========
  out <- data.frame(row.names = 1L)
  for (param in names(parameter_values)) {
    value <- parameter_values[[param]]
    if (is.atomic(value) &&
      length(value) == 1L &&
      (is.null(names(value)) || !any(nzchar(names(value))))) {
      out[[param]] <- value
    } else {
      out[[param]] <- I(list(value))
    }
  }

  # Return the one-row parameter data frame ===================================
  out
}

# Resolve geometry-bearing component slots available for simulation overrides.
#' @noRd
.simulation_component_slots <- function(object) {
  # Return the supported component slots present on this scatterer ============
  intersect(
    methods::slotNames(object),
    c("body", "bladder", "backbone", "shell", "fluid")
  )
}

# Resolve convenience aliases that map onto structured reforge targets.
#' @noRd
.simulation_reforge_alias_groups <- function(object) {
  # Resolve the current reforge method signature for this scatterer ===========
  valid_reforge_params <- .discover_reforge_params(class(object))
  out <- list()

  # Map body-dimension aliases onto body_target when available ================
  if ("body_target" %in% valid_reforge_params) {
    out$body_target <- list(
      aliases = c(
        length_body = "length",
        width_body = "width",
        height_body = "height",
        radius_body = "radius"
      ),
      isometric = if ("isometric_body" %in% valid_reforge_params) {
        "isometric_body"
      } else {
        NULL
      },
      legacy_conflicts = intersect(c("length", "radius"), valid_reforge_params)
    )
  }

  # Map swimbladder aliases onto swimbladder_target when available ============
  if ("swimbladder_target" %in% valid_reforge_params) {
    out$swimbladder_target <- list(
      aliases = c(
        length_swimbladder = "length",
        width_swimbladder = "width",
        height_swimbladder = "height",
        length_bladder = "length",
        width_bladder = "width",
        height_bladder = "height"
      ),
      isometric = if ("isometric_swimbladder" %in% valid_reforge_params) {
        "isometric_swimbladder"
      } else {
        NULL
      },
      legacy_conflicts = character(0)
    )
  }

  # Map backbone aliases onto backbone_target when available ==================
  if ("backbone_target" %in% valid_reforge_params) {
    out$backbone_target <- list(
      aliases = c(
        length_backbone = "length",
        width_backbone = "width",
        height_backbone = "height",
        radius_backbone = "radius"
      ),
      isometric = if ("isometric_backbone" %in% valid_reforge_params) {
        "isometric_backbone"
      } else {
        NULL
      },
      legacy_conflicts = character(0)
    )
  }

  # Return the supported alias groups for this scatterer ======================
  out
}

# Resolve one scalar numeric value used by a simulation convenience alias.
#' @noRd
.simulation_alias_scalar <- function(value, alias_name) {
  # Validate that the alias resolved to one numeric draw ======================
  if (!is.numeric(value) || length(value) != 1L || is.na(value)) {
    stop(
      "Simulation alias '", alias_name,
      "' must resolve to one non-missing numeric value per realization.",
      call. = FALSE
    )
  }

  # Return the scalar alias value =============================================
  as.numeric(value)
}

# Normalize convenience aliases onto the active reforge() argument set.
#' @noRd
.normalize_simulation_reforge_parameters <- function(object, parameter_values) {
  # Start from the explicitly supported reforge parameters ====================
  valid_reforge_params <- .discover_reforge_params(class(object))
  reforge_parameters <- parameter_values[
    names(parameter_values) %in% valid_reforge_params
  ]
  alias_groups <- .simulation_reforge_alias_groups(object)

  # Merge component-dimension aliases into structured target vectors ==========
  for (target_name in names(alias_groups)) {
    spec <- alias_groups[[target_name]]
    alias_names <- intersect(names(parameter_values), names(spec$aliases))
    legacy_names <- intersect(names(reforge_parameters), spec$legacy_conflicts)
    if (length(alias_names) == 0 && length(legacy_names) == 0) {
      next
    }

    # Reject ambiguous mixtures of explicit targets and convenience aliases ===
    if (target_name %in% names(reforge_parameters)) {
      stop(
        "Specify either '", target_name, "' or its convenience aliases (",
        paste(sprintf("'%s'", alias_names), collapse = ", "),
        "), not both.",
        call. = FALSE
      )
    }

    # Reject duplicate dimension definitions across alias and legacy inputs ===
    target_sources <- c(
      stats::setNames(alias_names, unname(spec$aliases[alias_names])),
      stats::setNames(legacy_names, legacy_names)
    )
    target_dims <- names(target_sources)
    if (anyDuplicated(target_dims)) {
      duplicated_dims <- unique(target_dims[duplicated(target_dims)])
      duplicated_sources <- unique(
        unname(target_sources[target_dims %in% duplicated_dims])
      )
      stop(
        "Simulation inputs for '", target_name,
        "' duplicate one or more dimensions: ",
        paste(sprintf("'%s'", duplicated_dims), collapse = ", "),
        ". Conflicting parameters: ",
        paste(sprintf("'%s'", duplicated_sources), collapse = ", "),
        ".",
        call. = FALSE
      )
    }

    # Resolve the structured target vector from alias and legacy values =======
    alias_value <- if (length(alias_names) > 0) {
      stats::setNames(
        vapply(
          alias_names,
          function(alias_name) {
            .simulation_alias_scalar(parameter_values[[alias_name]], alias_name)
          },
          numeric(1)
        ),
        unname(spec$aliases[alias_names])
      )
    } else {
      numeric(0)
    }
    legacy_value <- if (length(legacy_names) > 0) {
      stats::setNames(
        vapply(
          legacy_names,
          function(legacy_name) {
            .simulation_alias_scalar(reforge_parameters[[legacy_name]], legacy_name)
          },
          numeric(1)
        ),
        legacy_names
      )
    } else {
      numeric(0)
    }
    target_value <- c(alias_value, legacy_value)
    reforge_parameters[[target_name]] <- target_value

    # Drop merged legacy arguments after promoting them to the target vector ==
    if (length(legacy_names) > 0) {
      reforge_parameters <- reforge_parameters[
        !names(reforge_parameters) %in% legacy_names
      ]
    }

    # Default multi-axis convenience aliases to anisotropic scaling ===========
    if (length(target_value) > 1 &&
      !is.null(spec$isometric) &&
      !spec$isometric %in% names(reforge_parameters)) {
      reforge_parameters[[spec$isometric]] <- FALSE
    }
  }

  # Return the normalized reforge argument bundle =============================
  reforge_parameters
}

# Apply direct simulation-parameter overrides to matching scatterer components.
#' @noRd
.apply_simulation_parameter_overrides <- function(object, parameter_values) {
  # Resolve the component slots that can accept direct parameter overrides =====
  component_slots <- .simulation_component_slots(object)

  # Apply matching parameters to each geometry-bearing component ==============
  for (param in names(parameter_values)) {
    for (slot_name in component_slots) {
      component_value <- methods::slot(object, slot_name)
      if (!is.list(component_value) || !param %in% names(component_value)) {
        next
      }

      component_value[[param]] <- parameter_values[[param]]
      methods::slot(object, slot_name) <- component_value

      if ("components" %in% methods::slotNames(object)) {
        component_registry <- methods::slot(object, "components")
        if (slot_name %in% names(component_registry)) {
          component_registry[[slot_name]] <- component_value
          methods::slot(object, "components") <- component_registry
        }
      }
    }
  }

  # Return the updated scatterer ==============================================
  object
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
  parameter_values <- .simulation_grid_row_values(simulation_grid, grid_index)
  # Create working copy of scattering object ===================================
  working_object <- object
  # Reforge object, if parameters require ======================================
  # ---- Normalize explicit and convenience reforge parameters ++++++++++++++++
  reforge_parameters <- .normalize_simulation_reforge_parameters(
    working_object,
    parameter_values
  )
  # ---- Reforge the working object when geometry parameters were supplied ++++
  if (length(reforge_parameters) > 0) {
    reforge_args <- c(list(object = working_object), reforge_parameters)
    working_object <- do.call(reforge, reforge_args)
  }
  # Override parameter definitions where appropriate ===========================
  working_object <- .apply_simulation_parameter_overrides(
    working_object,
    parameter_values
  )
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
      parameter_df <- .simulation_parameter_row_df(parameter_values)
      parameter_df <- parameter_df[rep(1L, nrow(df)), , drop = FALSE]
      data.frame(
        model = rep(mod_name, nrow(df)),
        parameter_df,
        df,
        row.names = NULL,
        check.names = FALSE
      )
    }
  )
  names(model_dfs) <- names(model_results)
  # Return output ==============================================================
  model_dfs
}
