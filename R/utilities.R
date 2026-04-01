################################################################################
################################################################################
# GENERAL UTILITY FUNCTIONS
################################################################################
################################################################################
#' Calculate maximum radius from explicit value or ratio
#'
#' Internal helper function to standardize radius calculation logic across
#' shape creation functions.
#'
#' @param radius Explicit radius value (can be NULL)
#' @param length Body length
#' @param length_radius_ratio Length-to-radius ratio (can be NULL)
#' @return Calculated radius
#' @keywords internal
#' @noRd
.calculate_max_radius <- function(radius, length, length_radius_ratio) {
  # Return the explicit radius when it is already available ====================
  if (!is.null(radius)) {
    return(radius)
  }

  # Derive the radius from the supplied length-to-radius ratio =================
  if (!is.null(length_radius_ratio)) {
    return(length / length_radius_ratio)
  }

  # Stop when neither radius pathway is available ==============================
  stop(
    "Either 'radius' or 'length_radius_ratio' must be provided.",
    call. = FALSE
  )
}

#' Null-coalescing operator
#'
#' Returns first argument if not NULL and not NA, otherwise returns second.
#'
#' @param a First value
#' @param b Fallback value
#' @return a if not NULL/NA, otherwise b
#' @keywords internal
#' @noRd
`%||%` <- function(a, b) {
  # Prefer the left-hand value unless it is NULL or a scalar NA ===============
  if (is.null(a)) {
    return(b)
  }

  if (is.atomic(a) && length(a) == 1 && is.na(a)) {
    return(b)
  }

  a
}

#' Return real part if imaginary part is (numerically) zero, else return complex
#' @param x Complex value, vector, or matrix
#' @return Numeric if imaginary part is zero, else complex
#' @keywords internal
#' @noRd
`%R%` <- function(x, tol = .Machine$double.eps) {
  # Collapse purely real complex inputs back onto the numeric scale ============
  if (is.complex(x)) {
    if (all(abs(Im(x)) < tol)) {
      return(Re(x))
    } else {
      return(x)
    }
  } else {
    return(x)
  }
}

#' Resolve parameter value for simulation grid
#'
#' Internal helper to resolve parameter values in simulate_ts based on type
#' (batch, function, scalar, vector).
#'
#' @param param_name Parameter name
#' @param param_value Parameter value (scalar, vector, or function)
#' @param batch_by Batch parameter names
#' @param batch_values Batch parameter values
#' @param grid_size Simulation grid size
#' @param simulation_grid Simulation grid data frame
#' @return Resolved parameter vector
#' @keywords internal
#' @noRd
.simplify_simulation_draws <- function(values) {
  # Return atomic scalars as a simple vector when possible ====================
  scalar_atomic <- vapply(
    values,
    function(x) {
      is.atomic(x) &&
        length(x) == 1L &&
        (is.null(names(x)) || !any(nzchar(names(x))))
    },
    logical(1)
  )
  if (all(scalar_atomic)) {
    return(unlist(values, use.names = FALSE))
  }

  # Preserve structured draws as a list =======================================
  values
}

#' Resolve whether one simulation parameter should be preserved as structured
#' @param x One simulation parameter value.
#' @return Logical.
#' @keywords internal
#' @noRd
.is_structured_simulation_value <- function(x) {
  # Preserve named atomic values so reforge targets retain their dimensions ====
  if (is.atomic(x) && !is.null(names(x)) && any(nzchar(names(x)))) {
    return(TRUE)
  }

  # Preserve list-like and tabular objects as structured values ===============
  is.list(x) || is.matrix(x) || is.data.frame(x)
}

#' Normalize one batched simulation parameter into candidate values
#' @param param_name Parameter name.
#' @param param_value Parameter definition supplied to `simulate_ts()`.
#' @return Candidate values used to expand the batch grid.
#' @keywords internal
#' @noRd
.normalize_simulation_batch_values <- function(param_name, param_value) {
  # Evaluate generating functions once per batched parameter ===================
  value <- if (is.function(param_value)) {
    result <- param_value()
    if (length(result) == 0) {
      stop(
        paste0(
          "Batch parameter '", param_name, "' function must return at least 1 ",
          "valid value."
        )
      )
    }
    result
  } else {
    param_value
  }

  # Reject empty candidate sets before grid construction ======================
  if (length(value) == 0) {
    stop(
      "Batch parameter '", param_name, "' must supply at least 1 valid value.",
      call. = FALSE
    )
  }

  # Wrap structured single values so they survive batching intact =============
  if (.is_structured_simulation_value(value) && !is.list(value)) {
    return(list(value))
  }

  # Return the normalized batch candidates ====================================
  value
}

#' Resolve one non-batched simulation parameter across the full grid
#' @param param_name Parameter name.
#' @param param_value Parameter definition supplied to `simulate_ts()`.
#' @param grid_size Simulation grid size.
#' @return Resolved parameter values.
#' @keywords internal
#' @noRd
.resolve_simulation_parameter_values <- function(param_name,
                                                 param_value,
                                                 grid_size) {
  # Draw one fresh value per realization for stochastic generators =============
  if (is.function(param_value)) {
    draws <- lapply(seq_len(grid_size), function(i) param_value())
    return(.simplify_simulation_draws(draws))
  }

  # Preserve structured scalar values across the full simulation grid ==========
  if (.is_structured_simulation_value(param_value)) {
    if (!is.list(param_value) && !is.data.frame(param_value) && !is.matrix(param_value)) {
      return(rep(list(param_value), grid_size))
    }
    if (length(param_value) == 1L) {
      return(rep(param_value, grid_size))
    }
    if (length(param_value) == grid_size) {
      return(param_value)
    }
  }

  # Recycle scalar atomic inputs across the whole simulation grid =============
  if (length(param_value) == 1) {
    return(rep(param_value, grid_size))
  }

  # Accept vectors that already match the full simulation size ================
  if (length(param_value) == grid_size) {
    return(param_value)
  }

  # Reject ambiguous parameter lengths before model construction ==============
  stop(
    sprintf(
      "Length of parameter '%s' [%d] does not match number of realizations [%d].",
      param_name, length(param_value), grid_size
    ),
    call. = FALSE
  )
}

#' Resolve parameter value for simulation grid
#'
#' Internal helper to resolve parameter values in simulate_ts based on type
#' (batch, function, scalar, vector).
#'
#' @param param_name Parameter name
#' @param param_value Parameter value (scalar, vector, or function)
#' @param batch_by Batch parameter names
#' @param batch_values Batch parameter values
#' @param grid_size Simulation grid size
#' @param simulation_grid Simulation grid data frame
#' @return Resolved parameter vector
#' @keywords internal
#' @noRd
.resolve_param_value <- function(param_name, param_value, batch_by,
                                 batch_values, grid_size, simulation_grid) {
  # Resolve values from the active batch grid when batching is enabled =========
  if (!is.null(batch_by) && param_name %in% batch_by) {
    idx <- simulation_grid[[paste0(param_name, "_idx")]]
    return(batch_values[[param_name]][idx])
  }

  # Resolve non-batched values across the full simulation grid ================
  .resolve_simulation_parameter_values(
    param_name = param_name,
    param_value = param_value,
    grid_size = grid_size
  )
}

#' Resolve one slot access while traversing nested scatterer structures
#' @param layer Current traversal layer.
#' @param sub_layer Requested child element.
#' @param class_name Display name used in error messages.
#' @keywords internal
#' @noRd
.extract_slot_layer <- function(layer, sub_layer, class_name) {
  # Guard against invalid slot requests before attempting slot access ==========
  if (!methods::.hasSlot(layer, sub_layer)) {
    stop(sprintf("%s does not have slot '%s'.", class_name, sub_layer))
  }

  # Return the requested slot ==================================================
  methods::slot(layer, sub_layer)
}

#' Resolve one list element while traversing nested scatterer structures
#' @param layer Current traversal layer.
#' @param sub_layer Requested child element.
#' @keywords internal
#' @noRd
.extract_list_layer <- function(layer, sub_layer) {
  # Flag missing list elements without raising immediately =====================
  if (!sub_layer %in% names(layer)) {
    return(list(layer = layer, fail_state = TRUE))
  }

  # Return the requested list element =========================================
  list(layer = layer[[sub_layer]], fail_state = FALSE)
}

#' Resolve one matrix row or column while traversing nested structures
#' @param layer Current traversal layer.
#' @param sub_layer Requested child element.
#' @keywords internal
#' @noRd
.extract_matrix_layer <- function(layer, sub_layer) {
  # Prefer row names before column names to match the legacy accessor ==========
  is_rows <- sub_layer %in% row.names(layer)
  is_cols <- sub_layer %in% colnames(layer)

  if (is_rows) {
    return(list(layer = layer[sub_layer, ], fail_state = FALSE))
  }
  if (is_cols) {
    return(list(layer = layer[, sub_layer], fail_state = FALSE))
  }

  # Flag unresolved matrix lookups for the caller ==============================
  list(layer = layer, fail_state = TRUE)
}

#' Resolve one named-vector element while traversing nested structures
#' @param layer Current traversal layer.
#' @param sub_layer Requested child element.
#' @keywords internal
#' @noRd
.extract_vector_layer <- function(layer, sub_layer) {
  # Flag unresolved vector lookups for the caller ==============================
  if (!sub_layer %in% names(layer)) {
    return(list(layer = layer, fail_state = TRUE))
  }

  # Return the requested named-vector value ===================================
  list(layer = layer[sub_layer], fail_state = FALSE)
}

#' Advance one step through a nested scatterer/shape structure
#' @param layer Current traversal layer.
#' @param sub_layer Requested child element.
#' @keywords internal
#' @noRd
.extract_next_layer <- function(layer, sub_layer) {
  # Handle S4 scatterer slots ==================================================
  if (methods::is(layer, "Scatterer")) {
    return(list(
      layer = .extract_slot_layer(layer, sub_layer, "Scattering object"),
      fail_state = FALSE
    ))
  }

  # Handle S4 shape slots ======================================================
  if (methods::is(layer, "Shape")) {
    return(list(
      layer = .extract_slot_layer(layer, sub_layer, "Shape object"),
      fail_state = FALSE
    ))
  }

  # Handle recursive containers ================================================
  if (is.list(layer)) {
    return(.extract_list_layer(layer, sub_layer))
  }
  if (is.matrix(layer)) {
    return(.extract_matrix_layer(layer, sub_layer))
  }
  if (is.vector(layer)) {
    return(.extract_vector_layer(layer, sub_layer))
  }

  # Return the current layer unchanged when no traversal rule applies ==========
  list(layer = layer, fail_state = FALSE)
}

################################################################################
################################################################################
# Accessor functions
################################################################################
################################################################################
#' Extract nested components, slots, or matrix/vector fields from objects
#'
#' @description
#' Convenience accessor for reaching into `Scatterer` and `Shape` objects
#' without spelling out direct slot access repeatedly. `extract()` can also walk
#' through nested lists, matrices, and named vectors, which makes it useful for
#' pulling model outputs, component properties, or position-matrix fields from a
#' common interface. This is especially helpful after geometry manipulations
#' such as `brake()` and `reforge()`, where inspecting the stored body profile
#' is often the fastest way to confirm what changed.
#'
#' @param object Scatterer-class object.
#' @param feature Feature(s) of interest (e.g. body). This can either be a
#' scalar string, or a vector of names. When a vector is supplied, the function
#' recursively accesses the Scatterer object using the 'feature' vector as a
#' directory. For example, \code{feature = c("body", "rpos", "x")} would
#' extract the 'x' coordinate of the position matrix ('rpos') from the 'body'
#' scattering parameters.
#'
#' @return The extracted object, slot, list element, matrix row/column, or
#'   vector element identified by `feature`.
#'
#' @examples
#' obj <- fls_generate(
#'   shape = sphere(radius_body = 0.01, n_segments = 40),
#'   density_body = 1045,
#'   sound_speed_body = 1520
#' )
#'
#' extract(obj, "body")
#' extract(obj, c("body", "density"))
#' extract(obj, c("shape_parameters", "shape"))
#'
#' bent_obj <- brake(obj, radius_curvature = 5)
#' head(extract(bent_obj, c("body", "rpos", "z")))
#' extract(bent_obj, c("shape_parameters", "radius_curvature_ratio"))
#'
#' @keywords utility
#' @seealso [brake()], [reforge()], [translate_shape()], [reanchor_shape()],
#'   [inflate_shape()], [smooth_shape()], [resample_shape()], [flip_shape()],
#'   [offset_component()]
#' @rdname extract
#'
#' @export
extract <- function(object, feature) {
  # Initialize the traversal state used to walk the nested object ==============
  layer <- object
  fail_state <- FALSE
  accum_feature <- c()

  # Traverse one feature token at a time through the supported containers ======
  for (sub_layer in feature) {
    accum_feature <- c(accum_feature, sub_layer)
    next_layer <- .extract_next_layer(layer, sub_layer)

    if (isTRUE(next_layer$fail_state)) {
      fail_state <- TRUE
      break
    }

    layer <- next_layer$layer
  }

  # Report the deepest missing path when traversal fails =======================
  if (fail_state) {
    stop(
      sprintf("No feature '%s' ", sub_layer),
      sprintf(
        "found at the path '%s' ",
        paste(accum_feature, collapse = " -> ")
      ),
      "within the supplied Scatterer object."
    )
  } else {
    # Return the resolved nested value =========================================
    layer
  }
}
