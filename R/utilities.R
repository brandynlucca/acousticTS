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
  # Return the left-hand value unless it is NULL or NA =========================
  if (!is.null(a) && !is.na(a)) a else b
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
.resolve_param_value <- function(param_name, param_value, batch_by,
                                 batch_values, grid_size, simulation_grid) {
  # Resolve values from the active batch grid when batching is enabled =========
  if (!is.null(batch_by) && param_name %in% batch_by) {
    idx <- simulation_grid[[paste0(param_name, "_idx")]]
    return(batch_values[[param_name]][idx])
  }

  # Draw one fresh value per realization for stochastic generators =============
  if (is.function(param_value)) {
    return(replicate(grid_size, param_value()))
  }

  # Recycle scalar inputs across the whole simulation grid =====================
  if (length(param_value) == 1) {
    return(rep(param_value, grid_size))
  }

  # Accept vectors that already match the full simulation size =================
  if (length(param_value) == grid_size) {
    return(param_value)
  }

  # Reject ambiguous vector lengths before model construction ==================
  sim_type <- if (is.null(batch_by)) "realizations" else "batched realizations"
  stop(
    sprintf(
      "Length of parameter '%s' [%d] does not match number of %s [%d].",
      param_name, length(param_value), sim_type, grid_size
    ),
    call. = FALSE
  )
}

################################################################################
################################################################################
# Accessor functions
################################################################################
################################################################################
#' Extract nested components, slots, or matrix/vector fields from package objects
#'
#' @description
#' Convenience accessor for reaching into `Scatterer` and `Shape` objects
#' without spelling out direct slot access repeatedly. `extract()` can also walk
#' through nested lists, matrices, and named vectors, which makes it useful for
#' pulling model outputs, component properties, or position-matrix fields from a
#' common interface.
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
#' @keywords utility
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
    if (methods::is(layer, "Scatterer")) {
      if (!methods::.hasSlot(layer, sub_layer)) {
        stop(sprintf("Scattering object does not have slot '%s'.", sub_layer))
      }

      layer <- methods::slot(layer, sub_layer)
    } else if (methods::is(layer, "Shape")) {
      if (!methods::.hasSlot(layer, sub_layer)) {
        stop(sprintf("Shape object does not have slot '%s'.", sub_layer))
      }

      layer <- methods::slot(layer, sub_layer)
    } else if (is.list(layer)) {
      if (! sub_layer %in% names(layer)) {
        fail_state <- TRUE
        break
      }
      layer <- layer[[sub_layer]]
    } else if (is.matrix(layer)) {
      is_rows <- sub_layer %in% row.names(layer)
      is_cols <- sub_layer %in% colnames(layer)
      if (is_rows) {
        layer <- layer[sub_layer, ]
      } else if (is_cols) {
        layer <- layer[, sub_layer]
      } else {
        fail_state <- TRUE
        break
      }
    } else if (is.vector(layer)) {
      if (!sub_layer %in% names(layer)) {
        fail_state <- TRUE
        break
      }
      layer <- layer[sub_layer]
    }
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
