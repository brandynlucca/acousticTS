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
  if (!is.null(radius)) {
    return(radius)
  }

  if (!is.null(length_radius_ratio)) {
    return(length / length_radius_ratio)
  }

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
  if (!is.null(a) && !is.na(a)) a else b
}

#' Return real part if imaginary part is (numerically) zero, else return complex
#' @param x Complex value, vector, or matrix
#' @return Numeric if imaginary part is zero, else complex
#' @keywords internal
#' @noRd
`%R%` <- function(x, tol = .Machine$double.eps) {
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
  if (!is.null(batch_by) && param_name %in% batch_by) {
    idx <- simulation_grid[[paste0(param_name, "_idx")]]
    return(batch_values[[param_name]][idx])
  }

  if (is.function(param_value)) {
    return(replicate(grid_size, param_value()))
  }

  if (length(param_value) == 1) {
    return(rep(param_value, grid_size))
  }

  if (length(param_value) == grid_size) {
    return(param_value)
  }

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
#' Primary accessor function for extracting specific features and layers from
#' Scatterer objects
#'
#' @param object Scatterer-class object.
#' @param feature Feature(s) of interest (e.g. body). This can either be a
#' scalar string, or a vector of names. When a vector is supplied, the function
#' recursively accesses the Scatterer object using the 'feature' vector as a
#' directory. For example, \code{feature = c("body", "rpos", "x")} would
#' extract the 'x' coordinate of the position matrix ('rpos') from the 'body'
#' scattering parameters.
#'
#' @keywords utility
#' @rdname extract
#'
#' @export
extract <- function(object, feature) {
  layer <- object
  fail_state <- FALSE
  accum_feature <- c()

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
    layer
  }
}
