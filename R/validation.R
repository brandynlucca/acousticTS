################################################################################
################################################################################
# Validation support functions
################################################################################
################################################################################
#' Validate a named target-dimension vector used in reforge operations
#'
#' Internal helper that checks \code{target} is a named numeric vector whose
#' names are drawn from \code{valid_dims}.
#'
#' @param target Named numeric vector of target dimensions (e.g. length, width).
#' @param target_name String used in error messages to identify the argument.
#' @param valid_dims Character vector of acceptable dimension names.
#' @return \code{target} unchanged, or \code{NULL} if \code{target} is
#'   \code{NULL}.
#' @keywords internal
#' @noRd
.validate_dimensions_target <- function(target, target_name, valid_dims) {
  if (is.null(target)) return(NULL)
  if (!is.numeric(target)) {
    stop(paste0("'", target_name, "' must be numeric."), call. = FALSE)
  }
  if (is.null(names(target)) || any(names(target) == "")) {
    stop(paste0("'", target_name, "' must be a named numeric vector."),
         call. = FALSE)
  }
  invalid <- setdiff(names(target), valid_dims)
  if (length(invalid) > 0) {
    stop(
      sprintf(
        "'%s' has invalid dimension names: %s. Valid names are: %s.",
        target_name,
        paste0("'", invalid, "'", collapse = ", "),
        paste0("'", valid_dims, "'", collapse = ", ")
      ),
      call. = FALSE
    )
  }
  target
}
#' Validate dimensional parameter (target or scale)
#'
#' Internal helper to validate named numeric vectors used in reforge operations.
#'
#' @param dims The dimensions of the parameter defined for rescaling
#' @param dims_name Name of dimension
#' @param valid_dims Character vector of valid dimension names
#' @param isometry Boolean argument defining whether to isometrically rescale 
#' all dimensions in the shape
#' @param iso_name Name of dimension used for isometric rescaling when toggled
#' @return Validated parameter (invisibly)
#' @keywords internal
#' @noRd
.validate_dimension_scaling <- function(dims,
                                        dims_name,
                                        valid_dims,
                                        isometry,
                                        iso_name) {
  # No dimensions provided
  if (is.null(dims)) return(NULL)

  # Check if dimensions are all numeric
  if (!is.numeric(dims)) {
    stop("'", dims_name, "' must be numeric.", call. = FALSE)
  }

  # Handle scalar case where `isometry=TRUE`
  if (length(dims) == 1 && isometry &&
      (is.null(names(dims)) || all(names(dims) %in% valid_dims)) ) {
    return(stats::setNames(rep(dims, length(valid_dims)), valid_dims))
  }

  # Check for missing names
  if (is.null(names(dims))) {
    stop(
      sprintf(
        paste0(
          "'%s' must either be a scalar or a named vector with dimensions: %s."
        ),
        dims_name,
        paste0("'", valid_dims, "'", collapse=", ")
      ),
      call.=FALSE
    )
  }

  # Handle vector case where `isometry=TRUE`
  if (length(dims) > 1 && isometry) {
    # ---- Get the names
    input_names <- names(dims)
    stop(
      sprintf(
        paste0(
          "'%s' contains more than 1 dimension while '%s' is TRUE. When TRUE ",
          "'%s' must be a scalar value."
        ),
        dims_name,
        iso_name,
        dims_name
      ),
      call. = FALSE
    )
  }

  # Check for invalid names
  invalid_dims <- setdiff(names(dims), valid_dims)
  if (length(invalid_dims) > 0) {
    stop(
      sprintf(
        paste0(
          "'%s' has one or more invalid dimensions: %s. ",
          "Expected dimensions are: %s."
        ),
        dims_name,
        paste0("'", invalid_dims, "'", collapse=", "),
        paste0("'", valid_dims, "'", collapse=", ")
      ),
      call. = FALSE
    )
  }

  # Handle scalar case where `isometry=FALSE`
  if (!isometry) {
    if (is.null(names(dims))) {
      stop(
        sprintf(
          paste0(
            "When '%s' is FALSE, '%s' must be a named vector with at least ",
            "one of the following dimensions: %s."
          ),
          iso_name,
          dims_name,
          paste0("'", valid_dims, "'", collapse=", ")
        ),
        call. = FALSE
      )
    }
    # ---- Initialize vector
    output_dims <- stats::setNames(rep(1, length(valid_dims)), valid_dims)
    # ---- Override and return
    output_dims[names(dims)] <- dims
    return(output_dims)
  }

  dims
}
#' Validate elastic moduli inputs
#'
#' Internal helper to validate that sufficient elastic moduli are provided
#' for calculations.
#'
#' @param ... Named elastic parameters (K, E, G, nu)
#' @param min_required Minimum number of non-NULL parameters
#' @param param_name Name of parameter being calculated (for error message)
#' @return Character vector of provided parameter names
#' @keywords internal
#' @noRd
.validate_elastic_inputs <- function(..., min_required = 2, param_name = NULL) {
  inputs <- list(...)
  provided <- !vapply(inputs, is.null, logical(1))

  if (sum(provided) < min_required) {
    msg <- sprintf(
      "At least %d elasticity moduli values are required",
      min_required
    )
    if (!is.null(param_name)) {
      msg <- paste(msg, "to calculate", param_name)
    }
    stop(msg, ".", call. = FALSE)
  }

  names(provided)[provided]
}

#' Validate brake/bending parameters
#' @param body_df Dataframe containing shape information
#' @param radius_curvature Curvature defined by the radius of an overlapping 
#' osculating circle. This can either be an absolute or ratio quantity.
#' @param mode Domain of radius of curvature quantity that can either be a 
#' ratio or an absolute quantity
#' @keywords internal
#' @noRd
.validate_brake_params <- function(body_df, radius_curvature, mode) {
  if (!is.list(body_df) || is.null(body_df$rpos) || !is.matrix(body_df$rpos)) {
    stop("Body shape information must be a list with a matrix element 'rpos'.",
         call. = FALSE)
  }
  if (!is.numeric(radius_curvature) || radius_curvature <= 0) {
    stop("Radius of curvature must be a positive-only, real number.",
         call. = FALSE)
  }
  if (!mode %in% c("ratio", "measurement")) {
    stop("Radius-of-curvature 'mode' must be either 'ratio' or 'measurement'.",
         call. = FALSE)
  }
  invisible(TRUE)
}

#' Validate and extract model data for plotting
#' @param object Scatterer object
#' @param model_name Model name (e.g., "calibration", "dwba")
#' @return List with model, parameters, and shape data
#' @keywords internal
#' @noRd
.validate_and_extract_model <- function(object, model_name) {
  models <- extract(object, "model")
  if (length(models) == 0) {
    stop("ERROR: no model results detected in object.", call. = FALSE)
  }

  if (!model_name %in% names(models)) {
    stop(
      sprintf("Model '%s' not found. Available: %s",
              model_name, paste(names(models), collapse = ", ")),
      call. = FALSE
    )
  }

  list(
    model = models[[model_name]],
    parameters = extract(object, "model_parameters")[[model_name]],
    shape = extract(object, "shape_parameters") %||% extract(object, "body")
  )
}
