################################################################################
################################################################################
# SCATTERER CONSTRUCTION UTILITIES
################################################################################
################################################################################
#' Build a standard metadata list for scatterer constructors
#' @param ID Optional identifier supplied by the user.
#' @param default_id Default identifier used when `ID` is not supplied.
#' @param extra Optional named list of additional metadata fields.
#' @keywords internal
#' @noRd
.scatterer_metadata <- function(ID = NULL,
                                default_id = "UID",
                                extra = list()) {
  c(
    list(ID = ifelse(!is.null(ID), ID, default_id)),
    extra
  )
}

#' Normalize scatterer-constructor units to SI defaults
#' @param theta_units Angular units supplied to the constructor.
#' @param length_units Length units supplied to the constructor.
#' @param radius_units Radius units supplied to the constructor.
#' @param diameter_units Diameter units supplied to the constructor.
#' @param context Constructor label for warnings.
#' @keywords internal
#' @noRd
.normalize_scatterer_units <- function(theta_units = "radians",
                                       length_units = "m",
                                       radius_units = "m",
                                       diameter_units = "m",
                                       context = "scatterer") {
  list(
    theta_units = "radians",
    length_units = "m",
    radius_units = "m",
    diameter_units = "m"
  )
}

#' Validate contrast-versus-absolute material inputs for one component
#' @param contrast Density or sound-speed contrast input.
#' @param absolute Absolute density or sound-speed input.
#' @param contrast_name Contrast argument name.
#' @param absolute_name Absolute-property argument name.
#' @param component_label Human-readable component label.
#' @param require_one Logical; whether at least one of the pair must be
#'   supplied.
#' @keywords internal
#' @noRd
.validate_material_pair <- function(contrast,
                                    absolute,
                                    contrast_name,
                                    absolute_name,
                                    component_label = "component",
                                    require_one = TRUE) {
  if (!is.null(contrast) && !is.null(absolute)) {
    stop(
      "Cannot specify both ", contrast_name, " and ", absolute_name,
      " for the ", component_label, ".",
      call. = FALSE
    )
  }

  if (require_one && is.null(contrast) && is.null(absolute)) {
    stop(
      "Supply either ", contrast_name, " or ", absolute_name,
      " for the ", component_label, ".",
      call. = FALSE
    )
  }

  invisible(TRUE)
}

#' Validate density and sound-speed inputs for one scatterer component
#' @param g Density contrast.
#' @param density Absolute density.
#' @param h Sound-speed contrast.
#' @param sound_speed Absolute sound speed.
#' @param component_label Human-readable component label.
#' @param require_one Logical; whether both property pairs must be supplied.
#' @keywords internal
#' @noRd
.validate_component_material_inputs <- function(g = NULL,
                                                density = NULL,
                                                h = NULL,
                                                sound_speed = NULL,
                                                component_label = "component",
                                                require_one = TRUE) {
  .validate_material_pair(
    contrast = g,
    absolute = density,
    contrast_name = paste0("g_", component_label),
    absolute_name = paste0("density_", component_label),
    component_label = component_label,
    require_one = require_one
  )
  .validate_material_pair(
    contrast = h,
    absolute = sound_speed,
    contrast_name = paste0("h_", component_label),
    absolute_name = paste0("sound_speed_", component_label),
    component_label = component_label,
    require_one = require_one
  )

  invisible(TRUE)
}

#' Resolve a canonical single-component shape input for scatterer constructors
#' @param shape Shape object or legacy shape specifier.
#' @param arguments Evaluated constructor arguments.
#' @keywords internal
#' @noRd
.resolve_scatterer_shape_input <- function(shape, arguments) {
  if (methods::is(shape, "Shape")) {
    return(shape)
  }

  .resolve_shape(shape, arguments)
}

#' Resolve an explicit component shape or fall back to a legacy arbitrary shape
#' @param shape_obj Explicit `Shape` object, if supplied.
#' @param suffix Component suffix such as `"body"` or `"bladder"`.
#' @param x Along-axis coordinates.
#' @param w Width coordinates.
#' @param zU Upper profile coordinates.
#' @param zL Lower profile coordinates.
#' @keywords internal
#' @noRd
.resolve_profile_component_shape <- function(shape_obj = NULL,
                                             suffix = "body",
                                             x = NULL,
                                             w = NULL,
                                             zU = NULL,
                                             zL = NULL) {
  if (methods::is(shape_obj, "Shape")) {
    return(shape_obj)
  }

  args <- list(
    x = x,
    w = w,
    zU = zU,
    zL = zL
  )
  names(args) <- paste0(names(args), "_", suffix)

  do.call(arbitrary, args)
}

#' Internal geometry-schema registry
#' @keywords internal
#' @noRd
.geometry_contract_schema <- function() {
  list(
    shape_column_major = list(
      storage = "n_points x n_fields",
      x = c("x", "x_body", "x_bladder"),
      w = c("w", "w_body", "w_bladder", "y", "y_body", "y_bladder"),
      zU = c("zU", "zU_body", "zU_bladder"),
      zL = c("zL", "zL_body", "zL_bladder"),
      radius = c("a", "radius", "radius_body", "radius_bladder")
    ),
    profile_row_major = list(
      storage = "n_fields x n_points",
      x = c("x", "x_body", "x_bladder"),
      w = c("w", "w_body", "w_bladder", "y", "y_body", "y_bladder"),
      zU = c("zU", "zU_body", "zU_bladder"),
      zL = c("zL", "zL_body", "zL_bladder"),
      radius = c("a", "radius", "radius_body", "radius_bladder")
    )
  )
}

#' Validate that a position matrix satisfies the internal geometry contract
#' @param position_matrix Numeric position matrix.
#' @param storage One of `"shape_column_major"` or `"profile_row_major"`.
#' @param context Character label used in error messages.
#' @keywords internal
#' @noRd
.validate_geometry_contract <- function(position_matrix,
                                        storage = c(
                                          "shape_column_major",
                                          "profile_row_major"
                                        ),
                                        context = "Geometry") {
  storage <- match.arg(storage)

  if (!is.matrix(position_matrix)) {
    stop(context, " must be stored as a matrix.", call. = FALSE)
  }

  schema <- .geometry_contract_schema()[[storage]]
  axis_names <- if (storage == "shape_column_major") {
    colnames(position_matrix)
  } else {
    rownames(position_matrix)
  }
  axis_label <- if (storage == "shape_column_major") "column" else "row"

  if (is.null(axis_names)) {
    stop(
      context, " is missing ", axis_label,
      " names required by the internal geometry contract.",
      call. = FALSE
    )
  }

  if (!any(schema$x %in% axis_names)) {
    stop(
      context, " must define an x-axis ", axis_label,
      " (one of: ", paste(schema$x, collapse = ", "), ").",
      call. = FALSE
    )
  }

  has_zU <- any(schema$zU %in% axis_names)
  has_zL <- any(schema$zL %in% axis_names)
  if (xor(has_zU, has_zL)) {
    stop(
      context, " must define both upper and lower height coordinates ",
      "(zU and zL) whenever either is present.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}

#' Resolve the number of intervals represented by a position matrix
#' @param position_matrix Numeric position matrix.
#' @param row_major Logical; whether nodes are stored in columns.
#' @keywords internal
#' @noRd
.shape_segment_count <- function(position_matrix, row_major = FALSE) {
  if (!is.matrix(position_matrix)) {
    stop("'position_matrix' must be a matrix.", call. = FALSE)
  }

  if (row_major) {
    return(ncol(position_matrix) - 1L)
  }

  nrow(position_matrix) - 1L
}

#' Resolve a user-facing shape label from a shape object or request
#' @param shape_input Shape object.
#' @param requested_shape Optional original shape request.
#' @keywords internal
#' @noRd
.scatterer_shape_name <- function(shape_input, requested_shape = NULL) {
  if (is.character(requested_shape) && requested_shape == "arbitrary") {
    return("Arbitrary")
  }

  paste0(class(shape_input))
}

#' Resolve a representative scalar radius for scatterer metadata
#' @param shape_input Optional shape object.
#' @param position_matrix Optional position matrix.
#' @param shape_parameters Optional shape-parameter list.
#' @param error_context Context label used in fallback errors.
#' @keywords internal
#' @noRd
.shape_max_radius <- function(shape_input = NULL,
                              position_matrix = NULL,
                              shape_parameters = NULL,
                              error_context = "shape") {
  if (!is.null(shape_input)) {
    position_matrix <- acousticTS::extract(shape_input, "position_matrix")
    shape_parameters <- acousticTS::extract(shape_input, "shape_parameters")
  }

  if (!is.null(shape_parameters$radius) &&
      length(shape_parameters$radius) == 1 &&
      !is.na(shape_parameters$radius)) {
    return(as.numeric(shape_parameters$radius)[1])
  }

  max(
    .shape_radius_profile(
      position_matrix = position_matrix,
      shape_parameters = shape_parameters,
      error_context = error_context
    ),
    na.rm = TRUE
  )
}

#' Append class-specific shape metadata used by scatterer constructors
#' @param shape_params Base shape-parameter list.
#' @param shape_input Shape object.
#' @param requested_shape Optional original shape request.
#' @keywords internal
#' @noRd
.append_shape_specific_parameters <- function(shape_params,
                                              shape_input,
                                              requested_shape = NULL) {
  input_params <- acousticTS::extract(shape_input, "shape_parameters")

  if (methods::is(shape_input, "ProlateSpheroid")) {
    shape_params$semimajor_length <- input_params$semimajor_length
    shape_params$semiminor_length <- input_params$semiminor_length
  }

  if (methods::is(shape_input, "Cylinder")) {
    shape_params$radius_curvature_ratio <- input_params$radius_curvature_ratio
    if ("taper_order" %in% names(input_params)) {
      shape_params$taper_order <- input_params$taper_order
    }
  }

  if (methods::is(shape_input, "Sphere")) {
    if ("radius_shape" %in% names(input_params)) {
      shape_params$radius_shape <- input_params$radius_shape
    }
  }

  shape_params$shape <- .scatterer_shape_name(
    shape_input,
    requested_shape = requested_shape
  )

  shape_params
}

#' Build the common shape-parameter block used by several scatterer constructors
#' @param shape_input Shape object.
#' @param requested_shape Optional original shape request.
#' @param error_context Context label for radius resolution.
#' @param extra_units Optional named list of unit fields.
#' @keywords internal
#' @noRd
.shape_common_parameters <- function(shape_input,
                                     requested_shape = NULL,
                                     error_context = "shape",
                                     extra_units = list()) {
  position_matrix <- acousticTS::extract(shape_input, "position_matrix")
  .validate_geometry_contract(
    position_matrix,
    storage = "shape_column_major",
    context = paste0(error_context, " position matrix")
  )

  params <- list(
    length = .shape_length(position_matrix = position_matrix),
    radius = .shape_max_radius(
      shape_input = shape_input,
      error_context = error_context
    ),
    n_segments = .shape_segment_count(position_matrix)
  )

  params <- .append_shape_specific_parameters(
    params,
    shape_input,
    requested_shape = requested_shape
  )

  c(params, extra_units)
}

#' Define the canonical x/w/zU/zL schema for profile-style components
#' @param suffix Component suffix such as `"body"` or `"bladder"`.
#' @keywords internal
#' @noRd
.shape_component_schema <- function(suffix = "body") {
  list(
    rows = c(
      paste0("x_", suffix),
      paste0("w_", suffix),
      paste0("zU_", suffix),
      paste0("zL_", suffix)
    ),
    x = c(paste0("x_", suffix), "x_body", "x_bladder", "x"),
    w = c(
      paste0("w_", suffix),
      "w_body", "w_bladder", "w",
      paste0("y_", suffix), "y_body", "y_bladder", "y"
    ),
    zU = c(paste0("zU_", suffix), "zU_body", "zU_bladder", "zU"),
    zL = c(paste0("zL_", suffix), "zL_body", "zL_bladder", "zL")
  )
}

#' Build a canonical row-major profile matrix from component axes
#' @param x Along-body coordinates.
#' @param w Width coordinates.
#' @param zU Upper profile coordinates.
#' @param zL Lower profile coordinates.
#' @param suffix Component suffix such as `"body"` or `"bladder"`.
#' @keywords internal
#' @noRd
.canonical_profile_matrix <- function(x, w, zU, zL, suffix = "body") {
  schema <- .shape_component_schema(suffix)
  rpos <- rbind(x, w, zU, zL)
  rownames(rpos) <- schema$rows
  rpos
}

#' Extract a named column from a shape position matrix
#' @param position_matrix Shape position matrix.
#' @param candidates Candidate column names in priority order.
#' @param default Optional default vector returned when no candidate exists.
#' @keywords internal
#' @noRd
.extract_shape_component_column <- function(position_matrix,
                                           candidates,
                                           default = NULL) {
  for (nm in candidates) {
    if (!is.null(colnames(position_matrix)) && nm %in% colnames(position_matrix)) {
      return(position_matrix[, nm])
    }
  }

  if (!is.null(default)) {
    return(default)
  }

  position_matrix[, 1]
}

#' Extract a named row from a row-major profile matrix
#' @param position_matrix Row-major position matrix.
#' @param candidates Candidate row names in priority order.
#' @param default Optional default vector returned when no candidate exists.
#' @keywords internal
#' @noRd
.extract_shape_component_row <- function(position_matrix,
                                         candidates,
                                         default = NULL) {
  for (nm in candidates) {
    if (!is.null(rownames(position_matrix)) && nm %in% rownames(position_matrix)) {
      return(position_matrix[nm, ])
    }
  }

  if (!is.null(default)) {
    return(default)
  }

  position_matrix[1, ]
}

#' Convert a shape object into canonical row-major profile coordinates
#' @param shape_obj Shape object.
#' @param suffix Component suffix such as `"body"` or `"bladder"`.
#' @keywords internal
#' @noRd
.as_component_rpos <- function(shape_obj, suffix = "body") {
  position_matrix <- acousticTS::extract(shape_obj, "position_matrix")
  .validate_geometry_contract(
    position_matrix,
    storage = "shape_column_major",
    context = paste0(toupper(substring(suffix, 1, 1)), substring(suffix, 2), " shape")
  )
  schema <- .shape_component_schema(suffix)

  x <- .extract_shape_component_column(position_matrix, schema$x)
  w <- .extract_shape_component_column(
    position_matrix,
    schema$w,
    default = rep(0, length(x))
  )
  zU <- .extract_shape_component_column(
    position_matrix,
    schema$zU,
    default = rep(NA_real_, length(x))
  )
  zL <- .extract_shape_component_column(
    position_matrix,
    schema$zL,
    default = rep(NA_real_, length(x))
  )

  .canonical_profile_matrix(x, w, zU, zL, suffix = suffix)
}

#' Extract canonical x/w/zU/zL fields from a row-major profile matrix
#' @param position_matrix Row-major profile matrix.
#' @keywords internal
#' @noRd
.profile_fields <- function(position_matrix) {
  .validate_geometry_contract(
    position_matrix,
    storage = "profile_row_major",
    context = "Profile matrix"
  )

  list(
    x = .extract_shape_component_row(
      position_matrix,
      c("x", "x_body", "x_bladder")
    ),
    w = .extract_shape_component_row(
      position_matrix,
      c("w", "w_body", "w_bladder", "y", "y_body", "y_bladder"),
      default = rep(0, ncol(position_matrix))
    ),
    zU = .extract_shape_component_row(
      position_matrix,
      c("zU", "zU_body", "zU_bladder"),
      default = rep(NA_real_, ncol(position_matrix))
    ),
    zL = .extract_shape_component_row(
      position_matrix,
      c("zL", "zL_body", "zL_bladder"),
      default = rep(NA_real_, ncol(position_matrix))
    )
  )
}

#' Validate that a component profile contains at least one interval
#' @param rpos Row-major position matrix.
#' @param component_name Human-readable component label.
#' @keywords internal
#' @noRd
.validate_component_rpos <- function(rpos, component_name = "Component") {
  if (ncol(rpos) < 2) {
    stop(
      component_name,
      " shape must have at least two points (>=1 segment).",
      call. = FALSE
    )
  }

  invisible(rpos)
}

#' Build a fluid-like profile component from a shape object
#' @param shape_obj Shape object.
#' @param suffix Component suffix such as `"body"` or `"bladder"`.
#' @param theta Component orientation.
#' @param density Optional absolute density.
#' @param sound_speed Optional absolute sound speed.
#' @param g Optional density contrast.
#' @param h Optional sound-speed contrast.
#' @param extra Optional named list of extra component fields.
#' @keywords internal
#' @noRd
.build_fluid_profile_component <- function(shape_obj,
                                           suffix = "body",
                                           theta = pi / 2,
                                           density = NULL,
                                           sound_speed = NULL,
                                           g = NULL,
                                           h = NULL,
                                           extra = list()) {
  rpos <- .as_component_rpos(shape_obj, suffix = suffix)
  .validate_component_rpos(
    rpos,
    component_name = paste0(
      toupper(substring(suffix, 1, 1)),
      substring(suffix, 2)
    )
  )

  c(
    list(
      rpos = rpos,
      sound_speed = sound_speed,
      density = density,
      g = g,
      h = h,
      theta = theta
    ),
    extra
  )
}

#' Build a row-major fluid component from a shape object
#' @param shape_obj Shape object.
#' @param theta Component orientation.
#' @param density Optional absolute density.
#' @param sound_speed Optional absolute sound speed.
#' @param g Optional density contrast.
#' @param h Optional sound-speed contrast.
#' @param radius Optional radius metadata override.
#' @param extra Optional named list of extra component fields.
#' @keywords internal
#' @noRd
.build_row_major_fluid_component <- function(shape_obj,
                                             theta = pi / 2,
                                             density = NULL,
                                             sound_speed = NULL,
                                             g = NULL,
                                             h = NULL,
                                             radius = NULL,
                                             extra = list()) {
  position_matrix <- acousticTS::extract(shape_obj, "position_matrix")
  shape_parameters <- acousticTS::extract(shape_obj, "shape_parameters")

  .validate_geometry_contract(
    position_matrix,
    storage = "shape_column_major",
    context = "Shape position matrix"
  )

  if (is.null(radius)) {
    radius <- shape_parameters$radius
  }

  c(
    list(
      rpos = t(position_matrix),
      radius = radius,
      theta = theta,
      g = g,
      h = h,
      density = density,
      sound_speed = sound_speed
    ),
    extra
  )
}

#' Build a column-major fluid component from a shape object
#' @inheritParams .build_row_major_fluid_component
#' @keywords internal
#' @noRd
.build_column_major_fluid_component <- function(shape_obj,
                                                theta = pi / 2,
                                                density = NULL,
                                                sound_speed = NULL,
                                                g = NULL,
                                                h = NULL,
                                                radius = NULL,
                                                extra = list()) {
  position_matrix <- acousticTS::extract(shape_obj, "position_matrix")
  shape_parameters <- acousticTS::extract(shape_obj, "shape_parameters")

  .validate_geometry_contract(
    position_matrix,
    storage = "shape_column_major",
    context = "Shape position matrix"
  )

  if (is.null(radius)) {
    radius <- shape_parameters$radius
  }

  c(
    list(
      rpos = position_matrix,
      radius = radius,
      theta = theta,
      g = g,
      h = h,
      density = density,
      sound_speed = sound_speed
    ),
    extra
  )
}

#' Build a row-major elastic component from a shape object
#' @param shape_obj Shape object.
#' @param theta Component orientation.
#' @param density Optional absolute density.
#' @param sound_speed_longitudinal Optional longitudinal wave speed.
#' @param sound_speed_transversal Optional transversal wave speed.
#' @param radius Optional radius metadata override.
#' @param position_matrix_override Optional position matrix override.
#' @param extra Optional named list of extra component fields.
#' @keywords internal
#' @noRd
.build_row_major_elastic_component <- function(shape_obj,
                                               theta = pi / 2,
                                               density = NULL,
                                               sound_speed_longitudinal = NULL,
                                               sound_speed_transversal = NULL,
                                               radius = NULL,
                                               position_matrix_override = NULL,
                                               extra = list()) {
  position_matrix <- if (!is.null(position_matrix_override)) {
    position_matrix_override
  } else {
    acousticTS::extract(shape_obj, "position_matrix")
  }
  shape_parameters <- acousticTS::extract(shape_obj, "shape_parameters")

  .validate_geometry_contract(
    position_matrix,
    storage = "shape_column_major",
    context = "Elastic shape position matrix"
  )

  if (is.null(radius)) {
    radius <- shape_parameters$radius
  }

  c(
    list(
      rpos = t(position_matrix),
      radius = radius,
      theta = theta,
      density = density,
      sound_speed_longitudinal = sound_speed_longitudinal,
      sound_speed_transversal = sound_speed_transversal
    ),
    extra
  )
}
