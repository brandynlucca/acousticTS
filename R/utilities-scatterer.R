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
  # Build the standard metadata list with a resolved identifier ================
  c(
    list(ID = ifelse(!is.null(ID), ID, default_id)),
    extra
  )
}

#' Warn that one scatterer-constructor pathway is compatibility-only
#' @param message Warning text.
#' @keywords internal
#' @noRd
.warn_scatterer_constructor_compatibility <- function(message) {
  # Emit the compatibility warning immediately and return invisibly ============
  warning(message, call. = FALSE, immediate. = TRUE)
  invisible(NULL)
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
  # Collect the supplied unit arguments and their SI defaults ==================
  unit_args <- list(
    theta_units = theta_units %||% "radians",
    length_units = length_units %||% "m",
    radius_units = radius_units %||% "m",
    diameter_units = diameter_units %||% "m"
  )
  si_defaults <- list(
    theta_units = "radians",
    length_units = "m",
    radius_units = "m",
    diameter_units = "m"
  )

  # Warn on legacy non-SI inputs while forcing the SI contract =================
  for (nm in names(unit_args)) {
    if (!identical(unit_args[[nm]], si_defaults[[nm]])) {
      .warn_scatterer_constructor_compatibility(
        paste0(
          "'", nm, "' in ", context, " constructors is deprecated and ignored.",
          " Supply scatterer geometry in meters and orientations in radians."
        )
      )
    }
  }

  # Return the normalized constructor-unit bundle ==============================
  list(
    theta_units = "radians",
    length_units = "m",
    radius_units = "m",
    diameter_units = "m"
  )
}

#' Detect whether constructor arguments contain explicit profile coordinates
#' @param arguments Evaluated constructor arguments.
#' @keywords internal
#' @noRd
.has_explicit_profile_coordinates <- function(arguments) {
  # Exit early when the constructor call has no named coordinate inputs ========
  coordinate_names <- names(arguments)
  if (is.null(coordinate_names)) {
    return(FALSE)
  }

  has_x <- grepl("^x_", coordinate_names)
  if (!any(has_x)) {
    return(FALSE)
  }

  # Detect whether any x-style coordinate was actually supplied ================
  any(vapply(arguments[has_x], function(x) !is.null(x), logical(1)))
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
  # Reject simultaneous contrast and absolute-property inputs ==================
  if (!is.null(contrast) && !is.null(absolute)) {
    stop(
      "Cannot specify both ", contrast_name, " and ", absolute_name,
      " for the ", component_label, ".",
      call. = FALSE
    )
  }

  # Enforce that one input pathway is supplied when required ===================
  if (require_one && is.null(contrast) && is.null(absolute)) {
    stop(
      "Supply either ", contrast_name, " or ", absolute_name,
      " for the ", component_label, ".",
      call. = FALSE
    )
  }

  # Return invisibly after successful validation ===============================
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
  # Validate the density-style inputs for the requested component ==============
  .validate_material_pair(
    contrast = g,
    absolute = density,
    contrast_name = paste0("g_", component_label),
    absolute_name = paste0("density_", component_label),
    component_label = component_label,
    require_one = require_one
  )
  # Validate the sound-speed-style inputs for the requested component ==========
  .validate_material_pair(
    contrast = h,
    absolute = sound_speed,
    contrast_name = paste0("h_", component_label),
    absolute_name = paste0("sound_speed_", component_label),
    component_label = component_label,
    require_one = require_one
  )

  # Return invisibly after successful validation ===============================
  invisible(TRUE)
}

#' Resolve a canonical single-component shape input for scatterer constructors
#' @param shape Shape object or legacy shape specifier.
#' @param arguments Evaluated constructor arguments.
#' @keywords internal
#' @noRd
.resolve_scatterer_shape_input <- function(shape, arguments) {
  # Accept an explicit Shape object without further coercion ===================
  if (methods::is(shape, "Shape")) {
    return(shape)
  }

  # Detect whether the constructor supplied explicit profile coordinates =======
  has_profile_coordinates <- .has_explicit_profile_coordinates(arguments)
  if (is.null(shape)) {
    if (has_profile_coordinates) {
      return(.resolve_shape("arbitrary", arguments))
    }
    stop(
      "Supply 'shape' as a pre-built Shape object, or provide explicit ",
      "profile coordinates such as x_body/zU_body/zL_body.",
      call. = FALSE
    )
  }

  # Warn when the legacy character-dispatch pathway is used ====================
  if (is.character(shape) &&
      !(identical(shape, "arbitrary") && has_profile_coordinates)) {
    .warn_scatterer_constructor_compatibility(
      paste0(
        "Character-based 'shape' dispatch in scatterer constructors is ",
        "deprecated. Build the geometry first with `sphere()`, `cylinder()`, ",
        "`prolate_spheroid()`, `oblate_spheroid()`, or `arbitrary()`, then ",
        "pass the resulting Shape object."
      )
    )
  }

  # Resolve the requested shape through the shared shape factory ===============
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
  # Accept an explicit component shape when one is already available ===========
  if (methods::is(shape_obj, "Shape")) {
    return(shape_obj)
  }

  # Rebuild the component shape from explicit profile coordinates ==============
  args <- list(
    x = x,
    w = w,
    zU = zU,
    zL = zL
  )
  names(args) <- paste0(names(args), "_", suffix)

  # Return the reconstructed arbitrary component shape =========================
  do.call(arbitrary, args)
}

#' Internal geometry-schema registry
#' @keywords internal
#' @noRd
.geometry_contract_schema <- function() {
  # Return the shared geometry schema used across constructor workflows ========
  list(
    shape_column_major = list(
      storage = "n_points x n_fields",
      x = c(
        "x", "x_body", "x_bladder", "x_shell", "x_fluid", "x_backbone"
      ),
      w = c(
        "w", "w_body", "w_bladder", "w_shell", "w_fluid", "w_backbone",
        "y", "y_body", "y_bladder", "y_shell", "y_fluid", "y_backbone"
      ),
      zU = c(
        "zU", "zU_body", "zU_bladder", "zU_shell", "zU_fluid",
        "zU_backbone"
      ),
      zL = c(
        "zL", "zL_body", "zL_bladder", "zL_shell", "zL_fluid",
        "zL_backbone"
      ),
      radius = c(
        "a", "radius", "radius_body", "radius_bladder", "radius_shell",
        "radius_fluid", "radius_backbone"
      )
    ),
    profile_row_major = list(
      storage = "n_fields x n_points",
      x = c(
        "x", "x_body", "x_bladder", "x_shell", "x_fluid", "x_backbone"
      ),
      w = c(
        "w", "w_body", "w_bladder", "w_shell", "w_fluid", "w_backbone",
        "y", "y_body", "y_bladder", "y_shell", "y_fluid", "y_backbone"
      ),
      zU = c(
        "zU", "zU_body", "zU_bladder", "zU_shell", "zU_fluid",
        "zU_backbone"
      ),
      zL = c(
        "zL", "zL_body", "zL_bladder", "zL_shell", "zL_fluid",
        "zL_backbone"
      ),
      radius = c(
        "a", "radius", "radius_body", "radius_bladder", "radius_shell",
        "radius_fluid", "radius_backbone"
      )
    )
  )
}

#' Resolve a canonical geometry axis from a named position matrix
#' @param position_matrix Numeric position matrix.
#' @param axis One of `"x"`, `"w"`, `"zU"`, `"zL"`, or `"radius"`.
#' @param row_major Logical; whether axis names are stored on rows.
#' @param default Optional fallback vector returned when no named axis exists.
#' @param context Context label used in validation messages.
#' @keywords internal
#' @noRd
.geometry_axis_values <- function(position_matrix,
                                  axis = c("x", "w", "zU", "zL", "radius"),
                                  row_major = FALSE,
                                  default = NULL,
                                  context = "Geometry") {
  # Resolve the requested axis name and validate the matrix input ==============
  axis <- match.arg(axis)

  if (!is.matrix(position_matrix)) {
    stop(context, " must be stored as a matrix.", call. = FALSE)
  }

  storage <- if (row_major) "profile_row_major" else "shape_column_major"
  axis_names <- if (row_major) rownames(position_matrix) else
    colnames(position_matrix)

  # Prefer named axes from the internal geometry contract when available =======
  if (!is.null(axis_names)) {
    .validate_geometry_contract(
      position_matrix,
      storage = storage,
      context = context
    )

    schema <- .geometry_contract_schema()[[storage]]
    axis_idx <- match(schema[[axis]], axis_names, nomatch = 0)
    axis_idx <- axis_idx[axis_idx > 0]
    if (length(axis_idx) > 0) {
      values <- if (row_major) {
        position_matrix[axis_idx[1], ]
      } else {
        position_matrix[, axis_idx[1]]
      }
      return(as.numeric(values))
    }
  }

  # Fall back to the supplied default or the leading x axis ====================
  if (!is.null(default)) {
    return(as.numeric(default))
  }

  if (axis == "x") {
    return(as.numeric(if (row_major) position_matrix[1, ] else
      position_matrix[, 1]))
  }

  # Return NULL when the requested optional axis is unavailable ================
  NULL
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
  # Resolve the requested storage contract and validate the input matrix =======
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

  # Require named axes before applying the geometry schema checks ==============
  if (is.null(axis_names)) {
    stop(
      context, " is missing ", axis_label,
      " names required by the internal geometry contract.",
      call. = FALSE
    )
  }

  # Require an x axis and paired vertical envelopes when present ===============
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

  # Return invisibly after successful contract validation ======================
  invisible(TRUE)
}

#' Resolve the number of intervals represented by a position matrix
#' @param position_matrix Numeric position matrix.
#' @param row_major Logical; whether nodes are stored in columns.
#' @keywords internal
#' @noRd
.shape_segment_count <- function(position_matrix, row_major = FALSE) {
  # Validate the position matrix before counting its intervals =================
  if (!is.matrix(position_matrix)) {
    stop("'position_matrix' must be a matrix.", call. = FALSE)
  }

  if (row_major) {
    return(ncol(position_matrix) - 1L)
  }

  # Return the interval count for column-major geometry ========================
  nrow(position_matrix) - 1L
}

#' Resolve a user-facing shape label from a shape object or request
#' @param shape_input Shape object.
#' @param requested_shape Optional original shape request.
#' @keywords internal
#' @noRd
.scatterer_shape_name <- function(shape_input, requested_shape = NULL) {
  # Preserve the explicit arbitrary label requested by the user ================
  if (is.character(requested_shape) && requested_shape == "arbitrary") {
    return("Arbitrary")
  }

  # Otherwise use the concrete Shape class name ================================
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
  # Fall back to the supplied shape object when one is available ===============
  if (!is.null(shape_input)) {
    position_matrix <- acousticTS::extract(shape_input, "position_matrix")
    shape_parameters <- acousticTS::extract(shape_input, "shape_parameters")
  }

  # Prefer a scalar radius stored directly on the shape parameters =============
  if (!is.null(shape_parameters$radius) &&
      length(shape_parameters$radius) == 1 &&
      !is.na(shape_parameters$radius)) {
    return(as.numeric(shape_parameters$radius)[1])
  }

  # Otherwise recover the maximum nodewise radius from the profile =============
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
  # Recover the source shape parameters used for class-specific metadata =======
  input_params <- acousticTS::extract(shape_input, "shape_parameters")

  # Append any semiaxis metadata carried by spheroidal shapes ==================
  if (methods::is(shape_input, "ProlateSpheroid")) {
    shape_params$semimajor_length <- input_params$semimajor_length
    shape_params$semiminor_length <- input_params$semiminor_length
    if (!is.null(input_params$radius) && length(input_params$radius) > 1) {
      shape_params$radius_shape <- input_params$radius
    }
  }

  if (methods::is(shape_input, "OblateSpheroid")) {
    shape_params$semimajor_length <- input_params$semimajor_length
    shape_params$semiminor_length <- input_params$semiminor_length
    if (!is.null(input_params$radius) && length(input_params$radius) > 1) {
      shape_params$radius_shape <- input_params$radius
    }
  }

  # Append any cylinder-specific taper or curvature bookkeeping ================
  if (methods::is(shape_input, "Cylinder")) {
    shape_params$radius_curvature_ratio <- input_params$radius_curvature_ratio
    if ("taper_order" %in% names(input_params)) {
      shape_params$taper_order <- input_params$taper_order
    }
  }

  # Preserve any nodewise radius/diameter profiles when available =============
  if ("radius_shape" %in% names(input_params)) {
    shape_params$radius_shape <- input_params$radius_shape
  }
  if ("diameter_shape" %in% names(input_params)) {
    shape_params$diameter_shape <- input_params$diameter_shape
  }

  # Resolve the user-facing shape label and return the updated list ============
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
  # Validate the input shape geometry against the shared contract ==============
  position_matrix <- acousticTS::extract(shape_input, "position_matrix")
  .validate_geometry_contract(
    position_matrix,
    storage = "shape_column_major",
    context = paste0(error_context, " position matrix")
  )

  # Build the base length/radius/segment metadata for the constructor ==========
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

  # Return the full shape-parameter bundle including units =====================
  c(params, extra_units)
}

#' Define the canonical x/w/zU/zL schema for profile-style components
#' @param suffix Component suffix such as `"body"` or `"bladder"`.
#' @keywords internal
#' @noRd
.shape_component_schema <- function(suffix = "body") {
  # Return the canonical row-major schema for one named component ==============
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

#' Validate segmented material inputs against a component segment count
#' @param values Input value supplied for a segmented property.
#' @param expected_segments Expected number of segments.
#' @param argument_name Argument name used in error messages.
#' @keywords internal
#' @noRd
.validate_segmented_property_length <- function(values,
                                                expected_segments,
                                                argument_name) {
  # Ignore NULL and scalar inputs, which are allowed for segmented properties ==
  if (is.null(values) || length(values) <= 1) {
    return(invisible(TRUE))
  }

  # Require vector inputs to match the expected segment count ==================
  if (length(values) != expected_segments) {
    stop(
      sprintf(
        paste0(
          "Vector input for '%s' with %d elements does not match the expected ",
          "number of segments (%d)."
        ),
        argument_name,
        length(values),
        expected_segments
      ),
      call. = FALSE
    )
  }

  # Return invisibly after successful validation ===============================
  invisible(TRUE)
}

#' Validate all segmented material inputs for a fluid-like scatterer body
#' @param shape_obj Shape object used to define the body geometry.
#' @param g Density contrast.
#' @param density Absolute density.
#' @param h Sound-speed contrast.
#' @param sound_speed Absolute sound speed.
#' @keywords internal
#' @noRd
.validate_fls_segmented_material_lengths <- function(shape_obj,
                                                     g = NULL,
                                                     density = NULL,
                                                     h = NULL,
                                                     sound_speed = NULL) {
  # Recover the expected body segment count from the supplied shape ============
  n_segments <- .shape_segment_count(
    acousticTS::extract(shape_obj, "position_matrix")
  )

  .validate_segmented_property_length(g, n_segments, "g_body")
  .validate_segmented_property_length(density, n_segments, "density_body")
  .validate_segmented_property_length(h, n_segments, "h_body")
  .validate_segmented_property_length(
    sound_speed,
    n_segments,
    "sound_speed_body"
  )

  # Return invisibly after successful validation ===============================
  invisible(TRUE)
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
  # Build the canonical row-major profile matrix for one component =============
  schema <- .shape_component_schema(suffix)
  rpos <- rbind(x, w, zU, zL)
  rownames(rpos) <- schema$rows
  # Return the named profile matrix ============================================
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
  # Return the first matching named column from the candidate list =============
  for (nm in candidates) {
    if (!is.null(colnames(position_matrix)) && nm %in%
        colnames(position_matrix)) {
      return(position_matrix[, nm])
    }
  }

  if (!is.null(default)) {
    return(default)
  }

  # Fall back to the leading column when no candidate matched ==================
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
  # Return the first matching named row from the candidate list ================
  for (nm in candidates) {
    if (!is.null(rownames(position_matrix)) && nm %in%
        rownames(position_matrix)) {
      return(position_matrix[nm, ])
    }
  }

  if (!is.null(default)) {
    return(default)
  }

  # Fall back to the leading row when no candidate matched =====================
  position_matrix[1, ]
}

#' Translate a shape position matrix without changing its topology
#' @param position_matrix Column-major position matrix.
#' @param x_offset Along-axis translation.
#' @param y_offset Lateral translation.
#' @param z_offset Vertical translation applied to z/zU/zL coordinates.
#' @keywords internal
#' @noRd
.translate_shape_position_matrix <- function(position_matrix,
                                            x_offset = 0,
                                            y_offset = 0,
                                            z_offset = 0) {
  # Copy the shape matrix before applying any rigid translation ================
  translated <- position_matrix

  # Apply the requested offsets to any recognized coordinate columns ===========
  if (!is.null(colnames(translated))) {
    if ("x" %in% colnames(translated)) translated[, "x"] <-
        translated[, "x"] + x_offset
    if ("y" %in% colnames(translated)) translated[, "y"] <-
        translated[, "y"] + y_offset
    if ("w" %in% colnames(translated)) translated[, "w"] <-
        translated[, "w"] + y_offset
    if ("z" %in% colnames(translated)) translated[, "z"] <-
        translated[, "z"] + z_offset
    if ("zU" %in% colnames(translated)) translated[, "zU"] <-
        translated[, "zU"] + z_offset
    if ("zL" %in% colnames(translated)) translated[, "zL"] <-
        translated[, "zL"] + z_offset
  }

  # Return the translated shape matrix =========================================
  translated
}

#' Convert a shape object into canonical row-major profile coordinates
#' @param shape_obj Shape object.
#' @param suffix Component suffix such as `"body"` or `"bladder"`.
#' @keywords internal
#' @noRd
.as_component_rpos <- function(shape_obj, suffix = "body") {
  # Recover the source shape matrix and a readable validation label ============
  position_matrix <- acousticTS::extract(shape_obj, "position_matrix")
  context <- paste0(
    toupper(substring(suffix, 1, 1)),
    substring(suffix, 2),
    " shape"
  )

  # Resolve the canonical x/w/zU/zL component axes from the shape ==============
  x <- .geometry_axis_values(
    position_matrix,
    axis = "x",
    context = context
  )
  w <- .geometry_axis_values(
    position_matrix,
    axis = "w",
    default = rep(0, length(x)),
    context = context
  )
  zU <- .geometry_axis_values(
    position_matrix,
    axis = "zU",
    default = rep(NA_real_, length(x)),
    context = context
  )
  zL <- .geometry_axis_values(
    position_matrix,
    axis = "zL",
    default = rep(NA_real_, length(x)),
    context = context
  )

  # Return the canonical row-major component profile ===========================
  .canonical_profile_matrix(x, w, zU, zL, suffix = suffix)
}

#' Extract canonical x/w/zU/zL fields from a row-major profile matrix
#' @param position_matrix Row-major profile matrix.
#' @keywords internal
#' @noRd
.profile_fields <- function(position_matrix) {
  # Resolve the canonical axes from a row-major profile matrix =================
  context <- "Profile matrix"
  x <- .geometry_axis_values(
    position_matrix,
    axis = "x",
    row_major = TRUE,
    context = context
  )

  # Return the named field bundle used by profile-based helpers ================
  list(
    x = x,
    w = .geometry_axis_values(
      position_matrix,
      axis = "w",
      row_major = TRUE,
      default = rep(0, length(x)),
      context = context
    ),
    zU = .geometry_axis_values(
      position_matrix,
      axis = "zU",
      row_major = TRUE,
      default = rep(NA_real_, length(x)),
      context = context
    ),
    zL = .geometry_axis_values(
      position_matrix,
      axis = "zL",
      row_major = TRUE,
      default = rep(NA_real_, length(x)),
      context = context
    )
  )
}

#' Validate that a component profile contains at least one interval
#' @param rpos Row-major position matrix.
#' @param component_name Human-readable component label.
#' @keywords internal
#' @noRd
.validate_component_rpos <- function(rpos, component_name = "Component") {
  # Require at least one interval in the component profile =====================
  if (ncol(rpos) < 2) {
    stop(
      component_name,
      " shape must have at least two points (>=1 segment).",
      call. = FALSE
    )
  }

  # Return the profile invisibly after validation ==============================
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
  # Convert the shape into the canonical row-major component profile ===========
  rpos <- .as_component_rpos(shape_obj, suffix = suffix)
  .validate_component_rpos(
    rpos,
    component_name = paste0(
      toupper(substring(suffix, 1, 1)),
      substring(suffix, 2)
    )
  )

  # Assemble the fluid-like component fields and any extra metadata ============
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
  # Recover and validate the source shape geometry =============================
  position_matrix <- acousticTS::extract(shape_obj, "position_matrix")
  shape_parameters <- acousticTS::extract(shape_obj, "shape_parameters")

  .validate_geometry_contract(
    position_matrix,
    storage = "shape_column_major",
    context = "Shape position matrix"
  )

  # Fall back to the shape radius metadata when no override is supplied ========
  if (is.null(radius)) {
    radius <- shape_parameters$radius
  }

  # Assemble the row-major fluid component payload =============================
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
  # Recover and validate the source shape geometry =============================
  position_matrix <- acousticTS::extract(shape_obj, "position_matrix")
  shape_parameters <- acousticTS::extract(shape_obj, "shape_parameters")

  .validate_geometry_contract(
    position_matrix,
    storage = "shape_column_major",
    context = "Shape position matrix"
  )

  # Fall back to the shape radius metadata when no override is supplied ========
  if (is.null(radius)) {
    radius <- shape_parameters$radius
  }

  # Assemble the column-major fluid component payload ==========================
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
  # Recover the active shape geometry, honoring any explicit override ==========
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

  # Fall back to the shape radius metadata when no override is supplied ========
  if (is.null(radius)) {
    radius <- shape_parameters$radius
  }

  # Assemble the row-major elastic component payload ===========================
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

#' Complete the shell elastic-modulus set from any two supplied values
#' @param E Young's modulus.
#' @param G Shear modulus.
#' @param K Bulk modulus.
#' @param nu Poisson's ratio.
#' @keywords internal
#' @noRd
.complete_elastic_moduli <- function(E = NULL,
                                     G = NULL,
                                     K = NULL,
                                     nu = NULL) {
  # Collect the elastic inputs that were actually supplied =====================
  elastic_params <- Filter(
    Negate(is.null),
    list(E = E, nu = nu, G = G, K = K)
  )

  # Iteratively fill in the missing elastic constants from the known ones ======
  repeat {
    param_len <- length(elastic_params)
    if (is.null(elastic_params$nu)) {
      elastic_params$nu <- tryCatch(
        pois(elastic_params$K, elastic_params$E, elastic_params$G),
        error = function(e) NULL
      )
    }
    if (is.null(elastic_params$K)) {
      elastic_params$K <- tryCatch(
        bulk(elastic_params$E, elastic_params$G, elastic_params$nu),
        error = function(e) NULL
      )
    }
    if (is.null(elastic_params$E)) {
      elastic_params$E <- tryCatch(
        young(elastic_params$K, elastic_params$G, elastic_params$nu),
        error = function(e) NULL
      )
    }
    if (is.null(elastic_params$G)) {
      elastic_params$G <- tryCatch(
        shear(elastic_params$K, elastic_params$E, elastic_params$nu),
        error = function(e) NULL
      )
    }
    if (is.null(elastic_params$lambda)) {
      elastic_params$lambda <- tryCatch(
        lame(
          elastic_params$K, elastic_params$E, elastic_params$G,
          elastic_params$nu
        ),
        error = function(e) NULL
      )
    }

    if (length(elastic_params) == param_len) {
      break
    }
  }

  # Return the completed elastic-modulus set ===================================
  elastic_params
}

#' Normalize shape radius/diameter metadata for ESS-style reporting
#' @param params Shape-parameter list.
#' @keywords internal
#' @noRd
.normalize_shape_diameter <- function(params) {
  # Collapse vector radius metadata while preserving the original profile ======
  if (!is.null(params$radius) && length(params$radius) > 1 &&
      !all(is.na(params$radius))) {
    if (is.null(params$radius_shape)) {
      params$radius_shape <- params$radius
    }
    params$radius <- max(params$radius, na.rm = TRUE)
  }

  # Collapse vector diameter metadata while preserving the original profile ====
  if (!is.null(params$diameter) && length(params$diameter) > 1 &&
      !all(is.na(params$diameter))) {
    if (is.null(params$diameter_shape)) {
      params$diameter_shape <- params$diameter
    }
    params$diameter <- max(params$diameter, na.rm = TRUE)
  }

  # Reconstruct a scalar diameter when only radius-style metadata exists =======
  if (is.null(params$diameter)) {
    if (!is.null(params$diameter_shape)) {
      params$diameter <- if (length(params$diameter_shape) == 1) {
        params$diameter_shape
      } else {
        max(params$diameter_shape, na.rm = TRUE)
      }
    } else if (!is.null(params$radius)) {
      params$diameter <- params$radius * 2
    }
  }

  # Return the normalized reporting metadata ===================================
  params
}

#' Derive an approximate internal-fluid shape from an outer shell shape
#' @param shape_obj Outer shell Shape object.
#' @param shell_thickness Shell thickness in the same length units as the shape.
#' @keywords internal
#' @noRd
.inner_shape_from_shell <- function(shape_obj, shell_thickness) {
  # Exit early when no internal-fluid reconstruction was requested =============
  if (is.null(shell_thickness)) {
    return(NULL)
  }
  if (!is.numeric(shell_thickness) || length(shell_thickness) != 1 ||
      !is.finite(shell_thickness) || shell_thickness <= 0) {
    stop("'shell_thickness' must be a single positive number.", call. = FALSE)
  }

  # Recover the outer shell geometry and validate the thickness ================
  position_matrix <- acousticTS::extract(shape_obj, "position_matrix")
  shape_parameters <- acousticTS::extract(shape_obj, "shape_parameters")
  max_radius <- .shape_max_radius(shape_input = shape_obj,
                                  error_context = "ESS shell")
  if (!is.finite(max_radius) || shell_thickness >= max_radius) {
    stop(
      "'shell_thickness' must be smaller than the shell radius.",
      call. = FALSE
    )
  }

  # Use analytic inner-shape reconstructions for the canonical shape families ==
  class_name <- class(shape_obj)[1]
  if (identical(class_name, "Sphere")) {
    return(sphere(
      radius_body = shape_parameters$radius - shell_thickness,
      n_segments = shape_parameters$n_segments,
      diameter_units = shape_parameters$diameter_units %||% "m"
    ))
  }
  if (identical(class_name, "Cylinder")) {
    return(cylinder(
      length_body = shape_parameters$length,
      radius_body = max_radius - shell_thickness,
      taper = if (is.na(shape_parameters$taper_order)) NULL else {
        shape_parameters$taper_order
      },
      radius_curvature_ratio = if (is.na(
        shape_parameters$radius_curvature_ratio
      )) {
        NULL
      } else {
        shape_parameters$radius_curvature_ratio
      },
      n_segments = shape_parameters$n_segments,
      length_units = shape_parameters$length_units %||% "m"
    ))
  }
  if (identical(class_name, "ProlateSpheroid")) {
    return(prolate_spheroid(
      length_body = shape_parameters$length,
      radius_body = shape_parameters$semiminor_length - shell_thickness,
      n_segments = shape_parameters$n_segments,
      length_units = shape_parameters$length_units %||% "m"
    ))
  }
  if (identical(class_name, "OblateSpheroid")) {
    return(oblate_spheroid(
      length_body = shape_parameters$length,
      radius_body = shape_parameters$semimajor_length - shell_thickness,
      n_segments = shape_parameters$n_segments,
      length_units = shape_parameters$length_units %||% "m"
    ))
  }

  # Fall back to shrinking the explicit arbitrary profile node by node =========
  x <- .geometry_axis_values(
    position_matrix,
    axis = "x",
    context = "ESS shell position matrix"
  )
  radius_profile <- pmax(
    .shape_radius_profile(
      position_matrix = position_matrix,
      shape_parameters = shape_parameters,
      error_context = "ESS shell"
    ) - shell_thickness,
    0
  )
  width <- .geometry_axis_values(
    position_matrix,
    axis = "w",
    default = 2 * radius_profile,
    context = "ESS shell position matrix"
  )
  zU <- .geometry_axis_values(
    position_matrix,
    axis = "zU",
    default = radius_profile,
    context = "ESS shell position matrix"
  )
  zL <- .geometry_axis_values(
    position_matrix,
    axis = "zL",
    default = -radius_profile,
    context = "ESS shell position matrix"
  )
  z_center <- if (!is.null(colnames(position_matrix)) &&
      any(c("z", "z_body", "z_shell", "z_fluid") %in%
          colnames(position_matrix))) {
    .extract_shape_component_column(
      position_matrix,
      c("z", "z_body", "z_shell", "z_fluid"),
      default = (zU + zL) / 2
    )
  } else {
    (zU + zL) / 2
  }

  # Return the reconstructed inner arbitrary shape =============================
  arbitrary(
    x_body = x,
    w_body = pmax(width - 2 * shell_thickness, 0),
    zU_body = z_center + radius_profile,
    zL_body = z_center - radius_profile,
    radius_body = radius_profile,
    length_units = shape_parameters$length_units %||%
      shape_parameters$diameter_units %||%
      "m"
  )
}

#' Resolve ESS shell/fluid shapes from either explicit shapes or legacy inputs
#' @param shape Shape object or legacy shape identifier.
#' @param arguments Evaluated constructor arguments.
#' @param radius_shell Optional shell radius used by legacy canonical dispatch.
#' @param shell_thickness Optional shell thickness.
#' @param x_body Optional x coordinates for arbitrary legacy shells.
#' @param y_body Optional y coordinates for arbitrary legacy shells.
#' @param z_body Optional z coordinates for arbitrary legacy shells.
#' @keywords internal
#' @noRd
.resolve_ess_shape_components <- function(shape,
                                          arguments,
                                          radius_shell = NULL,
                                          shell_thickness = NULL,
                                          x_body = NULL,
                                          y_body = NULL,
                                          z_body = NULL) {
  # Resolve the outer shell shape from explicit or legacy inputs ===============
  shell_shape <- if (methods::is(shape, "Shape")) {
    shape
  } else if (is.null(shape) && !is.null(x_body)) {
    arbitrary(
      x_body = x_body,
      y_body = y_body,
      z_body = z_body,
      radius_body = radius_shell
    )
  } else if (identical(shape, "arbitrary")) {
    arbitrary(
      x_body = x_body,
      y_body = y_body,
      z_body = z_body,
      radius_body = radius_shell
    )
  } else if (is.null(shape)) {
    stop(
      "Supply 'shape' as a pre-built Shape object for `ess_generate()`, or ",
      "provide explicit shell profile coordinates such as ",
      "x_body/y_body/z_body.",
      call. = FALSE
    )
  } else {
    shell_call <- arguments
    shell_call$radius_body <- radius_shell
    .resolve_scatterer_shape_input(shape, shell_call)
  }

  # Return the shell and its derived internal-fluid surrogate ==================
  list(
    shell = shell_shape,
    fluid = .inner_shape_from_shell(shell_shape, shell_thickness)
  )
}
