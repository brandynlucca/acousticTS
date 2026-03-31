################################################################################
# FORGE FUNCTIONS FOR MANIPULATING SCATTERER SHAPES
################################################################################
################################################################################
# PRIMARY FORGE GENERATION FUNCTION
################################################################################
#' Resize or reparameterize a scatterer object
#'
#' @description
#' Generic function for rescaling or otherwise reparameterizing an existing
#' scatterer while preserving its class semantics. The class-specific methods
#' handle the component bookkeeping needed to keep the resulting object
#' structurally valid after lengths, widths, heights, shell thicknesses, or
#' related geometry descriptors are changed.
#'
#' @param object A scatterer object.
#' @param ... Additional arguments passed to specific methods.
#'
#' @details
#' `reforge()` is the package's main post-construction geometry-adjustment
#' tool. It is useful when a target should be modified in place rather than
#' rebuilt from scratch. The available method-specific arguments depend on the
#' scatterer class and typically include direct scale factors and/or target
#' dimensions for length, width, or height.
#'
#' @return A scatterer of the same broad class as `object`, rebuilt with the
#'   requested geometric changes.
#'
#' @seealso [extract()], [plot.Scatterer()], [brake()]
#'
#' @examples
#' obj <- fls_generate(
#'   shape = sphere(radius_body = 0.01, n_segments = 40),
#'   density_body = 1045,
#'   sound_speed_body = 1520
#' )
#'
#' bigger_obj <- reforge(obj, body_target = c(length = 0.03))
#' extract(bigger_obj, c("shape_parameters", "length"))
#' @export
setGeneric(
  "reforge",
  function(object, ...) {
    object
  }
)

#' Resolve scaling vector from direct scale or target dimensions
#' @param scale Optional direct scaling input.
#' @param target Optional named target-dimension input.
#' @param dims Named vector of current component dimensions.
#' @keywords internal
#' @noRd
.reforge_scale_vector <- function(scale = NULL, target = NULL, dims) {
  # Preserve the direct scaling pathway when no targets are supplied ===========
  if (is.null(target)) {
    return(list(suffix = "_scale", scale = scale))
  }

  # Convert named target dimensions into multiplicative scale factors ==========
  out <- c()
  for (nm in names(target)) {
    out[nm] <- target[nm] / dims[nm]
  }

  # Return the normalized scaling specification ================================
  list(suffix = "_target", scale = out)
}

#' Derive length/width/height from a row-major component position matrix
#' @param rpos Row-major position matrix.
#' @param length Optional externally supplied length override.
#' @keywords internal
#' @noRd
.reforge_component_dimensions <- function(rpos, length = NULL) {
  # Recover the width profile needed for scalar component summaries ============
  width_profile <- .shape_width_profile(
    position_matrix = rpos,
    row_major = TRUE
  )
  # Return length, width, and height for the current component =================
  c(
    length = if (!is.null(length)) {
      length
    } else {
      .shape_length(position_matrix = rpos, row_major = TRUE)
    },
    width = max(abs(width_profile), na.rm = TRUE),
    height = max(
      .shape_height_profile(position_matrix = rpos, row_major = TRUE),
      na.rm = TRUE
    )
  )
}

#' Resolve the first matching profile row index from candidate names
#' @param rpos Row-major position matrix.
#' @param candidates Candidate row names.
#' @keywords internal
#' @noRd
.reforge_profile_row_idx <- function(rpos, candidates) {
  # Return the first matching row index for the requested profile field ========
  idx <- match(candidates, rownames(rpos), nomatch = 0)
  idx <- idx[idx > 0]
  if (length(idx) == 0) {
    integer(0)
  } else {
    idx[1]
  }
}

#' Compute the local profile centerline from row-major geometry fields
#' @param fields Output of \code{.profile_fields()}.
#' @keywords internal
#' @noRd
.reforge_profile_centerline <- function(fields) {
  # Use the envelope midpoint whenever upper/lower profiles are available ======
  if (!is.null(fields$zU) && !is.null(fields$zL) &&
    !all(is.na(fields$zU)) && !all(is.na(fields$zL))) {
    return((fields$zU + fields$zL) / 2)
  }

  # Fall back to a flat centerline when no explicit envelopes exist ============
  rep(0, length(fields$x))
}

#' Compute the local half-height from row-major geometry fields
#' @param fields Output of \code{.profile_fields()}.
#' @keywords internal
#' @noRd
.reforge_profile_half_height <- function(fields) {
  # Recover half-height from the stored vertical envelopes =====================
  if (!is.null(fields$zU) && !is.null(fields$zL) &&
    !all(is.na(fields$zU)) && !all(is.na(fields$zL))) {
    return((fields$zU - fields$zL) / 2)
  }

  # Fall back to zeros when no vertical envelope is available ==================
  rep(0, length(fields$x))
}

#' Apply named axis scalings to a row-major component position matrix
#' @param rpos Row-major position matrix.
#' @param scales Named vector with length/width/height multipliers.
#' @keywords internal
#' @noRd
.reforge_apply_axis_scaling <- function(rpos, scales, scale_centerline = FALSE) {
  # Leave the component untouched when no scaling was requested ================
  if (is.null(scales)) {
    return(rpos)
  }

  # Validate the row-major component matrix before editing it ==================
  .validate_geometry_contract(
    rpos,
    storage = "profile_row_major",
    context = "Reforge component matrix"
  )

  # Resolve the centerline and envelope rows used by the scaling update =======
  z_idx <- .reforge_profile_row_idx(
    rpos,
    c("z", "z_body", "z_bladder", "z_shell", "z_fluid", "z_backbone")
  )
  zU_idx <- .reforge_profile_row_idx(
    rpos,
    .geometry_contract_schema()$profile_row_major$zU
  )
  zL_idx <- .reforge_profile_row_idx(
    rpos,
    .geometry_contract_schema()$profile_row_major$zL
  )

  # Rescale the axial coordinate relative to its leading anchor point ==========
  if (scales["length"] != 1) {
    x_idx <- .reforge_profile_row_idx(
      rpos,
      .geometry_contract_schema()$profile_row_major$x
    )
    anchor <- rpos[x_idx, 1]
    rpos[x_idx, ] <- anchor + (rpos[x_idx, ] - anchor) * scales["length"]

    # Rescale the stored centerline path when curved profiles request it =======
    if (isTRUE(scale_centerline) && length(z_idx) > 0) {
      z_anchor <- rpos[z_idx, 1]
      z_old <- rpos[z_idx, ]
      z_new <- z_anchor + (z_old - z_anchor) * scales["length"]
      z_delta <- z_new - z_old
      rpos[z_idx, ] <- z_new

      if (length(zU_idx) > 0) {
        rpos[zU_idx, ] <- rpos[zU_idx, ] + z_delta
      }
      if (length(zL_idx) > 0) {
        rpos[zL_idx, ] <- rpos[zL_idx, ] + z_delta
      }
    }
  }

  # Rescale the lateral width profile when present =============================
  if (scales["width"] != 1) {
    w_idx <- .reforge_profile_row_idx(
      rpos,
      .geometry_contract_schema()$profile_row_major$w
    )
    if (length(w_idx) > 0) {
      rpos[w_idx, ] <- rpos[w_idx, ] * scales["width"]
    }
  }

  # Rescale the vertical profile while preserving its local centerline =========
  if (scales["height"] != 1) {
    if (length(zU_idx) > 0 && length(zL_idx) > 0) {
      z_center <- if (length(z_idx) > 0) {
        rpos[z_idx, ]
      } else {
        (rpos[zU_idx, ] + rpos[zL_idx, ]) / 2
      }
      half_height <- (rpos[zU_idx, ] - rpos[zL_idx, ]) / 2
      half_height <- half_height * scales["height"]
      rpos[zU_idx, ] <- z_center + half_height
      rpos[zL_idx, ] <- z_center - half_height
      if (length(z_idx) > 0) {
        rpos[z_idx, ] <- z_center
      }
    } else if (length(z_idx) > 0) {
      rpos[z_idx, ] <- rpos[z_idx, ] * scales["height"]
    }
  }

  # Return the scaled component matrix =========================================
  rpos
}

#' Resample a row-major component position matrix to a new point count
#' @param rpos Row-major position matrix.
#' @param n_points Target number of points.
#' @keywords internal
#' @noRd
.reforge_resample_rows <- function(rpos, n_points) {
  # Preserve the original component when no resampling was requested ===========
  if (is.null(n_points)) {
    return(rpos)
  }

  # Resample the row-major profile to the requested node count =================
  .resample_rpos_rows(rpos, as.integer(n_points))
}

#' Express an internal component start position relative to its parent length
#' @param component_rpos Row-major internal-component position matrix.
#' @param parent_rpos Row-major parent-component position matrix.
#' @keywords internal
#' @noRd
.reforge_relative_start <- function(component_rpos, parent_rpos) {
  # Validate the internal and parent component matrices ========================
  .validate_geometry_contract(
    component_rpos,
    storage = "profile_row_major",
    context = "Internal component matrix"
  )
  .validate_geometry_contract(
    parent_rpos,
    storage = "profile_row_major",
    context = "Parent component matrix"
  )
  # Express the internal start point relative to the parent length =============
  component_x <- .shape_x(component_rpos, row_major = TRUE)
  parent_x <- .shape_x(parent_rpos, row_major = TRUE)
  (component_x[1] - min(parent_x, na.rm = TRUE)) /
    .shape_length(position_matrix = parent_rpos, row_major = TRUE)
}

#' Reposition an internal component to preserve its relative body start
#' @param component_rpos Row-major internal-component position matrix.
#' @param relative_start Relative start position along the parent length.
#' @param parent_rpos Row-major parent-component position matrix.
#' @keywords internal
#' @noRd
.reforge_shift_to_relative_start <- function(component_rpos,
                                             relative_start,
                                             parent_rpos) {
  # Validate the internal and parent component matrices ========================
  .validate_geometry_contract(
    component_rpos,
    storage = "profile_row_major",
    context = "Internal component matrix"
  )
  .validate_geometry_contract(
    parent_rpos,
    storage = "profile_row_major",
    context = "Parent component matrix"
  )
  # Convert the stored relative start back into an absolute x offset ===========
  parent_x <- .shape_x(parent_rpos, row_major = TRUE)
  component_x <- .shape_x(component_rpos, row_major = TRUE)
  new_start <- min(parent_x, na.rm = TRUE) +
    relative_start * .shape_length(
      position_matrix = parent_rpos,
      row_major = TRUE
    )
  shift_amount <- new_start - component_x[1]
  x_idx <- .reforge_profile_row_idx(
    component_rpos,
    .geometry_contract_schema()$profile_row_major$x
  )
  component_rpos[x_idx, ] <- component_rpos[x_idx, ] + shift_amount
  # Return the shifted internal component ======================================
  component_rpos
}

#' Express an internal component centerline offset relative to parent thickness
#' @param component_rpos Row-major internal-component position matrix.
#' @param parent_rpos Row-major parent-component position matrix.
#' @keywords internal
#' @noRd
.reforge_relative_vertical_offset <- function(component_rpos, parent_rpos) {
  # Validate the internal and parent component matrices ========================
  .validate_geometry_contract(
    component_rpos,
    storage = "profile_row_major",
    context = "Internal component matrix"
  )
  .validate_geometry_contract(
    parent_rpos,
    storage = "profile_row_major",
    context = "Parent component matrix"
  )

  # Interpolate the parent centerline/thickness onto the component x-grid ======
  component_fields <- .profile_fields(component_rpos)
  parent_fields <- .profile_fields(parent_rpos)
  interp <- function(x, y) {
    stats::approx(x, y, xout = component_fields$x, rule = 2)$y
  }

  parent_center <- interp(
    parent_fields$x,
    .reforge_profile_centerline(parent_fields)
  )
  parent_half_height <- interp(
    parent_fields$x,
    .reforge_profile_half_height(parent_fields)
  )
  component_center <- .reforge_profile_centerline(component_fields)

  # Store a robust median offset ratio so it can be reconstructed later =======
  valid <- is.finite(parent_half_height) &
    (abs(parent_half_height) > sqrt(.Machine$double.eps))
  if (!any(valid)) {
    return(0)
  }

  stats::median(
    (component_center[valid] - parent_center[valid]) /
      parent_half_height[valid],
    na.rm = TRUE
  )
}

#' Reposition an internal component to preserve relative vertical placement
#' @param component_rpos Row-major internal-component position matrix.
#' @param relative_offset Relative centerline offset scaled by parent half-height.
#' @param parent_rpos Row-major parent-component position matrix.
#' @keywords internal
#' @noRd
.reforge_shift_to_relative_vertical_offset <- function(component_rpos,
                                                       relative_offset,
                                                       parent_rpos) {
  # Validate the internal and parent component matrices ========================
  .validate_geometry_contract(
    component_rpos,
    storage = "profile_row_major",
    context = "Internal component matrix"
  )
  .validate_geometry_contract(
    parent_rpos,
    storage = "profile_row_major",
    context = "Parent component matrix"
  )

  # Skip no-op calls cleanly ==================================================
  if (is.null(relative_offset) || !is.finite(relative_offset)) {
    return(component_rpos)
  }

  # Interpolate the parent centerline/thickness onto the component x-grid ======
  component_fields <- .profile_fields(component_rpos)
  parent_fields <- .profile_fields(parent_rpos)
  interp <- function(x, y) {
    stats::approx(x, y, xout = component_fields$x, rule = 2)$y
  }

  parent_center <- interp(
    parent_fields$x,
    .reforge_profile_centerline(parent_fields)
  )
  parent_half_height <- interp(
    parent_fields$x,
    .reforge_profile_half_height(parent_fields)
  )
  current_center <- .reforge_profile_centerline(component_fields)
  target_center <- parent_center + relative_offset * parent_half_height
  delta <- target_center - current_center

  # Shift every vertical coordinate row by the reconstructed centerline delta ==
  z_idx <- .reforge_profile_row_idx(component_rpos, c("z", "z_body", "z_bladder"))
  zU_idx <- .reforge_profile_row_idx(
    component_rpos,
    .geometry_contract_schema()$profile_row_major$zU
  )
  zL_idx <- .reforge_profile_row_idx(
    component_rpos,
    .geometry_contract_schema()$profile_row_major$zL
  )

  if (length(z_idx) > 0) {
    component_rpos[z_idx, ] <- component_rpos[z_idx, ] + delta
  }
  if (length(zU_idx) > 0) {
    component_rpos[zU_idx, ] <- component_rpos[zU_idx, ] + delta
  }
  if (length(zL_idx) > 0) {
    component_rpos[zL_idx, ] <- component_rpos[zL_idx, ] + delta
  }

  # Return the shifted internal component =====================================
  component_rpos
}

#' Validate whether an internal component exceeds its containing body bounds
#' @param rpos_b Row-major body position matrix.
#' @param rpos_i Row-major internal-component position matrix.
#' @param component_label Human-readable label for the internal component.
#' @param action Containment action: warn, error, or ignore.
#' @keywords internal
#' @noRd
.reforge_check_internal_containment <- function(rpos_b,
                                                rpos_i,
                                                component_label = "Swimbladder",
                                                action = c("warn", "error", "ignore")) {
  # Normalize the requested containment policy ================================
  action <- match.arg(action)

  # Validate the body and internal profile matrices ============================
  .validate_geometry_contract(
    rpos_b,
    storage = "profile_row_major",
    context = "Body profile matrix"
  )
  .validate_geometry_contract(
    rpos_i,
    storage = "profile_row_major",
    context = "Internal profile matrix"
  )

  # Interpolate both profiles onto a shared axial grid =========================
  body_fields <- .profile_fields(rpos_b)
  inner_fields <- .profile_fields(rpos_i)
  x_grid <- seq(
    max(min(body_fields$x), min(inner_fields$x)),
    min(max(body_fields$x), max(inner_fields$x)),
    length.out = 200
  )

  interp <- function(x, y) stats::approx(x, y, xout = x_grid, rule = 2)$y

  # Compare the interpolated envelopes along the shared grid ===================
  body_y <- interp(body_fields$x, body_fields$w)
  body_zU <- interp(body_fields$x, body_fields$zU)
  body_zL <- interp(body_fields$x, body_fields$zL)

  inner_y <- interp(inner_fields$x, inner_fields$w)
  inner_zU <- interp(inner_fields$x, inner_fields$zU)
  inner_zL <- interp(inner_fields$x, inner_fields$zL)

  contained <- (inner_y <= body_y) & (inner_y >= -body_y) &
    (inner_zU <= body_zU) & (inner_zL >= body_zL)

  # Apply the configured containment policy ===================================
  if (!all(contained)) {
    message_text <- paste0(component_label, " exceeds body bounds at some positions.")
    if (identical(action, "error")) {
      stop(message_text, call. = FALSE)
    }
    if (identical(action, "warn")) {
      warning(message_text, call. = FALSE)
    }
  }

  invisible(all(contained))
}

#' Resizing function for swimbladdered targets
#' @param object SBF-class object.
#' @param body_scale Proportional scaling to the body length, width, and height
#' dimensions. When a single value is supplied, all dimensions are scaled using
#' the same scaling factor. Otherwise, this input must be a named numeric
#' vector.
#' @param body_target Target dimensions (m) for the body length, width, and
#' height dimensions. This input must be a named numeric vector.
#' @param swimbladder_scale Proportional scaling to the swimbladder length,
#' width, and height dimensions. When a single value is supplied, all
#' dimensions are scaled using the same scaling factor. Otherwise, this input
#' must be a named numeric vector.
#' @param swimbladder_target Target dimensions (m) for the swimbladder length,
#' width, and height dimensions. This input must be a named numeric vector.
#' @param swimbladder_inflation_factor Proportional swimbladder volume where
#' the swimbladder x-axis origin and terminus are both held constant.
#' @param maintain_ratio Maintain size ratio between body and
#' swimbladder.
#' @param isometric_body Logical; maintain isometric scaling for body.
#' @param isometric_swimbladder Logical; maintain isometric scaling for bladder.
#' @param n_segments_body Number of segments along the body.
#' @param n_segments_swimbladder Number of segments along the bladder.
#' @param containment Containment policy for internal geometry checks. Use
#'   `"warn"` to keep the current warning behavior, `"error"` to fail fast for
#'   invalid internal geometries, or `"ignore"` to skip containment checks.
#' @keywords internal
#' @export
setMethod(
  "reforge",
  signature(object = "SBF"),
  function(
    object,
    body_scale = NULL,
    body_target = NULL,
    swimbladder_scale = NULL,
    swimbladder_target = NULL,
    isometric_body = TRUE,
    isometric_swimbladder = TRUE,
    maintain_ratio = TRUE,
    swimbladder_inflation_factor = 1.0,
    n_segments_body = NULL,
    n_segments_swimbladder = NULL,
    containment = c("warn", "error", "ignore")
  ) {
    ############################################################################
    # Validation ===============================================================
    containment <- match.arg(containment)
    if (is.null(body_scale) && is.null(swimbladder_scale) &&
      is.null(body_target) && is.null(swimbladder_target) &&
      is.null(n_segments_body) && is.null(n_segments_swimbladder) &&
      swimbladder_inflation_factor == 1.0) {
      stop(
        "Must specify at least one scaling, target, inflation factor, ",
        "or segment count parameter."
      )
    }
    if (!is.null(body_scale) && !is.null(body_target)) {
      stop("Specify only one of body_scale or body_target, not both.")
    }
    if (!is.null(swimbladder_scale) && !is.null(swimbladder_target)) {
      stop(
        "Specify only one of swimbladder_scale or swimbladder_target, not both."
      )
    }
    if ((!is.null(body_scale) || !is.null(body_target)) &&
      (!is.null(swimbladder_scale) || !is.null(swimbladder_target)) &&
      maintain_ratio) {
      maintain_ratio <- FALSE
      message(
        "Hidden State Change Warning:\n ",
        "Multiple axes specified for the body and swimbladder: ",
        "'maintain_ratio' will be ignored for those axes."
      )
    }
    ############################################################################
    # Extract components =======================================================
    body <- acousticTS::extract(object, "body")
    bladder <- acousticTS::extract(object, "bladder")
    shape <- acousticTS::extract(object, "shape_parameters")
    rpos_b <- body$rpos
    rpos_sb <- bladder$rpos
    ############################################################################
    # Calculate current dimensions =============================================
    body_dims <- .reforge_component_dimensions(
      rpos_b,
      length = shape$body$length
    )
    bladder_dims <- .reforge_component_dimensions(
      rpos_sb,
      length = shape$bladder$length %||%
        .shape_length(position_matrix = rpos_sb, row_major = TRUE)
    )
    ############################################################################
    # Calculate swimbladder origin relative to body position ===================
    bladder_relative_start <- .reforge_relative_start(rpos_sb, rpos_b)
    bladder_relative_vertical <- .reforge_relative_vertical_offset(
      rpos_sb,
      rpos_b
    )
    ############################################################################
    # Process target parameters ================================================
    # body_target <- validate_target(body_target, "body_target")
    body_target <- .validate_dimensions_target(
      body_target,
      "body_target",
      c("length", "width", "height")
    )
    body_scale_lst <- .reforge_scale_vector(body_scale, body_target, body_dims)
    swimbladder_target <- .validate_dimensions_target(
      swimbladder_target,
      "swimbladder_target",
      c("length", "width", "height")
    )
    swimbladder_scale_lst <- .reforge_scale_vector(
      swimbladder_scale,
      swimbladder_target,
      bladder_dims
    )
    ############################################################################
    # Process scaling parameters ===============================================
    body_scales <- .validate_dimension_scaling(
      dims = body_scale_lst$scale,
      dims_name = paste0("body", body_scale_lst$suffix),
      valid_dims = c("length", "width", "height"),
      isometry = isometric_body,
      iso_name = "isometric_body"
    )
    bladder_scales <- .validate_dimension_scaling(
      dims = swimbladder_scale_lst$scale,
      dims_name = paste0("swimbladder", swimbladder_scale_lst$suffix),
      valid_dims = c("length", "width", "height"),
      isometry = isometric_swimbladder,
      iso_name = "isometric_swimbladder"
    )
    ############################################################################
    # Apply ratio maintenance logic ============================================
    if (maintain_ratio) {
      if (!is.null(body_scales) && is.null(bladder_scales)) {
        bladder_scales <- body_scales
      } else if (is.null(body_scales) && !is.null(bladder_scales)) {
        body_scales <- bladder_scales
      }
    }
    ############################################################################
    # Interpolate segments first (before scaling) ==============================
    rpos_b <- .reforge_resample_rows(rpos_b, n_segments_body)
    rpos_sb <- .reforge_resample_rows(rpos_sb, n_segments_swimbladder)
    ############################################################################
    # Apply scaling ============================================================
    rpos_b <- .reforge_apply_axis_scaling(rpos_b, body_scales)
    rpos_sb <- .reforge_apply_axis_scaling(rpos_sb, bladder_scales)
    ############################################################################
    # Adjust swimbladder position within scaled body if needed =================
    if (!is.null(body_scales) && body_scales["length"] != 1) {
      rpos_sb <- .reforge_shift_to_relative_start(
        rpos_sb,
        bladder_relative_start,
        rpos_b
      )
    }
    rpos_sb <- .reforge_shift_to_relative_vertical_offset(
      rpos_sb,
      bladder_relative_vertical,
      rpos_b
    )
    ############################################################################
    # Apply bladder inflation factor ===========================================
    if (swimbladder_inflation_factor != 1.0) {
      # Preserve relative position to body
      x_bladder_origin <- bladder$rpos[1, 1] / max(body$rpos[1, ])
      xsb_start <- x_bladder_origin * max(rpos_b[1, ])
      xsb_offset <- rpos_sb[1, 1] - xsb_start

      rpos_sb[1, ] <- rpos_sb[1, ] - xsb_offset
      rpos_sb[2, ] <- rpos_sb[2, ] * swimbladder_inflation_factor
      rpos_sb[3, ] <- rpos_sb[3, ] * swimbladder_inflation_factor
      if (nrow(rpos_sb) >= 4) {
        rpos_sb[4, ] <- rpos_sb[4, ] *
          swimbladder_inflation_factor
      }
      rpos_sb <- .reforge_shift_to_relative_vertical_offset(
        rpos_sb,
        bladder_relative_vertical,
        rpos_b
      )
    }
    ############################################################################
    # Validate swimbladder containment =========================================
    .reforge_check_internal_containment(
      rpos_b,
      rpos_sb,
      action = containment
    )
    ############################################################################
    # Update object ============================================================
    methods::slot(object, "body")$rpos <- rpos_b
    methods::slot(object, "bladder")$rpos <- rpos_sb
    methods::slot(object, "shape_parameters")$body$length <- .shape_length(
      position_matrix = rpos_b,
      row_major = TRUE
    )
    methods::slot(object, "shape_parameters")$bladder$length <- .shape_length(
      position_matrix = rpos_sb,
      row_major = TRUE
    )
    methods::slot(object, "shape_parameters")$body$n_segments <- ncol(rpos_b)
    methods::slot(object, "shape_parameters")$bladder$n_segments <-
      ncol(rpos_sb)
    return(object)
  }
)
################################################################################
#' Reforge BBF-class object.
#'
#' Resize a backboned-fish scatterer by rescaling the flesh body and elastic
#' backbone independently or together. The body follows the same profile-based
#' length/width/height scaling used by `SBF`, while the backbone retains the
#' same cylinder-style component bookkeeping and preserves its relative axial
#' start within the body when body length changes.
#'
#' @param object BBF-class object.
#' @param body_scale Proportional scaling to the body length, width, and height
#'   dimensions. When a single value is supplied, all dimensions are scaled
#'   uniformly. Otherwise, this input must be a named numeric vector.
#' @param body_target Target dimensions (m) for the body length, width, and
#'   height dimensions. This input must be a named numeric vector.
#' @param backbone_scale Proportional scaling to the backbone length, width, and
#'   height dimensions. When a single value is supplied, all dimensions are
#'   scaled uniformly. Otherwise, this input must be a named numeric vector.
#' @param backbone_target Target dimensions (m) for the backbone length, width,
#'   and height dimensions. This input must be a named numeric vector.
#' @param isometric_body Logical; maintain isometric scaling for body.
#' @param isometric_backbone Logical; maintain isometric scaling for backbone.
#' @param maintain_ratio Maintain size ratio between body and backbone.
#' @param n_segments_body Number of points along the body profile.
#' @param n_segments_backbone Number of points along the backbone profile.
#' @param containment Containment policy for internal geometry checks. Use
#'   `"warn"` to keep the current warning behavior, `"error"` to fail fast for
#'   invalid internal geometries, or `"ignore"` to skip containment checks.
#' @return Modified BBF-class object.
#' @keywords internal
#' @export
setMethod(
  "reforge",
  signature(object = "BBF"),
  function(
    object,
    body_scale = NULL,
    body_target = NULL,
    backbone_scale = NULL,
    backbone_target = NULL,
    isometric_body = TRUE,
    isometric_backbone = TRUE,
    maintain_ratio = TRUE,
    n_segments_body = NULL,
    n_segments_backbone = NULL,
    containment = c("warn", "error", "ignore")
  ) {
    ############################################################################
    # Validation ===============================================================
    containment <- match.arg(containment)
    if (is.null(body_scale) && is.null(backbone_scale) &&
      is.null(body_target) && is.null(backbone_target) &&
      is.null(n_segments_body) && is.null(n_segments_backbone)) {
      stop(
        "Must specify at least one scaling, target, or segment count parameter."
      )
    }
    if (!is.null(body_scale) && !is.null(body_target)) {
      stop("Specify only one of body_scale or body_target, not both.")
    }
    if (!is.null(backbone_scale) && !is.null(backbone_target)) {
      stop("Specify only one of backbone_scale or backbone_target, not both.")
    }
    if ((!is.null(body_scale) || !is.null(body_target)) &&
      (!is.null(backbone_scale) || !is.null(backbone_target)) &&
      maintain_ratio) {
      maintain_ratio <- FALSE
      message(
        "Hidden State Change Warning:\n ",
        "Multiple axes specified for the body and backbone: ",
        "'maintain_ratio' will be ignored for those axes."
      )
    }
    ############################################################################
    # Extract components =======================================================
    body <- acousticTS::extract(object, "body")
    backbone <- acousticTS::extract(object, "backbone")
    shape <- acousticTS::extract(object, "shape_parameters")
    rpos_b <- body$rpos
    rpos_bb <- backbone$rpos
    ############################################################################
    # Calculate current dimensions =============================================
    body_dims <- .reforge_component_dimensions(
      rpos_b,
      length = shape$body$length
    )
    backbone_dims <- .reforge_component_dimensions(
      rpos_bb,
      length = shape$backbone$length %||%
        .shape_length(position_matrix = rpos_bb, row_major = TRUE)
    )
    ############################################################################
    # Preserve the backbone's axial start relative to the body =================
    backbone_relative_start <- .reforge_relative_start(rpos_bb, rpos_b)
    backbone_relative_vertical <- .reforge_relative_vertical_offset(
      rpos_bb,
      rpos_b
    )
    ############################################################################
    # Process target parameters ================================================
    body_target <- .validate_dimensions_target(
      body_target,
      "body_target",
      c("length", "width", "height")
    )
    body_scale_lst <- .reforge_scale_vector(body_scale, body_target, body_dims)
    backbone_target <- .validate_dimensions_target(
      backbone_target,
      "backbone_target",
      c("length", "width", "height")
    )
    backbone_scale_lst <- .reforge_scale_vector(
      backbone_scale,
      backbone_target,
      backbone_dims
    )
    ############################################################################
    # Process scaling parameters ===============================================
    body_scales <- .validate_dimension_scaling(
      dims = body_scale_lst$scale,
      dims_name = paste0("body", body_scale_lst$suffix),
      valid_dims = c("length", "width", "height"),
      isometry = isometric_body,
      iso_name = "isometric_body"
    )
    backbone_scales <- .validate_dimension_scaling(
      dims = backbone_scale_lst$scale,
      dims_name = paste0("backbone", backbone_scale_lst$suffix),
      valid_dims = c("length", "width", "height"),
      isometry = isometric_backbone,
      iso_name = "isometric_backbone"
    )
    ############################################################################
    # Apply ratio maintenance logic ============================================
    if (maintain_ratio) {
      if (!is.null(body_scales) && is.null(backbone_scales)) {
        backbone_scales <- body_scales
      } else if (is.null(body_scales) && !is.null(backbone_scales)) {
        body_scales <- backbone_scales
      }
    }
    ############################################################################
    # Keep the backbone cross-section circular =================================
    if (!is.null(backbone_scales)) {
      width_scale <- unname(backbone_scales["width"])
      height_scale <- unname(backbone_scales["height"])
      width_changed <- !isTRUE(all.equal(width_scale, 1))
      height_changed <- !isTRUE(all.equal(height_scale, 1))

      if (width_changed && height_changed &&
        !isTRUE(all.equal(width_scale, height_scale))) {
        stop(
          "Backbone reforge must preserve a circular cross-section: ",
          "width and height scaling must match.",
          call. = FALSE
        )
      }

      if (width_changed && !height_changed) {
        backbone_scales["height"] <- width_scale
      } else if (!width_changed && height_changed) {
        backbone_scales["width"] <- height_scale
      }
    }
    ############################################################################
    # Interpolate segments first (before scaling) ==============================
    rpos_b <- .reforge_resample_rows(rpos_b, n_segments_body)
    rpos_bb <- .reforge_resample_rows(rpos_bb, n_segments_backbone)
    ############################################################################
    # Apply scaling ============================================================
    rpos_b <- .reforge_apply_axis_scaling(rpos_b, body_scales)
    rpos_bb <- .reforge_apply_axis_scaling(rpos_bb, backbone_scales)
    ############################################################################
    # Restore the backbone's relative body start after body-length changes =====
    if (!is.null(body_scales) && body_scales["length"] != 1) {
      rpos_bb <- .reforge_shift_to_relative_start(
        rpos_bb,
        backbone_relative_start,
        rpos_b
      )
    }
    rpos_bb <- .reforge_shift_to_relative_vertical_offset(
      rpos_bb,
      backbone_relative_vertical,
      rpos_b
    )
    ############################################################################
    # Validate backbone containment ============================================
    .reforge_check_internal_containment(
      rpos_b,
      rpos_bb,
      component_label = "Backbone",
      action = containment
    )
    ############################################################################
    # Refresh radius profiles ==================================================
    body_radius <- .shape_radius_profile(
      position_matrix = rpos_b,
      row_major = TRUE,
      error_context = "BBF body"
    )
    backbone_radius <- .shape_radius_profile(
      position_matrix = rpos_bb,
      row_major = TRUE,
      error_context = "BBF backbone"
    )
    ############################################################################
    # Update object ============================================================
    methods::slot(object, "body")$rpos <- rpos_b
    methods::slot(object, "body")$radius <- body_radius
    methods::slot(object, "backbone")$rpos <- rpos_bb
    methods::slot(object, "backbone")$radius <- backbone_radius
    methods::slot(object, "components")$backbone <- methods::slot(
      object,
      "backbone"
    )

    methods::slot(object, "shape_parameters")$body$length <- .shape_length(
      position_matrix = rpos_b,
      row_major = TRUE
    )
    methods::slot(object, "shape_parameters")$body$radius <- if (
      all(is.na(body_radius))
    ) {
      NA_real_
    } else {
      max(body_radius, na.rm = TRUE)
    }
    methods::slot(object, "shape_parameters")$body$n_segments <- ncol(rpos_b)

    methods::slot(object, "shape_parameters")$backbone$length <- .shape_length(
      position_matrix = rpos_bb,
      row_major = TRUE
    )
    methods::slot(object, "shape_parameters")$backbone$radius <- max(
      backbone_radius,
      na.rm = TRUE
    )
    methods::slot(object, "shape_parameters")$backbone$n_segments <- ncol(
      rpos_bb
    )

    return(object)
  }
)
################################################################################
#' Reforge GAS-class object
#'
#' Resize a gas-filled scatterer by applying an isometric scale factor or
#' specifying a target maximum radius.  Optionally re-discretize the body
#' representation to a new segment count.  The underlying shape (sphere,
#' prolate spheroid, cylinder, arbitrary, etc.) is preserved; scaling is
#' applied uniformly to all axes of the position matrix.
#'
#' @param object GAS-class object.
#' @param scale Single positive scalar applied isometrically to every axis of
#'   the position matrix.  Mutually exclusive with \code{radius_target}.
#' @param radius_target Target \emph{maximum} body radius (m).  The scale
#'   factor is derived as \code{radius_target / max(current_radius)}.  Mutually
#'   exclusive with \code{scale}.
#' @param n_segments New number of discrete segments.  All position-matrix
#'   columns are re-interpolated along the x-axis.
#' @return Modified GAS-class object.
#' @keywords internal
#' @export
setMethod(
  "reforge",
  signature(object = "GAS"),
  function(object,
           scale = NULL,
           radius_target = NULL,
           n_segments = NULL) {
    ############################################################################
    # Validation ===============================================================
    if (is.null(scale) && is.null(radius_target) && is.null(n_segments)) {
      stop(
        "Must specify at least one of: scale, radius_target, or n_segments.",
        call. = FALSE
      )
    }
    if (!is.null(scale) && !is.null(radius_target)) {
      stop("Specify only one of scale or radius_target, not both.",
        call. = FALSE
      )
    }
    if (!is.null(scale) &&
      (!is.numeric(scale) || length(scale) != 1 || scale <= 0)) {
      stop("'scale' must be a single positive number.", call. = FALSE)
    }
    if (!is.null(radius_target) &&
      (!is.numeric(radius_target) || length(radius_target) != 1 ||
        radius_target <= 0)) {
      stop("'radius_target' must be a single positive number.", call. = FALSE)
    }
    if (!is.null(n_segments) &&
      (!is.numeric(n_segments) || length(n_segments) != 1 || n_segments < 1)) {
      stop("'n_segments' must be a single positive integer.", call. = FALSE)
    }
    ############################################################################
    body <- acousticTS::extract(object, "body")
    shape <- acousticTS::extract(object, "shape_parameters")
    rpos <- body$rpos
    # radius may be a scalar (sphere) or a per-point vector (cylinder, etc.)
    current_radius <- shape$radius
    current_max_r <- max(current_radius, na.rm = TRUE)
    # Derive scale from radius_target if given =================================
    if (!is.null(radius_target)) scale <- radius_target / current_max_r
    ############################################################################
    # Resample segments first ==================================================
    # Interpolate every non-x column so that arbitrary shapes are handled.
    if (!is.null(n_segments)) {
      rpos <- .resample_rpos(rpos, as.integer(n_segments) + 1L)
      methods::slot(object, "shape_parameters")$n_segments <-
        as.integer(n_segments)
    }
    ############################################################################
    # Apply scale ==============================================================
    if (!is.null(scale)) {
      rpos <- rpos * scale
      new_radius <- current_radius * scale
      methods::slot(object, "body")$radius <- new_radius
      methods::slot(object, "shape_parameters")$radius <- new_radius
      methods::slot(object, "shape_parameters")$length <-
        max(rpos[, 1]) - min(rpos[, 1])
    }
    methods::slot(object, "body")$rpos <- rpos
    return(object)
  }
)
################################################################################
#' Reforge CAL-class object
#'
#' Resize a calibration sphere by applying an isometric scale factor or
#' specifying a target diameter.  Optionally re-discretize to a new segment
#' count.  CAL objects are always spheres, so the position matrix follows
#' the \code{sphere()} convention (n_points x 5: x, y, z, zU, zL).
#'
#' @param object CAL-class object.
#' @param scale Single positive scale factor applied isometrically. Mutually
#'   exclusive with \code{diameter_target}.
#' @param diameter_target Target sphere diameter (m). Derives the scale factor
#'   internally. Mutually exclusive with \code{scale}.
#' @param n_segments New number of discrete segments along the major axis.
#' @return Modified CAL-class object.
#' @keywords internal
#' @export
setMethod(
  "reforge",
  signature(object = "CAL"),
  function(object,
           scale = NULL,
           diameter_target = NULL,
           n_segments = NULL) {
    ############################################################################
    # Validation ===============================================================
    if (is.null(scale) && is.null(diameter_target) && is.null(n_segments)) {
      stop(
        "Must specify at least one of: scale, diameter_target, or n_segments.",
        call. = FALSE
      )
    }
    if (!is.null(scale) && !is.null(diameter_target)) {
      stop("Specify only one of scale or diameter_target, not both.",
        call. = FALSE
      )
    }
    if (!is.null(scale) &&
      (!is.numeric(scale) || length(scale) != 1 || scale <= 0)) {
      stop("'scale' must be a single positive number.", call. = FALSE)
    }
    if (!is.null(diameter_target) &&
      (!is.numeric(diameter_target) || length(diameter_target) != 1 ||
        diameter_target <= 0)) {
      stop("'diameter_target' must be a single positive number.", call. = FALSE)
    }
    if (!is.null(n_segments) &&
      (!is.numeric(n_segments) || length(n_segments) != 1 || n_segments < 1)) {
      stop("'n_segments' must be a single positive integer.", call. = FALSE)
    }
    ############################################################################
    body <- acousticTS::extract(object, "body")
    shape <- acousticTS::extract(object, "shape_parameters")
    rpos <- body$rpos
    # CAL is always a sphere; radius is a scalar stored in shape_parameters
    current_radius <- shape$radius_body %||% shape$radius
    # Derive scale from diameter_target if given ===============================
    if (!is.null(diameter_target)) {
      scale <- (diameter_target / 2) /
        current_radius
    }
    ############################################################################
    # Resample segments first ==================================================
    if (!is.null(n_segments)) {
      rpos <- .resample_rpos(rpos, as.integer(n_segments) + 1L)
      methods::slot(object, "shape_parameters")$n_segments <-
        as.integer(n_segments)
    }
    ############################################################################
    # Apply scale ==============================================================
    if (!is.null(scale)) {
      rpos <- rpos * scale
      new_radius <- current_radius * scale
      methods::slot(object, "body")$radius <- new_radius
      methods::slot(object, "body")$diameter <- new_radius * 2
      methods::slot(object, "shape_parameters")$radius <- new_radius
      methods::slot(object, "shape_parameters")$diameter <- new_radius * 2
    }
    methods::slot(object, "body")$rpos <- rpos
    return(object)
  }
)
################################################################################
#' Reforge ESS-class object
#'
#' Resize an elastic-shelled scatterer by applying an isometric scale factor or
#' specifying a target maximum shell radius.  Shell thickness can be updated
#' independently, which rescales the fluid body so that the maximum fluid radius
#' equals \code{new_max_shell_radius - shell_thickness} (matching the convention
#' in \code{\link{ess_generate}}).  The underlying shape (sphere, prolate
#' spheroid, cylinder, etc.) is preserved for both shell and fluid bodies;
#' scaling is applied uniformly to all axes of each position matrix.
#'
#' @param object ESS-class object.
#' @param scale Single positive scalar applied isometrically to the shell (and
#'   fluid body, if present).  Mutually exclusive with \code{radius_target}.
#' @param radius_target Target \emph{maximum} outer shell radius (m).  Scale
#'   factor derived as \code{radius_target / max(current_shell_radius)}.
#'   Mutually exclusive with \code{scale}.
#' @param shell_thickness New shell wall thickness (m).  The fluid body is
#'   rescaled so its maximum radius equals
#'   \code{new_max_shell_radius - shell_thickness}.  Can be combined with
#'   \code{scale}/\code{radius_target} or used alone.
#' @param n_segments New number of discrete segments.  All columns of both the
#'   shell and fluid position matrices are re-interpolated along the x-axis.
#' @return Modified ESS-class object.
#' @keywords internal
#' @export
setMethod(
  "reforge",
  signature(object = "ESS"),
  function(object,
           scale = NULL,
           radius_target = NULL,
           shell_thickness = NULL,
           n_segments = NULL) {
    ############################################################################
    # Validation ===============================================================
    if (is.null(scale) && is.null(radius_target) &&
      is.null(shell_thickness) && is.null(n_segments)) {
      stop(
        paste0(
          "Must specify at least one of: scale, radius_target, ",
          "shell_thickness, or n_segments."
        ),
        call. = FALSE
      )
    }
    if (!is.null(scale) && !is.null(radius_target)) {
      stop("Specify only one of scale or radius_target, not both.",
        call. = FALSE
      )
    }
    if (!is.null(scale) &&
      (!is.numeric(scale) || length(scale) != 1 || scale <= 0)) {
      stop("'scale' must be a single positive number.", call. = FALSE)
    }
    if (!is.null(radius_target) &&
      (!is.numeric(radius_target) || length(radius_target) != 1 ||
        radius_target <= 0)) {
      stop("'radius_target' must be a single positive number.", call. = FALSE)
    }
    if (!is.null(shell_thickness) &&
      (!is.numeric(shell_thickness) || length(shell_thickness) != 1 ||
        shell_thickness <= 0)) {
      stop("'shell_thickness' must be a single positive number.", call. = FALSE)
    }
    if (!is.null(n_segments) &&
      (!is.numeric(n_segments) || length(n_segments) != 1 ||
        n_segments < 1)) {
      stop("'n_segments' must be a single positive integer.", call. = FALSE)
    }
    ############################################################################
    shell <- acousticTS::extract(object, "shell")
    fluid <- acousticTS::extract(object, "fluid")
    shape <- acousticTS::extract(object, "shape_parameters")
    rpos_shell <- shell$rpos
    rpos_fluid <- fluid$rpos # may be NULL
    # radius may be a scalar (sphere) or a per-point vector (cylinder, etc.)
    curr_shell_r <- shape$shell$radius
    curr_shell_max_r <- max(curr_shell_r, na.rm = TRUE)
    curr_fluid_r <- shape$fluid$radius # may be NA or a vector
    curr_fluid_max_r <- if (!is.null(curr_fluid_r) && !all(is.na(curr_fluid_r))) {
      max(curr_fluid_r, na.rm = TRUE)
    } else {
      NA_real_
    }
    # Derive scale from radius_target if given =================================
    if (!is.null(radius_target)) scale <- radius_target / curr_shell_max_r
    ############################################################################
    # Resample segments first ==================================================
    # .resample_rpos() is defined in utilities-geometry.R
    if (!is.null(n_segments)) {
      n_new <- as.integer(n_segments) + 1L
      rpos_shell <- .resample_rpos(rpos_shell, n_new)
      if (!is.null(rpos_fluid)) {
        rpos_fluid <- .resample_rpos(rpos_fluid, n_new)
      }
      methods::slot(object, "shape_parameters")$n_segments <-
        as.integer(n_segments)
    }
    ############################################################################
    # Apply scale to shell =====================================================
    if (!is.null(scale)) {
      rpos_shell <- rpos_shell * scale
      new_shell_r <- curr_shell_r * scale
      new_shell_max_r <- curr_shell_max_r * scale
      methods::slot(object, "shell")$radius <- new_shell_r
      methods::slot(object, "shape_parameters")$shell$radius <- new_shell_r
      methods::slot(object, "shape_parameters")$shell$diameter <-
        new_shell_r * 2
      # Scale fluid if present ------------------------------------------------
      if (!is.null(rpos_fluid) && !is.na(curr_fluid_max_r)) {
        if (!is.null(shell_thickness)) {
          # Fluid scaled so max_fluid_radius = new_shell_max_r - shell_thickness
          # (same convention as ess_generate)
          new_fluid_max_r <- new_shell_max_r - shell_thickness
          if (new_fluid_max_r <= 0) {
            stop("shell_thickness exceeds new shell radius.", call. = FALSE)
          }
          fluid_scale <- new_fluid_max_r / curr_fluid_max_r
          rpos_fluid <- rpos_fluid * fluid_scale
          new_fluid_r <- curr_fluid_r * fluid_scale
          methods::slot(object, "shell")$shell_thickness <- shell_thickness
        } else {
          rpos_fluid <- rpos_fluid * scale
          new_fluid_r <- curr_fluid_r * scale
        }
        methods::slot(object, "fluid")$radius <- new_fluid_r
        methods::slot(object, "shape_parameters")$fluid$radius <- new_fluid_r
        methods::slot(object, "shape_parameters")$fluid$diameter <-
          new_fluid_r * 2
      }
    } else if (!is.null(shell_thickness)) {
      # Thickness-only update: rescale fluid, leave shell unchanged ============
      new_fluid_max_r <- curr_shell_max_r - shell_thickness
      if (new_fluid_max_r <= 0) {
        stop("shell_thickness exceeds shell radius.", call. = FALSE)
      }
      if (!is.null(rpos_fluid) && !is.na(curr_fluid_max_r)) {
        fluid_scale <- new_fluid_max_r / curr_fluid_max_r
        rpos_fluid <- rpos_fluid * fluid_scale
        new_fluid_r <- curr_fluid_r * fluid_scale
        methods::slot(object, "fluid")$radius <- new_fluid_r
        methods::slot(object, "shape_parameters")$fluid$radius <- new_fluid_r
        methods::slot(object, "shape_parameters")$fluid$diameter <-
          new_fluid_r * 2
      }
      methods::slot(object, "shell")$shell_thickness <- shell_thickness
    }
    methods::slot(object, "shell")$rpos <- rpos_shell
    if (!is.null(rpos_fluid)) {
      methods::slot(object, "fluid")$rpos <- rpos_fluid
    }
    return(object)
  }
)
################################################################################
#' Reforge FLS-class object.
#' @param object FLS-class object.
#' @param body_scale Proportional scaling to the body length and radius. When a
#'   single value is supplied, both dimensions are scaled together. Otherwise,
#'   this input must be a named numeric vector using `length` and/or `radius`.
#' @param body_target Target dimensions (m) for the body length and/or radius.
#'   This input must be a named numeric vector.
#' @param isometric_body Logical; maintain isometric scaling for body.
#' @param n_segments_body New number of segments along the body profile.
#' @param length Legacy alias for a new body length resize.
#' @param radius Legacy alias for a new maximum body radius.
#' @param length_radius_ratio_constant Legacy toggle controlling whether a
#'   length-only resize also rescales radius.
#' @param n_segments Legacy alias for `n_segments_body`.
#'
#' @keywords internal
#' @export
setMethod(
  "reforge",
  signature(object = "FLS"),
  function(object,
           body_scale = NULL,
           body_target = NULL,
           isometric_body = TRUE,
           n_segments_body = NULL,
           length = NULL,
           radius = NULL,
           length_radius_ratio_constant = TRUE,
           n_segments = NULL) {
    ############################################################################
    # Validation ===============================================================
    if (!is.null(body_scale) && !is.null(body_target)) {
      stop("Specify only one of body_scale or body_target, not both.",
        call. = FALSE
      )
    }
    if (!is.null(n_segments) && !is.null(n_segments_body)) {
      stop("Specify only one of n_segments or n_segments_body, not both.",
        call. = FALSE
      )
    }
    if ((!is.null(length) || !is.null(radius)) &&
      (!is.null(body_scale) || !is.null(body_target))) {
      stop(
        paste0(
          "Use either the legacy length/radius arguments or the new ",
          "body_scale/body_target arguments, not both."
        ),
        call. = FALSE
      )
    }
    if (is.null(body_scale) && is.null(body_target) &&
      is.null(length) && is.null(radius) &&
      is.null(n_segments) && is.null(n_segments_body)) {
      stop(
        paste0(
          "Must specify at least one of: body_scale, body_target, length, ",
          "radius, n_segments, or n_segments_body."
        ),
        call. = FALSE
      )
    }
    ############################################################################
    # Extract shape and body ===================================================
    shape <- acousticTS::extract(object, "shape_parameters")
    body <- acousticTS::extract(object, "body")
    rpos <- body$rpos
    ############################################################################
    # Resolve legacy aliases onto the shared body scale/target pathway =========
    if (!is.null(n_segments)) {
      n_segments_body <- n_segments
    }
    if (!is.null(length) || !is.null(radius)) {
      if (!is.null(length) && !is.null(radius)) {
        body_target <- c(length = length, radius = radius)
        isometric_body <- FALSE
      } else if (!is.null(length)) {
        body_target <- c(length = length)
        isometric_body <- isTRUE(length_radius_ratio_constant)
      } else if (!is.null(radius)) {
        body_target <- c(radius = radius)
        isometric_body <- FALSE
      }
    }
    ############################################################################
    # Normalize body scaling through the shared validation helpers =============
    use_arc_length <- !is.null(body$arc_length) ||
      (!is.null(body$radius_curvature_ratio) &&
        is.finite(body$radius_curvature_ratio))
    body_dims <- c(
      length = if (use_arc_length) {
        .shape_arc_length(position_matrix = rpos, row_major = TRUE)
      } else {
        shape$length %||%
          .shape_length(position_matrix = rpos, row_major = TRUE)
      },
      radius = shape$radius %||%
        max(.shape_radius_profile(rpos, row_major = TRUE), na.rm = TRUE)
    )
    body_target <- .validate_dimensions_target(
      body_target,
      "body_target",
      c("length", "radius")
    )
    body_scale_lst <- .reforge_scale_vector(body_scale, body_target, body_dims)
    body_scales <- .validate_dimension_scaling(
      dims = body_scale_lst$scale,
      dims_name = paste0("body", body_scale_lst$suffix),
      valid_dims = c("length", "radius"),
      isometry = isometric_body,
      iso_name = "isometric_body"
    )
    ############################################################################
    # Resample to the requested body segment count =============================
    if (!is.null(n_segments_body)) {
      rpos <- .reforge_resample_rows(rpos, as.integer(n_segments_body) + 1L)
      methods::slot(object, "shape_parameters")$n_segments <-
        as.integer(n_segments_body)
    }
    ############################################################################
    # Apply the normalized length/radius scaling using the shared helper =======
    if (!is.null(body_scales)) {
      axis_scales <- c(
        length = unname(body_scales["length"]),
        width = unname(body_scales["radius"]),
        height = unname(body_scales["radius"])
      )
      rpos <- .reforge_apply_axis_scaling(
        rpos,
        axis_scales,
        scale_centerline = use_arc_length
      )
    }
    ############################################################################
    # Flush working copies to slots ============================================
    radii <- .shape_radius_profile(
      position_matrix = rpos,
      row_major = TRUE,
      error_context = "FLS body"
    )
    methods::slot(object, "body")$rpos <- rpos
    methods::slot(object, "body")$radius <- radii
    if (use_arc_length) {
      methods::slot(object, "body")$arc_length <- .shape_arc_length(
        position_matrix = rpos,
        row_major = TRUE
      )
      methods::slot(object, "shape_parameters")$length <- .shape_arc_length(
        position_matrix = rpos,
        row_major = TRUE
      )
    } else {
      methods::slot(object, "shape_parameters")$length <- .shape_length(
        position_matrix = rpos,
        row_major = TRUE
      )
    }
    methods::slot(object, "shape_parameters")$radius <- max(radii, na.rm = TRUE)
    return(object)
    ############################################################################
    # Rescale length ===========================================================
    if (!is.null(length)) {
      new_scale <- length / shape$length
      if (length_radius_ratio_constant) {
        # Isometric: scale all axes (x, y/z centerline path, radius rows).
        rpos <- rpos * new_scale
        # Radius vector: follow length scale unless caller also supplied radius.
        if (is.null(radius)) {
          radii <- radii * new_scale
        } else {
          r_scale <- radius / shape$radius
          radii <- radii * r_scale
          # Correct the already-scaled radius rows in rpos so they match.
          correction <- r_scale / new_scale
          if (nrow(rpos) >= 4) {
            rpos[seq(4L, nrow(rpos)), ] <- rpos[seq(4L, nrow(rpos)), ] *
              correction
          }
        }
        methods::slot(object, "shape_parameters")$radius <- max(radii)
      } else {
        # Length-only: move the x-axis while leaving the radius rows intact.
        rpos[1L, ] <- rpos[1L, ] * new_scale
      }
      methods::slot(object, "shape_parameters")$length <-
        abs(diff(range(rpos[1L, ])))
    }
    ############################################################################
    # Rescale radius (standalone — only when length was not specified) =========
    if (!is.null(radius) && is.null(length)) {
      r_scale <- radius / shape$radius
      radii <- radii * r_scale
      if (nrow(rpos) >= 4) {
        rpos[seq(4L, nrow(rpos)), ] <- rpos[seq(4L, nrow(rpos)), ] * r_scale
      }
      methods::slot(object, "shape_parameters")$radius <- max(radii)
    }
    ############################################################################
    # Flush working copies to slots ============================================
    methods::slot(object, "body")$rpos <- rpos
    methods::slot(object, "body")$radius <- radii
    return(object)
  }
)

#' Extract the local method formals from an S4 \code{reforge()} definition
#' @param method_definition S4 \code{MethodDefinition} object.
#' @keywords internal
#' @noRd
.reforge_method_formals <- function(method_definition) {
  # Prefer direct method formals when the wrapper exposes them =================
  direct_formals <- setdiff(names(formals(method_definition)), c("object", "..."))
  if (length(direct_formals) > 0) {
    return(direct_formals)
  }

  # Otherwise inspect the wrapped .local definition created by S4 dispatch =====
  method_body <- body(method_definition@.Data)
  if (!is.language(method_body) || length(method_body) < 2) {
    return(character(0))
  }

  local_assignment <- method_body[[2]]
  if (!is.call(local_assignment) ||
    !identical(as.character(local_assignment[[1]]), "<-")) {
    return(character(0))
  }

  local_fun <- local_assignment[[3]]
  if (!is.function(local_fun)) {
    return(character(0))
  }

  setdiff(names(formals(local_fun)), "object")
}

#' Get reforge parameters from the current method signatures
#' @param object_class Character string of object class
#' @return Character vector of parameter names
#' @keywords internal
#' @noRd
.discover_reforge_params <- function(object_class) {
  # Resolve the first class in the hierarchy that exposes a reforge method =====
  for (cls in object_class) {
    if (!methods::hasMethod("reforge", cls)) {
      next
    }

    method_definition <- methods::selectMethod("reforge", cls)
    params <- .reforge_method_formals(method_definition)
    if (length(params) > 0) {
      return(params)
    }
  }

  # Return an empty set when no matching reforge method exists ================
  character(0)
}
