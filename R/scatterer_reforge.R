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
#' bigger_obj <- reforge(obj, length = 0.03)
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
  width_profile <- .shape_width_profile(position_matrix = rpos, row_major = TRUE)
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

#' Apply named axis scalings to a row-major component position matrix
#' @param rpos Row-major position matrix.
#' @param scales Named vector with length/width/height multipliers.
#' @keywords internal
#' @noRd
.reforge_apply_axis_scaling <- function(rpos, scales) {
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

  # Resolve row indices from the internal profile schema =======================
  row_idx <- function(candidates) {
    idx <- match(candidates, rownames(rpos), nomatch = 0)
    idx <- idx[idx > 0]
    if (length(idx) == 0) {
      integer(0)
    } else {
      idx[1]
    }
  }

  # Rescale the axial coordinate relative to its leading anchor point ==========
  if (scales["length"] != 1) {
    x_idx <- row_idx(.geometry_contract_schema()$profile_row_major$x)
    anchor <- rpos[x_idx, 1]
    rpos[x_idx, ] <- anchor + (rpos[x_idx, ] - anchor) * scales["length"]
  }

  # Rescale the lateral width profile when present =============================
  if (scales["width"] != 1) {
    w_idx <- row_idx(.geometry_contract_schema()$profile_row_major$w)
    if (length(w_idx) > 0) {
      rpos[w_idx, ] <- rpos[w_idx, ] * scales["width"]
    }
  }

  # Rescale the vertical profile while preserving its local centerline =========
  if (scales["height"] != 1) {
    z_idx <- row_idx(c("z", "z_body", "z_bladder"))
    zU_idx <- row_idx(.geometry_contract_schema()$profile_row_major$zU)
    zL_idx <- row_idx(.geometry_contract_schema()$profile_row_major$zL)

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
    relative_start * .shape_length(position_matrix = parent_rpos, row_major = TRUE)
  shift_amount <- new_start - component_x[1]
  x_idx <- match(.geometry_contract_schema()$profile_row_major$x, rownames(component_rpos),
    nomatch = 0
  )
  x_idx <- x_idx[x_idx > 0][1]
  component_rpos[x_idx, ] <- component_rpos[x_idx, ] + shift_amount
  # Return the shifted internal component ======================================
  component_rpos
}

#' Warn if an internal component exceeds its containing body bounds
#' @param rpos_b Row-major body position matrix.
#' @param rpos_i Row-major internal-component position matrix.
#' @keywords internal
#' @noRd
.reforge_check_internal_containment <- function(rpos_b, rpos_i) {
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

  if (!all(contained)) {
    warning("Swimbladder exceeds body bounds at some positions.")
  }
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
    n_segments_swimbladder = NULL
  ) {
    ############################################################################
    # Validation ===============================================================
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
      dims=body_scale_lst$scale,
      dims_name=paste0("body", body_scale_lst$suffix),
      valid_dims=c("length", "width", "height"),
      isometry=isometric_body,
      iso_name="isometric_body"
    )
    bladder_scales <- .validate_dimension_scaling(
      dims=swimbladder_scale_lst$scale,
      dims_name=paste0("swimbladder", swimbladder_scale_lst$suffix),
      valid_dims=c("length", "width", "height"),
      isometry=isometric_swimbladder,
      iso_name="isometric_swimbladder"
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
    }
    ############################################################################
    # Validate swimbladder containment =========================================
    .reforge_check_internal_containment(rpos_b, rpos_sb)
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
      stop("Specify only one of scale or radius_target, not both.", call. = FALSE)
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
    body  <- acousticTS::extract(object, "body")
    shape <- acousticTS::extract(object, "shape_parameters")
    rpos  <- body$rpos
    # radius may be a scalar (sphere) or a per-point vector (cylinder, etc.)
    current_radius <- shape$radius
    current_max_r  <- max(current_radius, na.rm = TRUE)
    # Derive scale from radius_target if given =================================
    if (!is.null(radius_target)) scale <- radius_target / current_max_r
    ############################################################################
    # Resample segments first ==================================================
    # Interpolate every non-x column so that arbitrary shapes are handled.
    if (!is.null(n_segments)) {
      rpos <- .resample_rpos(rpos, as.integer(n_segments) + 1L)
      methods::slot(object, "shape_parameters")$n_segments <- as.integer(n_segments)
    }
    ############################################################################
    # Apply scale ==============================================================
    if (!is.null(scale)) {
      rpos           <- rpos * scale
      new_radius     <- current_radius * scale
      methods::slot(object, "body")$radius             <- new_radius
      methods::slot(object, "shape_parameters")$radius  <- new_radius
      methods::slot(object, "shape_parameters")$length  <-
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
           call. = FALSE)
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
    body  <- acousticTS::extract(object, "body")
    shape <- acousticTS::extract(object, "shape_parameters")
    rpos  <- body$rpos
    # CAL is always a sphere; radius is a scalar stored in shape_parameters
    current_radius <- shape$radius_body %||% shape$radius
    # Derive scale from diameter_target if given ===============================
    if (!is.null(diameter_target)) scale <- (diameter_target / 2) / current_radius
    ############################################################################
    # Resample segments first ==================================================
    if (!is.null(n_segments)) {
      rpos <- .resample_rpos(rpos, as.integer(n_segments) + 1L)
      methods::slot(object, "shape_parameters")$n_segments <- as.integer(n_segments)
    }
    ############################################################################
    # Apply scale ==============================================================
    if (!is.null(scale)) {
      rpos           <- rpos * scale
      new_radius     <- current_radius * scale
      methods::slot(object, "body")$radius               <- new_radius
      methods::slot(object, "body")$diameter             <- new_radius * 2
      methods::slot(object, "shape_parameters")$radius   <- new_radius
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
        paste0("Must specify at least one of: scale, radius_target, ",
               "shell_thickness, or n_segments."),
        call. = FALSE
      )
    }
    if (!is.null(scale) && !is.null(radius_target)) {
      stop("Specify only one of scale or radius_target, not both.", call. = FALSE)
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
        (!is.numeric(n_segments) || length(n_segments) != 1 || n_segments < 1)) {
      stop("'n_segments' must be a single positive integer.", call. = FALSE)
    }
    ############################################################################
    shell      <- acousticTS::extract(object, "shell")
    fluid      <- acousticTS::extract(object, "fluid")
    shape      <- acousticTS::extract(object, "shape_parameters")
    rpos_shell <- shell$rpos
    rpos_fluid <- fluid$rpos           # may be NULL
    # radius may be a scalar (sphere) or a per-point vector (cylinder, etc.)
    curr_shell_r     <- shape$shell$radius
    curr_shell_max_r <- max(curr_shell_r, na.rm = TRUE)
    curr_fluid_r     <- shape$fluid$radius  # may be NA or a vector
    curr_fluid_max_r <- if (!is.null(curr_fluid_r) && !all(is.na(curr_fluid_r)))
      max(curr_fluid_r, na.rm = TRUE) else NA_real_
    # Derive scale from radius_target if given =================================
    if (!is.null(radius_target)) scale <- radius_target / curr_shell_max_r
    ############################################################################
    # Resample segments first ==================================================
    # .resample_rpos() is defined in utilities-geometry.R
    if (!is.null(n_segments)) {
      n_new      <- as.integer(n_segments) + 1L
      rpos_shell <- .resample_rpos(rpos_shell, n_new)
      if (!is.null(rpos_fluid))
        rpos_fluid <- .resample_rpos(rpos_fluid, n_new)
      methods::slot(object, "shape_parameters")$n_segments <- as.integer(n_segments)
    }
    ############################################################################
    # Apply scale to shell =====================================================
    if (!is.null(scale)) {
      rpos_shell      <- rpos_shell * scale
      new_shell_r     <- curr_shell_r * scale
      new_shell_max_r <- curr_shell_max_r * scale
      methods::slot(object, "shell")$radius                    <- new_shell_r
      methods::slot(object, "shape_parameters")$shell$radius   <- new_shell_r
      methods::slot(object, "shape_parameters")$shell$diameter <- new_shell_r * 2
      # Scale fluid if present ------------------------------------------------
      if (!is.null(rpos_fluid) && !is.na(curr_fluid_max_r)) {
        if (!is.null(shell_thickness)) {
          # Fluid scaled so max_fluid_radius = new_shell_max_r - shell_thickness
          # (same convention as ess_generate)
          new_fluid_max_r <- new_shell_max_r - shell_thickness
          if (new_fluid_max_r <= 0)
            stop("shell_thickness exceeds new shell radius.", call. = FALSE)
          fluid_scale <- new_fluid_max_r / curr_fluid_max_r
          rpos_fluid  <- rpos_fluid * fluid_scale
          new_fluid_r <- curr_fluid_r * fluid_scale
          methods::slot(object, "shell")$shell_thickness <- shell_thickness
        } else {
          rpos_fluid  <- rpos_fluid * scale
          new_fluid_r <- curr_fluid_r * scale
        }
        methods::slot(object, "fluid")$radius                    <- new_fluid_r
        methods::slot(object, "shape_parameters")$fluid$radius   <- new_fluid_r
        methods::slot(object, "shape_parameters")$fluid$diameter <- new_fluid_r * 2
      }
    } else if (!is.null(shell_thickness)) {
      # Thickness-only update: rescale fluid, leave shell unchanged ============
      new_fluid_max_r <- curr_shell_max_r - shell_thickness
      if (new_fluid_max_r <= 0)
        stop("shell_thickness exceeds shell radius.", call. = FALSE)
      if (!is.null(rpos_fluid) && !is.na(curr_fluid_max_r)) {
        fluid_scale <- new_fluid_max_r / curr_fluid_max_r
        rpos_fluid  <- rpos_fluid * fluid_scale
        new_fluid_r <- curr_fluid_r * fluid_scale
        methods::slot(object, "fluid")$radius                    <- new_fluid_r
        methods::slot(object, "shape_parameters")$fluid$radius   <- new_fluid_r
        methods::slot(object, "shape_parameters")$fluid$diameter <- new_fluid_r * 2
      }
      methods::slot(object, "shell")$shell_thickness <- shell_thickness
    }
    methods::slot(object, "shell")$rpos <- rpos_shell
    if (!is.null(rpos_fluid))
      methods::slot(object, "fluid")$rpos <- rpos_fluid
    return(object)
  }
)
################################################################################
#' Reforge FLS-class object.
#' @param object FLS-class object.
#' @param length  New body length resize.
#' @param radius New radius size
#' @param n_segments New number of segments
#' @param length_radius_ratio_constant Keep length-to-radius ratio based on new
#' length
#'
#' @keywords internal
#' @export
setMethod(
  "reforge",
  signature(object = "FLS"),
  function(object,
           length = NULL,
           radius = NULL,
           length_radius_ratio_constant = TRUE,
           n_segments = NULL) {
    ############################################################################
    # Validation ===============================================================
    if (is.null(length) && is.null(radius) && is.null(n_segments)) {
      stop(
        "Must specify at least one of: length, radius, or n_segments.",
        call. = FALSE
      )
    }
    ############################################################################
    # Extract shape and body ===================================================
    shape <- acousticTS::extract(object, "shape_parameters")
    body  <- acousticTS::extract(object, "body")
    # Working copies updated progressively so later blocks always see
    # the current (possibly resampled / scaled) state.
    rpos  <- body$rpos
    radii <- body$radius
    ############################################################################
    # Resample to new segment count ============================================
    # .resample_rpos_rows() is defined in utilities-geometry.R
    if (!is.null(n_segments)) {
      n_new  <- as.integer(n_segments) + 1L
      x_new  <- seq(rpos[1, 1], rpos[1, ncol(rpos)], length.out = n_new)
      radii  <- stats::approx(x = rpos[1, ], y = radii, xout = x_new)$y
      rpos   <- .resample_rpos_rows(rpos, n_new)
      methods::slot(object, "shape_parameters")$n_segments <- n_segments
    }
    ############################################################################
    # Rescale length ===========================================================
    if (!is.null(length)) {
      new_scale <- length / shape$length
      if (length_radius_ratio_constant) {
        # Isometric: scale all axes (x, y/z centerline path, radius rows).
        rpos  <- rpos * new_scale
        # Radius vector: follow length scale unless caller also supplied radius.
        if (is.null(radius)) {
          radii <- radii * new_scale
        } else {
          r_scale <- radius / shape$radius
          radii   <- radii * r_scale
          # Correct the already-scaled radius rows in rpos so they match.
          correction <- r_scale / new_scale
          if (nrow(rpos) >= 4)
            rpos[seq(4L, nrow(rpos)), ] <- rpos[seq(4L, nrow(rpos)), ] * correction
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
      radii   <- radii * r_scale
      if (nrow(rpos) >= 4)
        rpos[seq(4L, nrow(rpos)), ] <- rpos[seq(4L, nrow(rpos)), ] * r_scale
      methods::slot(object, "shape_parameters")$radius <- max(radii)
    }
    ############################################################################
    # Flush working copies to slots ============================================
    methods::slot(object, "body")$rpos   <- rpos
    methods::slot(object, "body")$radius <- radii
    return(object)
  }
)

#' Get reforge parameters from known method signatures
#' @param object_class Character string of object class
#' @return Character vector of parameter names
#' @keywords internal
#' @noRd
.discover_reforge_params <- function(object_class) {
  # Return the known public parameter set for each scatterer class =============
  switch(object_class,
    "FLS" = c(
      "length", "radius", "length_radius_ratio_constant",
      "n_segments"
    ),
    "SBF" = c(
      "body_scale", "swimbladder_scale", "body_target",
      "swimbladder_target", "swimbladder_inflation_factor",
      "isometric_body", "isometric_swimbladder", "maintain_ratio",
      "n_segments_body", "n_segments_swimbladder"
    ),
    "GAS" = c("scale", "radius_target", "n_segments"),
    "CAL" = c("scale", "diameter_target", "n_segments"),
    "ESS" = c("scale", "radius_target", "shell_thickness", "n_segments"),
    character(0) # Default for unknown classes
  )
}
