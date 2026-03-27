################################################################################
################################################################################
# GEOMETRY, PROFILE, AND RESAMPLING UTILITIES
################################################################################
################################################################################
#' Extract x-axis values from a position matrix
#' @param position_matrix Numeric position matrix.
#' @param row_major Logical; whether the x-axis is stored in the first row
#'   instead of the first column.
#' @keywords internal
#' @noRd
.shape_x <- function(position_matrix, row_major = FALSE) {
  # Validate the position-matrix input before resolving the x axis =============
  if (!is.matrix(position_matrix)) {
    stop("'position_matrix' must be a matrix.", call. = FALSE)
  }

  # Return the canonical x-axis values from the geometry contract ==============
  .geometry_axis_values(
    position_matrix,
    axis = "x",
    row_major = row_major,
    default = if (row_major) position_matrix[1, ] else position_matrix[, 1],
    context = "Position matrix"
  )
}

#' Return the index order corresponding to a canonical x-axis orientation
#' @param x Numeric vector of x positions.
#' @param decreasing Logical; whether to return descending x-order.
#' @keywords internal
#' @noRd
.canonicalize_x_order <- function(x, decreasing = FALSE) {
  # Validate the x-axis input before ordering it ===============================
  if (!is.numeric(x)) {
    stop("'x' must be numeric.", call. = FALSE)
  }
  # Return the stable x-order permutation ======================================
  order(x, seq_along(x), decreasing = decreasing, na.last = TRUE)
}

#' Reorder a position matrix to a canonical x-axis orientation
#' @param position_matrix Numeric position matrix.
#' @param row_major Logical; whether the x-axis is stored in the first row.
#' @param decreasing Logical; whether the returned matrix should be ordered from
#'   large x to small x.
#' @keywords internal
#' @noRd
.canonicalize_position_matrix <- function(position_matrix,
                                          row_major = FALSE,
                                          decreasing = FALSE) {
  # Resolve the canonical x-order permutation for the position matrix ==========
  ord <- .canonicalize_x_order(
    .shape_x(position_matrix, row_major = row_major),
    decreasing = decreasing
  )

  if (row_major) {
    position_matrix[, ord, drop = FALSE]
  } else {
    position_matrix[ord, , drop = FALSE]
  }
}

#' Resolve geometric length from a position matrix or body list
#' @param position_matrix Numeric position matrix.
#' @param body Body list containing an \code{rpos} matrix.
#' @param row_major Logical; whether the x-axis is stored in the first row.
#' @keywords internal
#' @noRd
.shape_length <- function(position_matrix = NULL,
                          body = NULL,
                          row_major = FALSE) {
  # Fall back to the body profile when a body list is supplied =================
  if (!is.null(body)) {
    position_matrix <- body$rpos
    row_major <- TRUE
  }

  # Return the total x-extent of the profile ===================================
  x_vals <- .shape_x(position_matrix, row_major = row_major)
  abs(diff(range(x_vals, na.rm = TRUE)))
}

#' Resolve nodewise radius information from a shape/body representation
#' @param position_matrix Numeric position matrix.
#' @param shape_parameters Shape-parameter list.
#' @param body Body list containing \code{rpos} and optionally \code{radius}.
#' @param row_major Logical; whether the x-axis is stored in the first row.
#' @param allow_scalar Logical; whether scalar radius inputs should be recycled
#'   to the node count.
#' @param error_context Character string used in the fallback error message.
#' @keywords internal
#' @noRd
.shape_radius_profile <- function(position_matrix = NULL,
                                  shape_parameters = NULL,
                                  body = NULL,
                                  row_major = FALSE,
                                  allow_scalar = TRUE,
                                  error_context = "shape") {
  # Resolve body-style inputs into the shared position-matrix pathway ==========
  if (!is.null(body)) {
    position_matrix <- body$rpos
    row_major <- TRUE
    explicit_radius <- body$radius
  } else {
    explicit_radius <- NULL
  }

  # Define the explicit-radius fallback used by several geometry sources =======
  n_nodes <- if (row_major) ncol(position_matrix) else nrow(position_matrix)

  resolve_explicit <- function(values) {
    if (is.null(values) || all(is.na(values))) {
      return(NULL)
    }
    values <- as.numeric(values)
    if (length(values) == n_nodes) {
      return(values)
    }
    if (allow_scalar && length(values) == 1) {
      return(rep(values, n_nodes))
    }
    NULL
  }

  # Prefer explicit radius values already stored on the body ===================
  radius_profile <- resolve_explicit(explicit_radius)
  if (!is.null(radius_profile)) {
    return(radius_profile)
  }

  if (!is.null(shape_parameters)) {
    # Fall back to radius information stored on the shape parameters ===========
    radius_profile <- resolve_explicit(shape_parameters$radius_shape)
    if (!is.null(radius_profile)) {
      return(radius_profile)
    }

    radius_profile <- resolve_explicit(shape_parameters$radius)
    if (!is.null(radius_profile)) {
      return(radius_profile)
    }
  }

  # Recover a named radius axis from the position matrix when available ========
  radius_axis <- .geometry_axis_values(
    position_matrix,
    axis = "radius",
    row_major = row_major,
    default = NULL,
    context = error_context
  )
  if (!is.null(radius_axis)) {
    return(abs(radius_axis))
  }

  # Derive radius from zU/zL profile envelopes when present ====================
  zU <- .geometry_axis_values(
    position_matrix,
    axis = "zU",
    row_major = row_major,
    default = NULL,
    context = error_context
  )
  zL <- .geometry_axis_values(
    position_matrix,
    axis = "zL",
    row_major = row_major,
    default = NULL,
    context = error_context
  )
  if (!is.null(zU) && !is.null(zL)) {
    return(abs(zU - zL) / 2)
  }

  # Fall back to half the width profile when only w is available ===============
  width <- .geometry_axis_values(
    position_matrix,
    axis = "w",
    row_major = row_major,
    default = NULL,
    context = error_context
  )
  if (!is.null(width)) {
    return(abs(width) / 2)
  }

  # Stop when no radius-like information can be recovered ======================
  stop(
    "Unable to resolve the nodewise radius profile for ",
    error_context,
    ".",
    call. = FALSE
  )
}

#' Resolve nodewise width information from a shape/body representation
#' @inheritParams .shape_radius_profile
#' @keywords internal
#' @noRd
.shape_width_profile <- function(position_matrix = NULL,
                                 shape_parameters = NULL,
                                 body = NULL,
                                 row_major = FALSE,
                                 allow_scalar = TRUE) {
  # Resolve body-style inputs into the shared position-matrix pathway ==========
  if (!is.null(body)) {
    position_matrix <- body$rpos
    row_major <- TRUE
  }

  # Prefer an explicit width axis when one exists ==============================
  width <- .geometry_axis_values(
    position_matrix,
    axis = "w",
    row_major = row_major,
    default = NULL,
    context = "width derivation"
  )
  if (!is.null(width)) {
    return(width)
  }

  # Otherwise derive width from the radius profile =============================
  2 * .shape_radius_profile(
    position_matrix = position_matrix,
    shape_parameters = shape_parameters,
    body = body,
    row_major = row_major,
    allow_scalar = allow_scalar,
    error_context = "width derivation"
  )
}

#' Resolve nodewise height information from a shape/body representation
#' @inheritParams .shape_radius_profile
#' @keywords internal
#' @noRd
.shape_height_profile <- function(position_matrix = NULL,
                                  shape_parameters = NULL,
                                  body = NULL,
                                  row_major = FALSE,
                                  allow_scalar = TRUE) {
  # Resolve body-style inputs into the shared position-matrix pathway ==========
  if (!is.null(body)) {
    position_matrix <- body$rpos
    row_major <- TRUE
  }

  # Prefer explicit upper and lower height envelopes when available ============
  zU <- .geometry_axis_values(
    position_matrix,
    axis = "zU",
    row_major = row_major,
    default = NULL,
    context = "height derivation"
  )
  zL <- .geometry_axis_values(
    position_matrix,
    axis = "zL",
    row_major = row_major,
    default = NULL,
    context = "height derivation"
  )
  if (!is.null(zU) && !is.null(zL)) {
    return(zU - zL)
  }

  # Otherwise derive height from the radius profile ============================
  2 * .shape_radius_profile(
    position_matrix = position_matrix,
    shape_parameters = shape_parameters,
    body = body,
    row_major = row_major,
    allow_scalar = allow_scalar,
    error_context = "height derivation"
  )
}

#' Resample a column-major position matrix to a new number of points
#'
#' Interpolates all non-x columns of a position matrix stored as
#' n_points x ncols (column 1 = x-axis) along the x-axis via linear
#' approximation. Used by the GAS, CAL, and ESS reforge methods.
#'
#' @param rpos Numeric matrix (n_points x ncol), column 1 is the x-axis.
#' @param n_new Target number of rows in the output matrix.
#' @return Resampled numeric matrix with dimensions n_new x ncol.
#' @keywords internal
#' @noRd
.resample_rpos <- function(rpos, n_new) {
  # Build the new x-axis grid spanning the original profile ====================
  x_old <- rpos[, 1]
  x_new <- seq(x_old[1], x_old[nrow(rpos)], length.out = n_new)
  rpos_out <- matrix(0, nrow = n_new, ncol = ncol(rpos))
  if (!is.null(colnames(rpos))) colnames(rpos_out) <- colnames(rpos)
  rpos_out[, 1] <- x_new
  # Interpolate each non-x column onto the new grid ============================
  for (j in seq_len(ncol(rpos))[-1]) {
    rpos_out[, j] <- stats::approx(x_old, rpos[, j], xout = x_new)$y
  }
  # Return the resampled column-major position matrix ==========================
  rpos_out
}

#' Resample a row-major position matrix to a new number of points
#'
#' Interpolates all non-x rows of a position matrix stored as
#' nrows x n_points (row 1 = x-axis) along the x-axis via linear
#' approximation. Used by the FLS reforge method.
#'
#' @param rpos Numeric matrix (nrow x n_points), row 1 is the x-axis.
#' @param n_new Target number of columns in the output matrix.
#' @return Resampled numeric matrix with dimensions nrow x n_new.
#' @keywords internal
#' @noRd
.resample_rpos_rows <- function(rpos, n_new) {
  # Build the new x-axis grid spanning the original profile ====================
  x_new <- seq(rpos[1, 1], rpos[1, ncol(rpos)], length.out = n_new)
  # Interpolate each non-x row onto the new grid ===============================
  rpos_out <- rbind(
    x_new,
    t(vapply(2:nrow(rpos), function(i) {
      stats::approx(rpos[1, ], rpos[i, ], xout = x_new)$y
    }, FUN.VALUE = numeric(n_new)))
  )
  if (!is.null(rownames(rpos))) rownames(rpos_out) <- rownames(rpos)
  # Return the resampled row-major position matrix =============================
  rpos_out
}

################################################################################
#' Support function for bending scatterer body shape and position matrix
#' @param object Dataframe or scatterer-class object
#' @param radius_curvature Radius of curvature that can be parameterized either
#' as a ratio relative to body length or actual measurement
#' @param mode Either "ratio" or "measurement"
#' @return A bent version of \code{object}, returned as the same broad object
#'   type with updated geometry and curvature metadata.
#' @keywords shape manipulator
#' @rdname brake
#' @export
brake <- function(object, radius_curvature, mode = "ratio") {
  # Dispatch the bend helper according to the input object type ================
  class_type <- typeof(object)
  output <- switch(class_type,
    list = brake_df(object, radius_curvature, mode),
    S4 = brake_scatterer(object, radius_curvature, mode)
  )
  return(output)
}

################################################################################
#' Support function for bending scatterer position matrix dataframe
#' @param body_df Dataframe object containing body shape information
#' @param radius_curvature Radius of curvature that can be parameterized either
#' as a ratio relative to body length or actual measurement
#' @param mode Either "ratio" or "measurement"
#'
#' @keywords internal
#' @rdname brake_df
#' @noRd
brake_df <- function(body_df, radius_curvature, mode = "ratio") {
  # Validate the body data frame and requested curvature inputs ================
  if (
    !is.list(body_df) || is.null(body_df$rpos) || !is.matrix(body_df$rpos)
  ) {
    stop("Body shape information must be a list with a matrix element 'rpos'.")
  }
  if (!is.numeric(radius_curvature) || radius_curvature <= 0) {
    stop("Radius of curvature must be a positive-only, real number.")
  }
  if (!mode %in% c("ratio", "measurement")) {
    stop("Radius-of-curvature 'mode' must be either 'ratio' or 'measurement'.")
  }

  # Recover the working geometry and normalize the curvature mode ==============
  rpos <- body_df$rpos
  if (ncol(rpos) < 2) {
    stop("Position matrix 'rpos' must have at least two columns.")
  }
  L <- max(rpos[1, ])
  n_segments <- ncol(rpos)

  radius_curvature_new <- switch(mode,
    ratio = radius_curvature,
    measurement = radius_curvature / L
  )

  # Warn when the arc angle becomes too coarse for the segment count ===========
  arc_angle_per_segment <- L / (radius_curvature_new * L * (n_segments - 1))
  if (arc_angle_per_segment > pi / 8) {
    warning(
      paste0(
        "Arc angle per segment [", round(arc_angle_per_segment, 5), "] ",
        "exceeds pi/8; increase 'n_segments' [", n_segments, "] ",
        "or decrease the effective radius of curvature ('radius_curvature') ",
        "relative to body length [", radius_curvature_new, "]. Otherwise, ",
        "model predictions may be unstable and/or unreliable."
      )
    )
  }

  # Construct the bent centerline geometry along the osculating circle =========
  gamma_max <- 0.5 / radius_curvature_new
  theta <- seq(-gamma_max, gamma_max, length.out = n_segments)
  x_new <- (radius_curvature_new * L) * sin(theta) + (L / 2)
  z_new <- (radius_curvature_new * L) * (1 - cos(theta))

  rpos_direction <- ifelse(which.max(body_df$rpos[1, ]) == 1, "REV", "FWD")
  x_direction <- switch(rpos_direction,
    REV = rev(x_new),
    FWD = x_new
  )
  z_direction <- switch(rpos_direction,
    REV = -rev(z_new),
    FWD = -z_new
  )

  # Recompute the arc length and curvature-center distances ====================
  arc_lengths <- vapply(seq_len(length(x_direction) - 1), function(i) {
    chord_length <- sqrt(
      (x_direction[i + 1] - x_direction[i])^2 +
        (z_direction[i + 1] - z_direction[i])^2
    )
    cl_theta <- 2 * asin(chord_length / (2 * (radius_curvature_new * L)))
    (radius_curvature_new * L) * cl_theta
  }, FUN.VALUE = numeric(1))
  total_arc_length <- sum(arc_lengths)

  Q <- c(L / 2, radius_curvature_new * L)
  Q_d <- sqrt(colSums((rbind(x_direction, z_direction) - Q)^2))
  if (any(body_df$radius > Q_d)) {
    warning(
      "One or more body segments [n=", sum(body_df$radius > Q_d), "] have a ",
      "thickness (radius) greater than their distance to the center of the ",
      "osculating circle used to curve the body shape. This means the ",
      "segment(s) will overlap the center of the arc, causing self-",
      "intersection or unrealistic geometry. Consider reducing the segment ",
      "thickness, increasing the radius of curvature to reduce the shape ",
      "bend, or decreasing the number of segments to avoid overlap."
    )
  }

  # Update the body data frame with the curved geometry ========================
  rpos[c(1, 3), ] <- rbind(x_direction, z_direction)
  body_df_new <- body_df
  body_df_new$rpos <- rpos
  body_df_new$arc_length <- total_arc_length
  body_df_new$radius_curvature_ratio <- radius_curvature_new
  return(body_df_new)
}

################################################################################
#' Support function for bending scatterer body shape scatterer object
#' @param object Scatterer-class object
#' @param radius_curvature Radius of curvature that can be parameterized either
#' as a ratio relative to body length or actual measurement
#' @param mode Either "ratio" or "measurement"
#' @rdname brake_scatterer
#' @importFrom methods slot<-
#' @keywords internal
#' @noRd
brake_scatterer <- function(object, radius_curvature, mode = "ratio") {
  # Apply the body-data-frame brake helper and write it back into the object ===
  body <- acousticTS::extract(object, "body")
  body_curved <- brake_df(body, radius_curvature, mode)
  methods::slot(object, "body") <- body_curved
  # Mirror the curvature metadata onto the stored shape parameters ============
  if ("shape_parameters" %in% methods::slotNames(object)) {
    shape_parameters <- methods::slot(object, "shape_parameters")
    shape_parameters$radius_curvature_ratio <-
      body_curved$radius_curvature_ratio
    methods::slot(object, "shape_parameters") <- shape_parameters
  }
  return(object)
}

################################################################################
#' Support rotation function for KRM (swimbladder)
#' @inheritParams body_rotation
#' @keywords internal
#' @noRd
bladder_rotation <- function(sum_rpos, rpos, theta, k_length) {
  # Rotate the swimbladder midpoint coordinates into the KRM frame =============
  v <- (sum_rpos[1, ] * cos(theta) + sum_rpos[3, ] * sin(theta)) / 2
  v <- matrix(
    data = rep(v, each = k_length),
    ncol = length(v),
    nrow = k_length
  )
  delta_u <- diff(rpos[1, ]) * sin(theta)
  # Return the rotated coordinate and segment-length terms =====================
  list(v = v, delta_u = delta_u)
}

################################################################################
#' Support rotating function for KRM (body)
#' @param sum_rpos Summed position matrix
#' @param rpos Position matrix
#' @param theta Orientation angle
#' @param k_length Length of wavenumber vector
#' @keywords internal
#' @noRd
body_rotation <- function(sum_rpos, rpos, theta, k_length) {
  # Expand the summed coordinates across the requested wavenumber grid =========
  axial_sum <- matrix(
    data = rep(sum_rpos[1, ], each = k_length),
    ncol = length(sum_rpos[1, ]),
    nrow = k_length
  )
  dorsal_sum <- matrix(
    data = rep(sum_rpos[3, ], each = k_length),
    ncol = length(sum_rpos[3, ]),
    nrow = k_length
  )
  ventral_sum <- matrix(
    data = rep(sum_rpos[4, ], each = k_length),
    ncol = length(sum_rpos[4, ]),
    nrow = k_length
  )
  # Rotate the dorsal and ventral body surfaces into the KRM frame =============
  vbU <- (axial_sum * cos(theta) + dorsal_sum * sin(theta)) / 2
  vbL <- (axial_sum * cos(theta) + ventral_sum * sin(theta)) / 2
  delta_u <- diff(rpos[1, ]) * sin(theta)
  # Return the rotated body coordinates and segment-length terms ===============
  list(vbU = vbU, vbL = vbL, delta_u = delta_u)
}

################################################################################
#' Discretize vector into separate intervals of different length
#' @param x1 Desired vector/interval
#' @param x0 Original or initial vector/interval
#' @noRd
segmentize <- function(x1, x0) {
  # Identify each interval that needs to be repartitioned ======================
  segments_to_update <- Map(function(start, end) {
    points_between <- which(x1 < start & x1 > end)

    if (length(points_between) > 0) {
      new_values <- seq(
        from = start,
        to = end,
        length.out = length(points_between) + 2
      )[2:(length(points_between) + 1)]
      list(indices = points_between, values = new_values)
    } else {
      NULL
    }
  }, x0[-length(x0)], x0[-1])

  segments_to_update <- segments_to_update[!vapply(
    segments_to_update,
    is.null,
    logical(1)
  )]

  # Update the target grid interval by interval ================================
  lapply(segments_to_update, function(update) {
    x1[update$indices] <<- update$values
  })

  # Return the repartitioned vector ============================================
  x1
}

################################################################################
#' Resolve nodewise body radius profile for DWBA-style models
#' @param body Body slot from an FLS-class object.
#' @keywords internal
#' @noRd
.dwba_body_radius <- function(body) {
  # Recover the nodewise radius profile used by DWBA-style models ==============
  .shape_radius_profile(
    body = body,
    row_major = TRUE,
    error_context = "DWBA/SDWBA"
  )
}

################################################################################
#' Rebuild canonical FLS shapes into the nodewise profile used by DWBA/SDWBA
#' @param object FLS-class object.
#' @keywords internal
#' @noRd
.as_dwba_profile <- function(object) {
  # Recover the body and shape metadata used to build the DWBA profile =========
  body <- acousticTS::extract(object, "body")
  shape <- acousticTS::extract(object, "shape_parameters")
  shape_class <- if ("shape" %in% names(shape)) shape$shape else NULL

  is_straight_centerline <- is.matrix(body$rpos) &&
    nrow(body$rpos) >= 3 &&
    all(abs(body$rpos[2, ]) <= sqrt(.Machine$double.eps)) &&
    all(abs(body$rpos[3, ]) <= sqrt(.Machine$double.eps))

  canonical_shape <- !is.null(shape_class) &&
    shape_class %in% c("Sphere", "ProlateSpheroid", "Cylinder")

  # Preserve noncanonical or already-curved shapes as explicit profiles ========
  if (!canonical_shape || !is_straight_centerline) {
    body$radius <- .shape_radius_profile(
      body = body,
      row_major = TRUE,
      error_context = "DWBA/SDWBA"
    )
    methods::slot(object, "body") <- body
    return(object)
  }

  # Rebuild the canonical node grid for the recognized shape family ============
  n_nodes <- shape$n_segments + 1
  decreasing_x <- body$rpos[1, 1] > body$rpos[1, ncol(body$rpos)]
  x_nodes <- body$rpos[1, ]
  y_nodes <- rep(0, n_nodes)
  z_nodes <- rep(0, n_nodes)
  x_center <- mean(range(body$rpos[1, ]))

  # Evaluate the canonical radius profile for the requested shape ==============
  radius_nodes <- switch(shape_class,
    Sphere = {
      radius_body <- as.numeric(shape$radius)[1]
      v <- seq(0, pi, length.out = n_nodes)
      x_nodes <- if (decreasing_x) {
        x_center + radius_body * cos(v)
      } else {
        x_center - radius_body * cos(v)
      }
      radius_body * sin(v)
    },
    ProlateSpheroid = {
      semi_major <- if ("semimajor_length" %in% names(shape) &&
        !is.null(shape$semimajor_length)) {
        as.numeric(shape$semimajor_length)[1]
      } else {
        as.numeric(shape$length)[1] / 2
      }
      semi_minor <- if ("semiminor_length" %in% names(shape) &&
        !is.null(shape$semiminor_length)) {
        as.numeric(shape$semiminor_length)[1]
      } else {
        as.numeric(shape$radius)[1]
      }
      v <- seq(0, pi, length.out = n_nodes)
      x_nodes <- if (decreasing_x) {
        x_center + semi_major * cos(v)
      } else {
        x_center - semi_major * cos(v)
      }
      semi_minor * sin(v)
    },
    Cylinder = {
      max_radius <- as.numeric(shape$radius)[1]
      x_nodes <- seq(
        from = body$rpos[1, 1],
        to = body$rpos[1, ncol(body$rpos)],
        length.out = n_nodes
      )
      taper_order <- if ("taper_order" %in% names(shape)) {
        shape$taper_order
      } else {
        NA_real_
      }
      if (is.na(taper_order)) {
        rep(max_radius, n_nodes)
      } else {
        x_n_axis <- seq(-1, 1, length.out = n_nodes)
        tapering <- sqrt(pmax(0, 1 - x_n_axis^taper_order))
        max_radius * tapering
      }
    }
  )

  # Write the rebuilt DWBA profile back into the object ========================
  body$rpos <- rbind(
    x = x_nodes,
    y = y_nodes,
    z = z_nodes,
    zU = z_nodes + radius_nodes,
    zL = z_nodes - radius_nodes
  )
  body$radius <- radius_nodes
  methods::slot(object, "body") <- body
  object
}

#' Backward-compatible DWBA profile helper alias
#' @inheritParams .as_dwba_profile
#' @keywords internal
#' @noRd
.dwba_profile_object <- function(object) {
  # Preserve the historical helper alias used by older code ====================
  .as_dwba_profile(object)
}

################################################################################
#' Rebuild canonical FLS shapes into the nodewise profile used by KRM
#' @param object FLS-class object.
#' @keywords internal
#' @noRd
.as_krm_profile <- function(object) {
  # Recover the body and shape metadata used to build the KRM profile ==========
  body <- acousticTS::extract(object, "body")
  shape <- acousticTS::extract(object, "shape_parameters")
  shape_class <- if ("shape" %in% names(shape)) shape$shape else NULL

  is_straight_centerline <- is.matrix(body$rpos) &&
    nrow(body$rpos) >= 3 &&
    all(abs(body$rpos[2, ]) <= sqrt(.Machine$double.eps)) &&
    all(abs(body$rpos[3, ]) <= sqrt(.Machine$double.eps))

  canonical_shape <- !is.null(shape_class) &&
    shape_class %in% c("Sphere", "ProlateSpheroid", "Cylinder")

  # Preserve noncanonical or already-curved shapes as explicit profiles ========
  if (!canonical_shape || !is_straight_centerline) {
    body$radius <- .shape_radius_profile(
      body = body,
      row_major = TRUE,
      error_context = "KRM"
    )
    ord <- .canonicalize_x_order(.shape_x(body$rpos, row_major = TRUE),
      decreasing = TRUE
    )
    body$rpos <- body$rpos[, ord, drop = FALSE]
    if (length(body$radius) == ncol(body$rpos)) {
      body$radius <- body$radius[ord]
    }
    methods::slot(object, "body") <- body
    return(object)
  }

  # Rebuild the canonical node grid for the recognized shape family ============
  n_nodes <- shape$n_segments + 1

  # Evaluate the canonical KRM radius profile for the requested shape ==========
  prof <- switch(shape_class,
    Sphere = {
      radius_body <- as.numeric(shape$radius)[1]
      v <- seq(-pi / 2, pi / 2, length.out = n_nodes)
      list(
        x = radius_body * (sin(v) + 1),
        radius = radius_body * cos(v)
      )
    },
    ProlateSpheroid = {
      semi_major <- if ("semimajor_length" %in% names(shape) &&
        !is.null(shape$semimajor_length)) {
        as.numeric(shape$semimajor_length)[1]
      } else {
        as.numeric(shape$length)[1] / 2
      }
      semi_minor <- if ("semiminor_length" %in% names(shape) &&
        !is.null(shape$semiminor_length)) {
        as.numeric(shape$semiminor_length)[1]
      } else {
        max(as.numeric(shape$radius), na.rm = TRUE)
      }
      v <- seq(-pi / 2, pi / 2, length.out = n_nodes)
      list(
        x = semi_major * (sin(v) + 1),
        radius = semi_minor * cos(v)
      )
    },
    Cylinder = {
      x_nodes <- seq(
        from = min(body$rpos[1, ]),
        to = max(body$rpos[1, ]),
        length.out = n_nodes
      )
      taper_order <- if ("taper_order" %in% names(shape)) {
        shape$taper_order
      } else {
        NA_real_
      }
      max_radius <- max(as.numeric(shape$radius), na.rm = TRUE)
      radius_nodes <- if (is.na(taper_order)) {
        rep(max_radius, n_nodes)
      } else {
        x_n_axis <- seq(-1, 1, length.out = n_nodes)
        max_radius * sqrt(pmax(1 - x_n_axis^taper_order, 0))
      }
      list(x = x_nodes, radius = radius_nodes)
    }
  )

  # Write the rebuilt KRM profile back into the object =========================
  body$rpos <- rbind(
    x = rev(prof$x),
    w = rev(prof$radius * 2),
    zU = rev(prof$radius),
    zL = rev(-prof$radius)
  )
  body$radius <- rev(prof$radius)
  methods::slot(object, "body") <- body
  object
}

#' Backward-compatible KRM profile helper alias
#' @inheritParams .as_krm_profile
#' @keywords internal
#' @noRd
.krm_profile_object <- function(object) {
  # Preserve the historical helper alias used by older code ====================
  .as_krm_profile(object)
}

################################################################################
#' Resample shape for SDWBA model with piecewise constant radius
#'
#' This function resamples the shape of a fluid-like scatterer (FLS) object for
#' use in stochastic distorted wave Born approximation (SDWBA) calculations.
#' The resampling preserves the overall shape of the scatterer while creating
#' a new representation with the specified number of segments. The radius
#' assignment uses a stepwise algorithm to maintain piecewise constant radius
#' values across segments.
#'
#' @param object FLS-class object to resample
#' @param n_segments Number of segments in the resampled shape
#'
#' @return FLS object with resampled shape
#' @keywords internal
#' @noRd
sdwba_resample <- function(object, n_segments) {
  # Validate the object class and rebuild its explicit DWBA profile ============
  if (!inherits(object, "FLS")) {
    stop("Object must be of class FLS")
  }

  object <- .as_dwba_profile(object)

  # Recover the original profile and build the new x grid ======================
  body <- extract(object, "body")
  orig_rpos <- body$rpos
  orig_x <- orig_rpos[1, ]
  n_orig <- length(orig_x)

  x_new <- seq(body$rpos[1, 1],
    body$rpos[1, dim(body$rpos)[2]],
    length.out = n_segments + 1
  )

  # Align resampled nodes to existing interior breakpoints when possible =======
  x_nearest <- vapply(orig_x, function(x) which.min(abs(x - x_new)), integer(1))
  x_new[x_nearest[2:(n_orig - 1)]] <- orig_x[2:(n_orig - 1)]

  # Repartition the new grid and interpolate the centerline coordinates ========
  x_new_seg <- segmentize(x_new, orig_x)
  new_rpos <- rbind(x = x_new_seg)

  rpos_interp <- apply(
    body$rpos[2:3, ],
    1,
    function(y) {
      stats::spline(
        x = body$rpos[1, ],
        y = y,
        xout = x_new_seg
      )$y
    }
  )
  new_rpos <- rbind(
    new_rpos,
    t(rpos_interp)
  )

  # Assign piecewise-constant radii over the new segment layout ================
  new_radius <- numeric(n_segments + 1)
  decreasing <- orig_x[1] > orig_x[length(orig_x)]
  lapply(1:(n_orig - 1), function(i) {
    if (decreasing) {
      indices <- which(x_new_seg <= orig_x[i]) + 1
      indices <- indices[indices <= length(new_radius)]
    } else {
      indices <- which(x_new_seg >= orig_x[i] &
        x_new_seg < orig_x[i + 1])
    }
    new_radius[indices] <<- body$radius[i + 1]
  })

  # Rebuild the profile envelopes and update the scatterer slots ===============
  new_rpos <- rbind(
    new_rpos,
    rbind(
      zL = new_rpos[3, ] - new_radius,
      zU = new_rpos[3, ] + new_radius
    )
  )

  object@body$rpos <- new_rpos
  object@body$radius <- new_radius
  object@shape_parameters$n_segments <- n_segments
  object@shape_parameters$length <- abs(diff(range(new_rpos[1, ])))

  object
}
