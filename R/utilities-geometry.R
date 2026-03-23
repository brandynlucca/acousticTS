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
  if (!is.matrix(position_matrix)) {
    stop("'position_matrix' must be a matrix.", call. = FALSE)
  }

  if (row_major) {
    if (!is.null(rownames(position_matrix))) {
      x_idx <- match(c("x", "x_body", "x_bladder"), rownames(position_matrix),
        nomatch = 0
      )
      x_idx <- x_idx[x_idx > 0]
      if (length(x_idx) > 0) {
        return(as.numeric(position_matrix[x_idx[1], ]))
      }
    }
    return(as.numeric(position_matrix[1, ]))
  }

  if (!is.null(colnames(position_matrix))) {
    x_idx <- match(c("x", "x_body", "x_bladder"), colnames(position_matrix),
      nomatch = 0
    )
    x_idx <- x_idx[x_idx > 0]
    if (length(x_idx) > 0) {
      return(as.numeric(position_matrix[, x_idx[1]]))
    }
  }

  as.numeric(position_matrix[, 1])
}

#' Return the index order corresponding to a canonical x-axis orientation
#' @param x Numeric vector of x positions.
#' @param decreasing Logical; whether to return descending x-order.
#' @keywords internal
#' @noRd
.canonicalize_x_order <- function(x, decreasing = FALSE) {
  if (!is.numeric(x)) {
    stop("'x' must be numeric.", call. = FALSE)
  }
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
  if (!is.null(body)) {
    position_matrix <- body$rpos
    row_major <- TRUE
  }

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
  if (!is.null(body)) {
    position_matrix <- body$rpos
    row_major <- TRUE
    explicit_radius <- body$radius
  } else {
    explicit_radius <- NULL
  }

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

  radius_profile <- resolve_explicit(explicit_radius)
  if (!is.null(radius_profile)) {
    return(radius_profile)
  }

  if (!is.null(shape_parameters)) {
    radius_profile <- resolve_explicit(shape_parameters$radius_shape)
    if (!is.null(radius_profile)) {
      return(radius_profile)
    }

    radius_profile <- resolve_explicit(shape_parameters$radius)
    if (!is.null(radius_profile)) {
      return(radius_profile)
    }
  }

  if (row_major) {
    rn <- rownames(position_matrix)
    if (!is.null(rn) && all(c("zU", "zL") %in% rn)) {
      return(as.numeric(abs(position_matrix["zU", ] - position_matrix["zL", ]) / 2))
    }
    if (!is.null(rn) && all(c("zU_body", "zL_body") %in% rn)) {
      return(as.numeric(abs(position_matrix["zU_body", ] - position_matrix["zL_body", ]) / 2))
    }
    if (!is.null(rn) && all(c("zU_bladder", "zL_bladder") %in% rn)) {
      return(as.numeric(abs(position_matrix["zU_bladder", ] - position_matrix["zL_bladder", ]) / 2))
    }
    if (!is.null(rn) && "w" %in% rn) {
      return(as.numeric(abs(position_matrix["w", ]) / 2))
    }
    if (!is.null(rn) && "w_body" %in% rn) {
      return(as.numeric(abs(position_matrix["w_body", ]) / 2))
    }
    if (!is.null(rn) && "w_bladder" %in% rn) {
      return(as.numeric(abs(position_matrix["w_bladder", ]) / 2))
    }
  } else {
    cn <- colnames(position_matrix)
    if (!is.null(cn) && "a" %in% cn) {
      return(as.numeric(abs(position_matrix[, "a"])))
    }
    if (!is.null(cn) && all(c("zU", "zL") %in% cn)) {
      return(as.numeric(abs(position_matrix[, "zU"] - position_matrix[, "zL"]) / 2))
    }
    if (!is.null(cn) && "w" %in% cn) {
      return(as.numeric(abs(position_matrix[, "w"]) / 2))
    }
  }

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
  if (!is.null(body)) {
    position_matrix <- body$rpos
    row_major <- TRUE
  }

  if (row_major) {
    rn <- rownames(position_matrix)
    if (!is.null(rn) && "w" %in% rn) {
      return(as.numeric(position_matrix["w", ]))
    }
    if (!is.null(rn) && "w_body" %in% rn) {
      return(as.numeric(position_matrix["w_body", ]))
    }
    if (!is.null(rn) && "w_bladder" %in% rn) {
      return(as.numeric(position_matrix["w_bladder", ]))
    }
  } else {
    cn <- colnames(position_matrix)
    if (!is.null(cn) && "w" %in% cn) {
      return(as.numeric(position_matrix[, "w"]))
    }
    if (!is.null(cn) && "y" %in% cn) {
      return(as.numeric(position_matrix[, "y"]))
    }
  }

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
  if (!is.null(body)) {
    position_matrix <- body$rpos
    row_major <- TRUE
  }

  if (row_major) {
    rn <- rownames(position_matrix)
    if (!is.null(rn) && all(c("zU", "zL") %in% rn)) {
      return(as.numeric(position_matrix["zU", ] - position_matrix["zL", ]))
    }
    if (!is.null(rn) && all(c("zU_body", "zL_body") %in% rn)) {
      return(as.numeric(position_matrix["zU_body", ] - position_matrix["zL_body", ]))
    }
    if (!is.null(rn) && all(c("zU_bladder", "zL_bladder") %in% rn)) {
      return(as.numeric(position_matrix["zU_bladder", ] - position_matrix["zL_bladder", ]))
    }
  } else {
    cn <- colnames(position_matrix)
    if (!is.null(cn) && all(c("zU", "zL") %in% cn)) {
      return(as.numeric(position_matrix[, "zU"] - position_matrix[, "zL"]))
    }
  }

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
  x_old <- rpos[, 1]
  x_new <- seq(x_old[1], x_old[nrow(rpos)], length.out = n_new)
  rpos_out <- matrix(0, nrow = n_new, ncol = ncol(rpos))
  if (!is.null(colnames(rpos))) colnames(rpos_out) <- colnames(rpos)
  rpos_out[, 1] <- x_new
  for (j in seq_len(ncol(rpos))[-1]) {
    rpos_out[, j] <- stats::approx(x_old, rpos[, j], xout = x_new)$y
  }
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
  x_new <- seq(rpos[1, 1], rpos[1, ncol(rpos)], length.out = n_new)
  rpos_out <- rbind(
    x_new,
    t(vapply(2:nrow(rpos), function(i) {
      stats::approx(rpos[1, ], rpos[i, ], xout = x_new)$y
    }, FUN.VALUE = numeric(n_new)))
  )
  if (!is.null(rownames(rpos))) rownames(rpos_out) <- rownames(rpos)
  rpos_out
}

################################################################################
#' Support function for bending scatterer body shape and position matrix
#' @param object Dataframe or scatterer-class object
#' @param radius_curvature Radius of curvature that can be parameterized either
#' as a ratio relative to body length or actual measurement
#' @param mode Either "ratio" or "measurement"
#'
#' @keywords shape manipulator
#' @rdname brake
#' @export
brake <- function(object, radius_curvature, mode = "ratio") {
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
#' @noRd
brake_scatterer <- function(object, radius_curvature, mode = "ratio") {
  body <- acousticTS::extract(object, "body")
  body_curved <- brake_df(body, radius_curvature, mode)
  slot(object, "body") <- body_curved
  return(object)
}

################################################################################
#' Support rotation function for KRM (swimbladder)
#' @inheritParams body_rotation
#' @export
bladder_rotation <- function(sum_rpos, rpos, theta, k_length) {
  v <- (sum_rpos[1, ] * cos(theta) + sum_rpos[3, ] * sin(theta)) / 2
  v <- matrix(
    data = rep(v, each = k_length),
    ncol = length(v),
    nrow = k_length
  )
  delta_u <- diff(rpos[1, ]) * sin(theta)
  list(v = v, delta_u = delta_u)
}

################################################################################
#' Support rotating function for KRM (body)
#' @param sum_rpos Summed position matrix
#' @param rpos Position matrix
#' @param theta Orientation angle
#' @param k_length Length of wavenumber vector
#' @export
body_rotation <- function(sum_rpos, rpos, theta, k_length) {
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
  vbU <- (axial_sum * cos(theta) + dorsal_sum * sin(theta)) / 2
  vbL <- (axial_sum * cos(theta) + ventral_sum * sin(theta)) / 2
  delta_u <- diff(rpos[1, ]) * sin(theta)
  list(vbU = vbU, vbL = vbL, delta_u = delta_u)
}

################################################################################
#' Discretize vector into separate intervals of different length
#' @param x1 Desired vector/interval
#' @param x0 Original or initial vector/interval
#' @noRd
segmentize <- function(x1, x0) {
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

  lapply(segments_to_update, function(update) {
    x1[update$indices] <<- update$values
  })

  x1
}

################################################################################
#' Resolve nodewise body radius profile for DWBA-style models
#' @param body Body slot from an FLS-class object.
#' @keywords internal
#' @noRd
.dwba_body_radius <- function(body) {
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
  body <- acousticTS::extract(object, "body")
  shape <- acousticTS::extract(object, "shape_parameters")
  shape_class <- if ("shape" %in% names(shape)) shape$shape else NULL

  is_straight_centerline <- is.matrix(body$rpos) &&
    nrow(body$rpos) >= 3 &&
    all(abs(body$rpos[2, ]) <= sqrt(.Machine$double.eps)) &&
    all(abs(body$rpos[3, ]) <= sqrt(.Machine$double.eps))

  canonical_shape <- !is.null(shape_class) &&
    shape_class %in% c("Sphere", "ProlateSpheroid", "Cylinder")

  if (!canonical_shape || !is_straight_centerline) {
    body$radius <- .shape_radius_profile(
      body = body,
      row_major = TRUE,
      error_context = "DWBA/SDWBA"
    )
    methods::slot(object, "body") <- body
    return(object)
  }

  n_nodes <- shape$n_segments + 1
  decreasing_x <- body$rpos[1, 1] > body$rpos[1, ncol(body$rpos)]
  x_nodes <- body$rpos[1, ]
  y_nodes <- rep(0, n_nodes)
  z_nodes <- rep(0, n_nodes)
  x_center <- mean(range(body$rpos[1, ]))

  radius_nodes <- switch(
    shape_class,
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
  .as_dwba_profile(object)
}

################################################################################
#' Rebuild canonical FLS shapes into the nodewise profile used by KRM
#' @param object FLS-class object.
#' @keywords internal
#' @noRd
.as_krm_profile <- function(object) {
  body <- acousticTS::extract(object, "body")
  shape <- acousticTS::extract(object, "shape_parameters")
  shape_class <- if ("shape" %in% names(shape)) shape$shape else NULL

  is_straight_centerline <- is.matrix(body$rpos) &&
    nrow(body$rpos) >= 3 &&
    all(abs(body$rpos[2, ]) <= sqrt(.Machine$double.eps)) &&
    all(abs(body$rpos[3, ]) <= sqrt(.Machine$double.eps))

  canonical_shape <- !is.null(shape_class) &&
    shape_class %in% c("Sphere", "ProlateSpheroid", "Cylinder")

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

  n_nodes <- shape$n_segments + 1

  prof <- switch(
    shape_class,
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
#' @noRd
sdwba_resample <- function(object, n_segments) {
  if (!inherits(object, "FLS")) {
    stop("Object must be of class FLS")
  }

  object <- .as_dwba_profile(object)

  body <- extract(object, "body")
  orig_rpos <- body$rpos
  orig_x <- orig_rpos[1, ]
  n_orig <- length(orig_x)

  x_new <- seq(body$rpos[1, 1],
    body$rpos[1, dim(body$rpos)[2]],
    length.out = n_segments + 1
  )

  x_nearest <- vapply(orig_x, function(x) which.min(abs(x - x_new)), integer(1))
  x_new[x_nearest[2:(n_orig - 1)]] <- orig_x[2:(n_orig - 1)]

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
