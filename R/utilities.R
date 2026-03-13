################################################################################
################################################################################
# UTILITY FUNCTIONS FOR VARIOUS DATA WRANGLING AND MANIPULATION OPERATIONS
################################################################################
################################################################################
# Internal Helper Functions
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

#' Derive density and sound speed contrasts from absolute properties
#'
#' Prefer explicit contrasts `g`/`h` when present; otherwise compute them
#' relative to the provided medium properties when absolute density and/or
#' sound speed are available.
#'
#' @param body List-like object containing optional `g`, `h`, `density`,
#'   and `sound_speed` entries.
#' @param medium_sound_speed Surrounding-medium sound speed.
#' @param medium_density Surrounding-medium density.
#' @return Named list with elements `g` and `h` (may be NA if insufficient
#'   information was provided).
#' @keywords internal
#' @noRd
.derive_contrasts <- function(body, medium_sound_speed, medium_density) {
  h <- if (!is.null(body$h)) {
    body$h
  } else if (!is.null(body$sound_speed)) {
    body$sound_speed / medium_sound_speed
  } else {
    NA
  }

  g <- if (!is.null(body$g)) {
    body$g
  } else if (!is.null(body$density)) {
    body$density / medium_density
  } else {
    NA
  }

  list(g = g, h = h)
}

#' Extract common initialization components
#'
#' Internal helper to extract shape, body, and medium parameters from scatterer
#' objects during model initialization.
#'
#' @param object Scatterer object
#' @param sound_speed_sw Sound speed in seawater (m/s)
#' @param density_sw Density of seawater (kg/m^3)
#' @return List with shape, body, and medium components
#' @keywords internal
#' @noRd
.init_common <- function(object, sound_speed_sw = 1500, density_sw = 1026) {
  list(
    shape = extract(object, "shape_parameters"),
    body = extract(object, "body"),
    medium = data.frame(
      sound_speed = sound_speed_sw,
      density = density_sw
    )
  )
}

#' Calculate all relevant wavenumbers
#'
#' Internal helper to calculate wavenumbers for multiple sound speeds.
#'
#' @param frequency Frequency vector (Hz)
#' @param sound_speeds Named list of sound speeds (m/s)
#' @return Named list of wavenumber vectors
#' @keywords internal
#' @noRd
.calc_wavenumbers <- function(frequency, sound_speeds) {
  stats::setNames(
    lapply(sound_speeds, function(c) wavenumber(frequency, c)),
    names(sound_speeds)
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

#' Extract material properties with fallback logic
#'
#' Internal helper to extract material properties from ESS components with
#' contrast-based fallback calculations.
#'
#' @param component Material component (shell, fluid, etc.)
#' @param ref_sound_speed Reference sound speed for contrast conversion
#' @param ref_density Reference density for contrast conversion
#' @return List of material properties
#' @keywords internal
#' @noRd
.extract_material_props <- function(component, ref_sound_speed = NULL,
                                    ref_density = NULL) {
  # Helper for property extraction with fallback
  get_prop <- function(direct, contrast, ref_val, default = NA) {
    if (!is.null(component[[direct]])) {
      return(component[[direct]])
    }
    if (!is.null(component[[contrast]]) && !is.null(ref_val)) {
      return(component[[contrast]] * ref_val)
    }
    return(default)
  }

  list(
    sound_speed = get_prop("sound_speed", "h", ref_sound_speed),
    density = get_prop("density", "g", ref_density),
    nu = component$nu %||% NA,
    K = component$K %||% NA,
    E = component$E %||% NA,
    G = component$G %||% NA,
    lambda = component$lambda %||% NA
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
  # Batch parameter case
  if (!is.null(batch_by) && param_name %in% batch_by) {
    idx <- simulation_grid[[paste0(param_name, "_idx")]]
    return(batch_values[[param_name]][idx])
  }

  # Function generator case
  if (is.function(param_value)) {
    return(replicate(grid_size, param_value()))
  }

  # Scalar case
  if (length(param_value) == 1) {
    return(rep(param_value, grid_size))
  }

  # Vector case (must match grid size)
  if (length(param_value) == grid_size) {
    return(param_value)
  }

  # Error case
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
  # Initialize top layer =======================================================
  layer <- object
  # Initialize failure state for error messaging ===============================
  fail_state <- FALSE
  accum_feature <- c()
  # ============================================================================
  # Recursively move through feature(s) -- treats as a directory path when 
  # provided a vector
  for (sub_layer in feature) {
    accum_feature <- c(accum_feature, sub_layer)
    if (is(layer, "Scatterer")) {
      # Scatterer object +++++++++++++++++++++++++++++++++++++++++++++++++++++++
      if (!methods::.hasSlot(layer, sub_layer)) {
        stop(sprintf("Scattering object does not have slot '%s'.", sub_layer))
      }
      
      layer <- methods::slot(layer, sub_layer)
    } else if (is.list(layer)) {
      # List +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      if (! sub_layer %in% names(layer)) {
        fail_state <- TRUE 
        break
      }
      
      layer <- layer[[sub_layer]]
    } else if (is.matrix(layer)) {
      # Matrix +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      # Check across row names +++++++++++++++++++++++++++++++++++++++++++++++++
      is_rows <- sub_layer %in% row.names(layer)
      # Check across column names ++++++++++++++++++++++++++++++++++++++++++++++
      is_cols <- sub_layer %in% colnames(layer)
      
      if (is_rows) {
        layer <- layer[sub_layer,]
      } else if (is_cols) {
        layer <- layer[,sub_layer]
      } else {
        fail_state <- TRUE
        break
      }
    } else if (is.vector(layer)) {
      # Named vector +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      if (!sub_layer %in% names(layer)) {
        fail_state <- TRUE
        break
      }
      
      layer <- layer[sub_layer]
    }
  }
  # Return output (error if failed, otherwise the defined layer) ===============
  if (fail_state) {
    stop(
      sprintf("No feature '%s' ", sub_layer),
      sprintf(
        "found at the path '%s' ", 
        paste(accum_feature, collapse=" -> ")
      ),
      "within the supplied Scatterer object."
    )
  } else {
    layer
  }
}
################################################################################
################################################################################
# Uniformly bend body shape
################################################################################
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
  # Check object type ==========================================================
  class_type <- typeof(input)
  # Bend shape =================================================================
  output <- switch(class_type,
    list = brake_df(object, radius_curvature, mode),
    S4 = brake_scatterer(object, radius_curvature, mode)
  )
  # Return output ==============================================================
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
  # Input validation ===========================================================
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
  # Extract object body shape ==================================================
  rpos <- body_df$rpos
  if (ncol(rpos) < 2) {
    stop("Position matrix 'rpos' must have at least two columns.")
  }
  L <- max(rpos[1, ])
  n_segments <- ncol(rpos)
  # Pull in value for use ======================================================
  radius_curvature_new <- switch(mode,
    ratio = radius_curvature,
    measurement = radius_curvature / L
  )
  # Check angle-per-segment ====================================================
  arc_angle_per_segment <- L / (radius_curvature_new * L * (n_segments - 1))
  # ---- Compare against threshold of 10 (arc) degrees in radians ++++++++++++++
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
  # Estimate the position angle of the bent cylinder ===========================
  # Geometry described in Stanton (1989) for total arc length ++++++++++++++++++
  gamma_max <- 0.5 / radius_curvature_new
  theta <- seq(-gamma_max, gamma_max, length.out = n_segments)
  # Rescale to the original body size ==========================================
  x_new <- (radius_curvature_new * L) * sin(theta) + (L / 2)
  z_new <- (radius_curvature_new * L) * (1 - cos(theta))
  # Determine position matrix direction ========================================
  rpos_direction <- ifelse(which.max(body_df$rpos[1, ]) == 1,
    "REV",
    "FWD"
  )
  x_direction <- switch(rpos_direction,
    REV = rev(x_new),
    FWD = x_new
  )
  z_direction <- switch(rpos_direction,
    REV = -rev(z_new),
    FWD = -z_new
  )
  # Calculate new arc length ===================================================
  arc_lengths <- vapply(seq_len(length(x_direction) - 1), function(i) {
    chord_length <- sqrt(
      (x_direction[i + 1] - x_direction[i])^2 +
        (z_direction[i + 1] - z_direction[i])^2
    )
    cl_theta <- 2 * asin(chord_length / (2 * (radius_curvature_new * L)))
    (radius_curvature_new * L) * cl_theta
  }, FUN.VALUE = numeric(1))
  total_arc_length <- sum(arc_lengths)
  # Check for segment thickness relative to curvature ==========================
  # ---- Get center coordinates of the osculating circle [Q] +++++++++++++++++++
  Q <- c(L / 2, radius_curvature_new * L)
  # ---- Calculate the distances +++++++++++++++++++++++++++++++++++++++++++++++
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
  # Update object ==============================================================
  rpos[c(1, 3), ] <- rbind(x_direction, z_direction)
  body_df_new <- body_df
  body_df_new$rpos <- rpos
  body_df_new$arc_length <- total_arc_length
  body_df_new$radius_curvature_ratio <- radius_curvature_new
  # Return object ==============================================================
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
  # Extract object body shape ==================================================
  body <- extract(object, "body")
  # Pull in value for use ======================================================
  body_curved <- brake_df(body, radius_curvature, mode)
  # Update object ==============================================================
  slot(object, "body") <- body_curved
  # Return object ==============================================================
  return(object)
}
################################################################################
#' Support rotation function for KRM (swimbladder)
#' @inheritParams body_rotation
#' @keywords internal
#' @noRd
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
#' @keywords internal
#' @noRd
body_rotation <- function(sum_rpos, rpos, theta, k_length) {
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
  vbU <- (dorsal_sum * cos(theta) + dorsal_sum * sin(theta)) / 2
  vbL <- (ventral_sum * cos(theta) + ventral_sum * sin(theta)) / 2
  delta_u <- diff(rpos[1, ]) * sin(theta)
  list(vbU = vbU, vbL = vbL, delta_u = delta_u)
}
################################################################################
#' Format data for the modal series solution model into the appropriate matrix
#' @param v Vector input.
#' @param limit Modal series limit.
#' @keywords internal
#' @export
#' @noRd
modal_matrix <- function(v, limit) {
  matrix(
    data = rep(v, each = limit + 1),
    ncol = length(v),
    nrow = limit + 1
  )
}
################################################################################
#' Discretize vector into separate intervals of different length
#' @param x1 Desired vector/interval
#' @param x0 Original or initial vector/interval
#' @noRd
segmentize <- function(x1, x0) {
  # For each consecutive pair of original coordinates
  segments_to_update <- Map(function(start, end) {
    # Find points between this pair
    points_between <- which(x1 < start & x1 > end)

    if (length(points_between) > 0) {
      # Redistribute them evenly between start and end
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

  # Apply all updates
  segments_to_update <- segments_to_update[!vapply(
    segments_to_update,
    is.null,
    logical(1)
  )]

  lapply(segments_to_update, function(update) {
    x1[update$indices] <<- update$values
  })

  # Return
  x1
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
  # Check input
  if (!inherits(object, "FLS")) {
    stop("Object must be of class FLS")
  }

  # Extract shape data
  body <- extract(object, "body")
  orig_rpos <- body$rpos
  orig_x <- orig_rpos[1, ]
  n_orig <- length(orig_x)

  # Create new x-coordinate positions (evenly spaced)
  x_new <- seq(body$rpos[1, 1],
    body$rpos[1, dim(body$rpos)[2]],
    length.out = n_segments + 1
  )

  # Find the closest x-coordinate value
  x_nearest <- vapply(orig_x, function(x) which.min(abs(x - x_new)), integer(1))

  # Update `x_new`
  x_new[x_nearest[2:(n_orig - 1)]] <- orig_x[2:(n_orig - 1)]

  # Segmentize the new coordinates
  x_new_seg <- segmentize(x_new, orig_x)
  # ---- Initialize the new position matrix
  new_rpos <- rbind(x = x_new_seg)

  # Interpolate coordinates using splines (vectorized)
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
  # ---- And transpose
  new_rpos <- rbind(
    new_rpos,
    t(rpos_interp)
  )

  # Initialize radius array
  new_radius <- numeric(n_segments + 1)
  decreasing <- orig_x[1] > orig_x[length(orig_x)]
  # ---- Update stepwise-assigned radius values [inplace operation]
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

  # Add the last component to the position matrix
  new_rpos <- rbind(
    new_rpos,
    rbind(
      zL = new_rpos[3, ] - new_radius,
      zU = new_rpos[3, ] + new_radius
    )
  )

  # Update object with new shape
  object@body$rpos <- new_rpos
  object@body$radius <- new_radius
  object@shape_parameters$n_segments <- n_segments
  object@shape_parameters$length <- abs(diff(range(new_rpos[1, ])))

  object
}