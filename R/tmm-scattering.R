################################################################################
# Transition matrix method (TMM) post-processing helpers
################################################################################

# Validate the stored TMM state before any core post-processing helper tries
# to use the retained modal blocks.
#' @noRd
.tmm_require_stored_blocks <- function(object) {
  # Recover the stored TMM model-parameter bundle from the scatterer ===========
  model_params <- acousticTS::extract(object, "model_parameters")$TMM
  if (is.null(model_params)) {
    stop(
      "No stored TMM model was found on the supplied object. Run ",
      "'target_strength(..., model = \"TMM\")' first.",
      call. = FALSE
    )
  }

  # Confirm that retained modal blocks are actually available ==================
  t_store <- model_params$parameters$t_matrix
  if (is.null(t_store) ||
      !length(t_store) ||
      all(vapply(t_store, is.null, logical(1)))) {
    stop(
      "Stored T-matrix blocks are required for this helper. Re-run ",
      "'target_strength(..., model = \"TMM\", store_t_matrix = TRUE)'.",
      call. = FALSE
    )
  }

  # Return the validated stored TMM state ======================================
  model_params
}

# The cylindrical retained path is now geometry-matched and intentionally
# limited to exact monostatic reuse. The helper is kept as a no-op so the
# public post-processing calls can keep one shared structure without warning on
# the supported cylinder workflow.
#' @noRd
.tmm_warn_exploratory_cylinder_blocks <- function(object, model_params) {
  # Keep the shared post-processing hook while suppressing cylinder warnings ===
  invisible(NULL)
}

# Wrap one angle onto [0, 2*pi).
#' @noRd
.tmm_wrap_angle_2pi <- function(angle) {
  # Wrap the input angle onto the principal 0-to-2pi interval ==================
  wrapped <- angle %% (2 * pi)
  ifelse(wrapped < 0, wrapped + 2 * pi, wrapped)
}

# Check whether the supplied receive direction is the exact monostatic
# opposite of the incident direction in the body-fixed convention.
#' @noRd
.tmm_is_monostatic_direction <- function(theta_body,
                                         phi_body,
                                         theta_scatter,
                                         phi_scatter,
                                         tol = 1e-8) {
  # Compute the exact body-fixed monostatic receive direction ==================
  expected_theta <- pi - theta_body
  expected_phi <- .tmm_wrap_angle_2pi(phi_body + pi)
  delta_phi <- abs(.tmm_wrap_angle_2pi(phi_scatter) - expected_phi)
  delta_phi <- min(delta_phi, 2 * pi - delta_phi)

  # Return whether the supplied receive angles match monostatic geometry =======
  abs(theta_scatter - expected_theta) <= tol && delta_phi <= tol
}

# Normalize a possibly-missing scalar angle onto the stored TMM defaults.
#' @noRd
.tmm_scalar_angle <- function(value, default, name) {
  # Fall back to the stored default when an angle is omitted ===================
  angle <- value %||% default
  if (!is.numeric(angle) || length(angle) != 1 || !is.finite(angle)) {
    stop(
      "'",
      name,
      "' must be a single finite angle in radians.",
      call. = FALSE
    )
  }
  # Return the validated scalar angle ==========================================
  as.numeric(angle)
}

# Normalize a numeric angle vector used for grid-style evaluations.
#' @noRd
.tmm_angle_vector <- function(value,
                              default = NULL,
                              n_default = NULL,
                              lower = NULL,
                              upper = NULL,
                              name) {
  # Build the default angle grid when the caller omitted one ===================
  angles <- value
  if (is.null(angles)) {
    if (is.null(default)) {
      stop("'", name, "' could not be resolved.", call. = FALSE)
    }
    if (!is.null(n_default)) {
      angles <- seq(default[1], default[2], length.out = n_default)
    } else {
      angles <- default
    }
  }

  # Validate the supplied vector and its allowed bounds ========================
  if (!is.numeric(angles) || !length(angles) || any(!is.finite(angles))) {
    stop(
      "'",
      name,
      "' must be a non-empty numeric vector of finite angles in radians.",
      call. = FALSE
    )
  }
  if (!is.null(lower) && any(angles < lower)) {
    stop("'", name, "' must be >= ", lower, " radians.", call. = FALSE)
  }
  if (!is.null(upper) && any(angles > upper)) {
    stop("'", name, "' must be <= ", upper, " radians.", call. = FALSE)
  }

  # Return the validated angle vector ==========================================
  as.numeric(angles)
}

# Build cell-style interval weights from a monotone angle grid so that
# user-supplied densities can be normalized into quadrature weights.
#' @noRd
.tmm_interval_weights <- function(x,
                                  lower = NULL,
                                  upper = NULL,
                                  name = "theta_body") {
  # Validate the monotone quadrature grid ======================================
  if (!is.numeric(x) || !length(x) || any(!is.finite(x))) {
    stop(
      "'",
      name,
      "' must be a non-empty numeric vector of finite values.",
      call. = FALSE
    )
  }
  if (length(x) == 1) {
    return(1)
  }

  # Build the midpoint cell edges and their interval widths ====================
  dx <- diff(x)
  if (any(dx <= 0)) {
    stop("'", name, "' must be strictly increasing.", call. = FALSE)
  }

  edges <- c(
    x[1] - dx[1] / 2,
    (x[-1] + x[-length(x)]) / 2,
    x[length(x)] + dx[length(dx)] / 2
  )
  if (!is.null(lower)) {
    edges[1] <- lower
  }
  if (!is.null(upper)) {
    edges[length(edges)] <- upper
  }

  # Return the cell widths used to normalize orientation densities =============
  diff(edges)
}

# Build a symmetric grid-edge vector for cell-based plotting and solid-angle
# quadrature on stored theta-phi scattering grids.
#' @noRd
.tmm_grid_edges <- function(values, lower, upper) {
  # Handle the degenerate one-cell case directly ===============================
  if (length(values) == 1) {
    return(c(lower, upper))
  }

  mids <- (values[-1] + values[-length(values)]) / 2
  # Return the midpoint-derived plotting or integration edges ==================
  c(lower, mids, upper)
}
# Re-evaluate the exact-family finite-cylinder monostatic response for one
# incident angle without rebuilding a dense spherical retained solve.
#' @noRd
.tmm_cylindrical_monostatic_f_bs <- function(acoustics_row,
                                             body_defaults,
                                             shape_parameters,
                                             boundary,
                                             theta_body) {
  # Recover the canonical finite-cylinder geometry and modal settings ==========
  radius_body <- max(as.numeric(shape_parameters$radius), na.rm = TRUE)
  length_body <- as.numeric(shape_parameters$length)[1]
  k1L <- length_body * acoustics_row$k_sw
  k1a <- acoustics_row$k_sw * sin(theta_body) * radius_body
  m_limit <- as.integer(acoustics_row$n_max)[1]
  nu <- neumann(0:m_limit)

  # Build the exact finite-cylinder modal coefficient vector ===================
  if (boundary %in% c("liquid_filled", "gas_filled")) {
    gh <- body_defaults$g_body * body_defaults$h_body
    k2a <- acoustics_row$k_sw * sin(theta_body) / body_defaults$h_body *
      radius_body
    Bm <- .fcms_bm_fluid(k1a, k2a, gh, nu, m_limit)
  } else if (boundary == "fixed_rigid") {
    Bm <- .fcms_bm_fixed_rigid(k1a, nu, m_limit)
  } else if (boundary == "pressure_release") {
    Bm <- .fcms_bm_pressure_release(k1a, nu, m_limit)
  } else {
    stop("Unsupported boundary for cylindrical TMM branch.", call. = FALSE)
  }

  if (is.matrix(Bm)) {
    Bm <- as.vector(Bm[, 1])
  }

  # Combine the axial sinc factor with the retained modal sum ==================
  axial_factor <- sin(k1L * cos(theta_body)) / (k1L * cos(theta_body))
  if (!is.finite(axial_factor)) {
    axial_factor <- 1
  }

  # Return the exact monostatic cylinder scattering amplitude ==================
  if (boundary %in% c("liquid_filled", "gas_filled")) {
    -length_body / pi * axial_factor * sum(Bm, na.rm = TRUE)
  } else {
    1i * length_body / pi * axial_factor * sum(Bm, na.rm = TRUE)
  }
}

# Evaluate the retained cylindrical branch. The current geometry-matched
# cylinder family only supports exact monostatic reuse; general-angle cylinder
# post-processing still needs a separate validated cylindrical operator.
#' @noRd
.tmm_scattering_cylindrical <- function(model_params,
                                        shape_parameters,
                                        theta_body,
                                        phi_body,
                                        theta_scatter,
                                        phi_scatter,
                                        frequency_idx = NULL) {
  # Enforce the current monostatic-only retained-cylinder scope ================
  if (!.tmm_is_monostatic_direction(
    theta_body = theta_body,
    phi_body = phi_body,
    theta_scatter = theta_scatter,
    phi_scatter = phi_scatter
  )) {
    stop(
      "Stored cylindrical TMM post-processing is currently available only for ",
      "the exact monostatic direction. General-angle cylinder scattering ",
      "still ",
      "needs a separate validated cylindrical operator.",
      call. = FALSE
    )
  }

  # Evaluate the exact finite-cylinder monostatic family at each frequency =====
  acoustics <- model_params$parameters$acoustics
  body_defaults <- model_params$body
  boundary <- model_params$parameters$boundary
  idx <- frequency_idx %||% seq_len(nrow(acoustics))

  vapply(
    idx,
    function(i) {
      .tmm_cylindrical_monostatic_f_bs(
        acoustics_row = acoustics[i, , drop = FALSE],
        body_defaults = body_defaults,
        shape_parameters = shape_parameters,
        boundary = boundary,
        theta_body = theta_body
      )
    },
    complex(1)
  )
}

# Evaluate the stored TMM blocks at one stored frequency over a set of
# incident/receive-angle combinations.
#' @noRd
.tmm_scattering_points <- function(model_params,
                                   frequency_idx,
                                   shape_parameters,
                                   theta_body,
                                   phi_body,
                                   theta_scatter,
                                   phi_scatter) {
  # Recover the stored TMM state for the requested frequency ===================
  parameters <- model_params$parameters
  acoustics <- parameters$acoustics
  n_eval <- length(theta_body)

  if (!all(c(length(phi_body), length(theta_scatter), length(phi_scatter)) ==
           n_eval)) {
    stop(
      "Stored TMM point evaluations require equal-length angle vectors.",
      call. = FALSE
    )
  }

  # Dispatch each angle tuple to the active retained-coordinate backend ========
  vapply(
    seq_len(n_eval),
    function(i) {
      if (parameters$coordinate_system == "spherical") {
        .tmm_scattering_spherical(
          t_store = parameters$t_matrix[frequency_idx],
          acoustics = acoustics[frequency_idx, , drop = FALSE],
          theta_body = theta_body[i],
          phi_body = phi_body[i],
          theta_scatter = theta_scatter[i],
          phi_scatter = phi_scatter[i]
        )[1]
      } else if (parameters$coordinate_system == "spheroidal") {
        .tmm_scattering_spheroidal(
          t_store = parameters$t_matrix[frequency_idx],
          acoustics = acoustics[frequency_idx, , drop = FALSE],
          parameters = parameters,
          theta_body = theta_body[i],
          phi_body = phi_body[i],
          theta_scatter = theta_scatter[i],
          phi_scatter = phi_scatter[i]
        )[1]
      } else if (parameters$coordinate_system == "cylindrical") {
        .tmm_scattering_cylindrical(
          model_params = model_params,
          shape_parameters = shape_parameters,
          theta_body = theta_body[i],
          phi_body = phi_body[i],
          theta_scatter = theta_scatter[i],
          phi_scatter = phi_scatter[i],
          frequency_idx = frequency_idx
        )[1]
      } else {
        stop("Unsupported stored TMM coordinate system.", call. = FALSE)
      }
    },
    complex(1)
  )
}

# Evaluate one stored TMM frequency over a full theta-phi receive grid while
# caching the block terms that do not vary with azimuth.
#' @noRd
.tmm_scattering_grid_matrix <- function(model_params,
                                        frequency_idx,
                                        shape_parameters,
                                        theta_body,
                                        phi_body,
                                        theta_scatter,
                                        phi_scatter) {
  # Recover the stored TMM frequency block and allocate the output grid ========
  parameters <- model_params$parameters
  acoustics <- parameters$acoustics
  t_store <- parameters$t_matrix[[frequency_idx]]
  f_scat <- matrix(
    0 + 0i,
    nrow = length(theta_scatter),
    ncol = length(phi_scatter)
  )

  # Evaluate the grid using the active retained-coordinate backend =============
  if (parameters$coordinate_system == "spherical") {
    mu_inc <- cos(theta_body)
    mu_scat <- cos(theta_scatter)
    k_sw <- acoustics$k_sw[frequency_idx]

    for (block in t_store) {
      # Reuse the theta-dependent block terms across all azimuth angles ========
      n_seq <- as.integer(block$n_seq)
      p_inc <- drop(.tmm_assoc_legendre_table(block$m, max(n_seq), mu_inc))
      a_inc <- .tmm_incident_plane_wave_coefficients(
        block$m,
        n_seq,
        mu_inc,
        p_inc
      )
      coeffs <- as.vector(block[["T"]] %*% a_inc)
      p_scat <- .tmm_assoc_legendre_table(block$m, max(n_seq), mu_scat)
      theta_term <- as.vector(
        p_scat[, seq_along(n_seq), drop = FALSE] %*%
          (((-1i)^(n_seq + 1)) * coeffs / k_sw)
      )
      azimuth_term <- cos(block$m * (phi_body - phi_scatter))
      f_scat <- f_scat + tcrossprod(theta_term, azimuth_term)
    }
  } else if (parameters$coordinate_system == "spheroidal") {
    f_scat <- prolate_spheroid_scattering_grid_from_tmatrix_cpp(
      acoustics = acoustics[frequency_idx, , drop = FALSE],
      t_matrix = t_store,
      theta_body = theta_body,
      phi_body = phi_body,
      theta_scatter = theta_scatter,
      phi_scatter = phi_scatter,
      precision = parameters$precision %||% "double"
    )
  } else if (parameters$coordinate_system == "cylindrical") {
    stop(
      "Stored cylindrical TMM grid evaluations are not available yet. ",
      "The current cylindrical retained operator only supports exact ",
      "monostatic reuse.",
      call. = FALSE
    )
  } else {
    stop("Unsupported stored TMM coordinate system.", call. = FALSE)
  }

  # Return the complex scattering grid =========================================
  f_scat
}

# Evaluate the stored spherical-coordinate blocks at an arbitrary scattering
# geometry without rebuilding the boundary solve.
#' @noRd
.tmm_scattering_spherical <- function(t_store,
                                      acoustics,
                                      theta_body,
                                      phi_body,
                                      theta_scatter,
                                      phi_scatter) {
  # Resolve the body-fixed incident and receive geometry =======================
  mu_inc <- cos(theta_body)
  mu_scat <- cos(theta_scatter)
  delta_phi <- phi_body - phi_scatter
  f_scat <- complex(length.out = length(t_store))

  # Evaluate each stored spherical frequency block in turn =====================
  for (i in seq_along(t_store)) {
    k_sw <- acoustics$k_sw[i]
    f_i <- 0 + 0i

    for (block in t_store[[i]]) {
      # Rebuild the incident and receive angular factors for this block ========
      n_seq <- as.integer(block$n_seq)
      p_inc <- drop(.tmm_assoc_legendre_table(block$m, max(n_seq), mu_inc))
      p_scat <- drop(.tmm_assoc_legendre_table(block$m, max(n_seq), mu_scat))
      a_inc <- .tmm_incident_plane_wave_coefficients(
        block$m,
        n_seq,
        mu_inc,
        p_inc
      )
      coeffs <- as.vector(block[["T"]] %*% a_inc)
      azimuth <- cos(block$m * delta_phi)

      f_i <- f_i + sum(coeffs * (((-1i)^(n_seq + 1)) * p_scat * azimuth / k_sw))
    }

    # Store the far-field amplitude on the package-wide normalization ==========
    f_scat[i] <- f_i
  }

  # Return the frequency-wise scattering amplitudes ============================
  f_scat
}

# Evaluate the stored spheroidal-coordinate blocks at an arbitrary scattering
# geometry using the retained prolate modal operator.
#' @noRd
.tmm_scattering_spheroidal <- function(t_store,
                                       acoustics,
                                       parameters,
                                       theta_body,
                                       phi_body,
                                       theta_scatter,
                                       phi_scatter) {
  # Delegate prolate general-angle evaluation to the compiled retained backend =
  prolate_spheroid_scattering_from_tmatrix_cpp(
    acoustics = acoustics,
    t_matrix = t_store,
    theta_body = theta_body,
    phi_body = phi_body,
    theta_scatter = theta_scatter,
    phi_scatter = phi_scatter,
    precision = parameters$precision %||% "double"
  )
}

#' Evaluate scattering from a stored TMM object
#'
#' @description
#' Evaluates the far-field scattering amplitude from a previously computed
#' `TMM` object using the stored transition-matrix blocks. This allows the same
#' retained
#' modal operator to be reused for arbitrary single-target incident and
#' receive-angle combinations without rebuilding the boundary solve.
#'
#' @param object Scatterer-object previously evaluated with
#'   `target_strength(..., model = "TMM", store_t_matrix = TRUE)`.
#' @param theta_body Incident polar angle (radians). Defaults to the stored TMM
#'   incident angle.
#' @param phi_body Incident azimuth angle (radians). Defaults to the stored TMM
#'   incident angle.
#' @param theta_scatter Receive polar angle (radians). Defaults to the exact
#'   monostatic direction, `pi - theta_body`.
#' @param phi_scatter Receive azimuth angle (radians). Defaults to the exact
#'   monostatic direction, `phi_body + pi`.
#'
#' @return A data frame with the frequency, complex scattering amplitude, the
#'   corresponding differential cross section, and its level in dB.
#'
#' @seealso [target_strength()], \code{\link{tmm_average_orientation}}
#' @export
tmm_scattering <- function(object,
                           theta_body = NULL,
                           phi_body = NULL,
                           theta_scatter = NULL,
                           phi_scatter = NULL) {
  # Recover the stored TMM state and default incident geometry =================
  model_params <- .tmm_require_stored_blocks(object)
  .tmm_warn_exploratory_cylinder_blocks(object, model_params)
  parameters <- model_params$parameters
  acoustics <- parameters$acoustics
  defaults <- model_params$body
  shape_parameters <- acousticTS::extract(object, "shape_parameters")

  # Resolve the requested incident and receive directions ======================
  theta_body <- .tmm_scalar_angle(theta_body, defaults$theta_body, "theta_body")
  phi_body <- .tmm_scalar_angle(phi_body, defaults$phi_body %||% pi, "phi_body")
  theta_scatter <- .tmm_scalar_angle(
    theta_scatter, pi - theta_body, "theta_scatter"
  )
  phi_scatter <- .tmm_scalar_angle(
    phi_scatter, phi_body + pi, "phi_scatter"
  )

  # Dispatch the retained evaluation through the active coordinate backend =====
  f_scat <- switch(
    parameters$coordinate_system,
    spherical = .tmm_scattering_spherical(
      t_store = parameters$t_matrix,
      acoustics = acoustics,
      theta_body = theta_body,
      phi_body = phi_body,
      theta_scatter = theta_scatter,
      phi_scatter = phi_scatter
    ),
    spheroidal = .tmm_scattering_spheroidal(
      t_store = parameters$t_matrix,
      acoustics = acoustics,
      parameters = parameters,
      theta_body = theta_body,
      phi_body = phi_body,
      theta_scatter = theta_scatter,
      phi_scatter = phi_scatter
    ),
    cylindrical = .tmm_scattering_cylindrical(
      model_params = model_params,
      shape_parameters = shape_parameters,
      theta_body = theta_body,
      phi_body = phi_body,
      theta_scatter = theta_scatter,
      phi_scatter = phi_scatter
    ),
    stop("Unsupported stored TMM coordinate system.", call. = FALSE)
  )

  # Return the scattering amplitude and derived cross-section summaries ========
  sigma_scat <- .sigma_bs(f_scat)
  data.frame(
    frequency = acoustics$frequency,
    f_scat = f_scat,
    sigma_scat = sigma_scat,
    sigma_scat_dB = db(sigma_scat)
  )
}

#' @noRd
.tmm_validate_orientation_n_theta <- function(n_theta) {
  if (!is.numeric(n_theta) ||
      length(n_theta) != 1 ||
      !is.finite(n_theta) ||
      n_theta < 1 ||
      n_theta %% 1 != 0) {
    stop("'n_theta' must be a single positive integer.", call. = FALSE)
  }

  as.integer(n_theta)
}

#' @noRd
.tmm_validate_orientation_interval <- function(lower, upper, distribution) {
  if (!is.numeric(lower) ||
      !is.numeric(upper) ||
      length(lower) != 1 ||
      length(upper) != 1 ||
      !is.finite(lower) ||
      !is.finite(upper) ||
      lower < 0 ||
      upper > pi ||
      lower >= upper) {
    stop(
      "'",
      distribution,
      "' requires a finite interval [lower, upper] inside [0, pi].",
      call. = FALSE
    )
  }

  list(lower = lower, upper = upper)
}

#' @noRd
.tmm_validate_orientation_normal <- function(mean_theta, sd_theta) {
  if (!is.numeric(mean_theta) ||
      length(mean_theta) != 1 ||
      !is.finite(mean_theta) ||
      !is.numeric(sd_theta) ||
      length(sd_theta) != 1 ||
      !is.finite(sd_theta) ||
      sd_theta <= 0) {
    stop(
      "'mean_theta' and 'sd_theta' must be finite numeric scalars ",
      "and 'sd_theta' must be > 0.",
      call. = FALSE
    )
  }

  list(mean_theta = mean_theta, sd_theta = sd_theta)
}

#' @noRd
.tmm_orientation_quadrature <- function(theta_body, weights) {
  theta_body <- .tmm_angle_vector(
    theta_body,
    lower = 0,
    upper = pi,
    name = "theta_body"
  )

  if (is.null(weights)) {
    weights <- rep(1 / length(theta_body), length(theta_body))
  } else {
    if (!is.numeric(weights) ||
        length(weights) != length(theta_body) ||
        any(!is.finite(weights)) ||
        any(weights < 0) ||
        sum(weights) <= 0) {
      stop(
        "'weights' must be a non-negative numeric vector with the same ",
        "length as 'theta_body'.",
        call. = FALSE
      )
    }
    weights <- weights / sum(weights)
  }

  list(theta_body = theta_body, weights = weights)
}

#' @noRd
.tmm_orientation_pdf <- function(theta_body, pdf) {
  theta_body <- .tmm_angle_vector(
    theta_body,
    lower = 0,
    upper = pi,
    name = "theta_body"
  )

  density_values <- if (is.function(pdf)) {
    as.numeric(pdf(theta_body))
  } else {
    pdf
  }

  list(
    theta_body = theta_body,
    weights = .tmm_distribution_weights(
      theta_body,
      density_values,
      name = "theta_body"
    )
  )
}

#' @noRd
.tmm_orientation_uniform <- function(lower, upper, n_theta) {
  bounds <- .tmm_validate_orientation_interval(
    lower = lower,
    upper = upper,
    distribution = "uniform"
  )
  theta_body <- seq(bounds$lower, bounds$upper, length.out = n_theta)

  list(
    theta_body = theta_body,
    weights = .tmm_distribution_weights(
      theta_body = theta_body,
      density_values = rep(1, length(theta_body)),
      lower = bounds$lower,
      upper = bounds$upper,
      name = "theta_body"
    )
  )
}

#' @noRd
.tmm_orientation_normal_density <- function(mean_theta, sd_theta, n_theta) {
  normal_args <- .tmm_validate_orientation_normal(mean_theta, sd_theta)
  theta_body <- seq(0, pi, length.out = n_theta)

  list(
    theta_body = theta_body,
    weights = .tmm_distribution_weights(
      theta_body = theta_body,
      density_values = stats::dnorm(
        theta_body,
        mean = normal_args$mean_theta,
        sd = normal_args$sd_theta
      ),
      lower = 0,
      upper = pi,
      name = "theta_body"
    )
  )
}

#' @noRd
.tmm_orientation_truncated_normal <- function(mean_theta,
                                              sd_theta,
                                              lower,
                                              upper,
                                              n_theta) {
  normal_args <- .tmm_validate_orientation_normal(mean_theta, sd_theta)
  bounds <- .tmm_validate_orientation_interval(
    lower = lower,
    upper = upper,
    distribution = "truncated_normal"
  )
  theta_body <- seq(bounds$lower, bounds$upper, length.out = n_theta)

  list(
    theta_body = theta_body,
    weights = .tmm_distribution_weights(
      theta_body = theta_body,
      density_values = stats::dnorm(
        theta_body,
        mean = normal_args$mean_theta,
        sd = normal_args$sd_theta
      ),
      lower = bounds$lower,
      upper = bounds$upper,
      name = "theta_body"
    )
  )
}

#' @noRd
.tmm_resolve_orientation_phi <- function(phi_body, n_angles) {
  phi_body <- rep_len(phi_body, n_angles)
  if (!is.numeric(phi_body) || any(!is.finite(phi_body))) {
    stop(
      "'phi_body' must be a finite numeric scalar or vector.",
      call. = FALSE
    )
  }

  phi_body
}

#' Evaluate a 2D scattering grid from a stored TMM object
#'
#' @description
#' Reuses the stored T-matrix blocks to evaluate the far-field scattering
#' response over a two-dimensional receive-angle grid at one stored frequency.
#' This is useful for bistatic scattering maps, heatmaps, and polar-style
#' visualizations without rebuilding the retained modal solve. In the current
#' package build, this helper is available for the spherical and spheroidal
#' stored branches. Stored cylinders intentionally stop at exact monostatic
#' reuse until a validated retained cylinder angular operator is added.
#'
#' @param object Scatterer-object previously evaluated with
#'   `target_strength(..., model = "TMM", store_t_matrix = TRUE)`.
#' @param frequency Stored frequency (Hz) to evaluate. Required when the object
#'   contains more than one stored frequency.
#' @param theta_body Incident polar angle (radians). Defaults to the stored TMM
#'   incident angle.
#' @param phi_body Incident azimuth angle (radians). Defaults to the stored TMM
#'   incident angle.
#' @param theta_scatter Optional vector of receive polar angles (radians).
#'   Defaults to an evenly spaced grid on `[0, pi]`.
#' @param phi_scatter Optional vector of receive azimuth angles (radians).
#'   Defaults to an evenly spaced grid on `[0, 2*pi]`.
#' @param n_theta Number of default polar-angle grid points when
#'   `theta_scatter` is not supplied.
#' @param n_phi Number of default azimuth grid points when `phi_scatter` is not
#'   supplied.
#'
#' @return A list containing the stored frequency, the incident angles used to
#'   build the grid, the receive-angle vectors, and matrices for the complex
#'   scattering amplitude, differential scattering cross section, and its level
#'   in dB.
#'
#' @seealso \code{\link{tmm_scattering}}, \code{\link{tmm_average_orientation}}
#' @export
tmm_scattering_grid <- function(object,
                                frequency = NULL,
                                theta_body = NULL,
                                phi_body = NULL,
                                theta_scatter = NULL,
                                phi_scatter = NULL,
                                n_theta = 91,
                                n_phi = 181) {
  # Recover the stored TMM state and choose the requested frequency ============
  model_params <- .tmm_require_stored_blocks(object)
  .tmm_warn_exploratory_cylinder_blocks(object, model_params)
  parameters <- model_params$parameters
  defaults <- model_params$body
  acoustics <- parameters$acoustics
  shape_parameters <- acousticTS::extract(object, "shape_parameters")
  idx <- .tmm_plot_frequency_index(frequency, acoustics$frequency)

  if (!is.numeric(n_theta) || length(n_theta) != 1 || !is.finite(n_theta) ||
      n_theta < 2 || n_theta %% 1 != 0) {
    stop("'n_theta' must be a single integer >= 2.", call. = FALSE)
  }
  if (!is.numeric(n_phi) || length(n_phi) != 1 || !is.finite(n_phi) ||
      n_phi < 2 || n_phi %% 1 != 0) {
    stop("'n_phi' must be a single integer >= 2.", call. = FALSE)
  }

  # Resolve the incident direction and receive-angle grids =====================
  theta_body <- .tmm_scalar_angle(theta_body, defaults$theta_body, "theta_body")
  phi_body <- .tmm_scalar_angle(phi_body, defaults$phi_body %||% pi, "phi_body")
  theta_scatter <- .tmm_angle_vector(
    theta_scatter,
    default = c(0, pi),
    n_default = n_theta,
    lower = 0,
    upper = pi,
    name = "theta_scatter"
  )
  phi_scatter <- .tmm_angle_vector(
    phi_scatter,
    default = c(0, 2 * pi),
    n_default = n_phi,
    lower = 0,
    upper = 2 * pi,
    name = "phi_scatter"
  )

  # Evaluate the retained operator over the requested scattering grid ==========
  f_scat <- matrix(
    0 + 0i,
    nrow = length(theta_scatter),
    ncol = length(phi_scatter),
    dimnames = list(NULL, NULL)
  )
  f_scat <- .tmm_scattering_grid_matrix(
    model_params = model_params,
    frequency_idx = idx,
    shape_parameters = shape_parameters,
    theta_body = theta_body,
    phi_body = phi_body,
    theta_scatter = theta_scatter,
    phi_scatter = phi_scatter
  )

  # Return the grid and its derived scattering summaries =======================
  sigma_scat <- .sigma_bs(f_scat)
  sigma_scat_dB <- db(sigma_scat)

  list(
    frequency = acoustics$frequency[idx],
    theta_body = theta_body,
    phi_body = phi_body,
    theta_scatter = theta_scatter,
    phi_scatter = phi_scatter,
    f_scat = f_scat,
    sigma_scat = sigma_scat,
    sigma_scat_dB = sigma_scat_dB
  )
}

# Resolve one stored frequency for helpers that need a single retained TMM
# frequency, using the nearest stored match when needed.
#' @noRd
.tmm_plot_frequency_index <- function(frequency, available) {
  # Require an explicit frequency when more than one is stored =================
  if (is.null(frequency)) {
    if (length(available) == 1) {
      return(1L)
    }
    stop(
      "A stored TMM object with multiple frequencies requires a scalar ",
      "'frequency' input.",
      call. = FALSE
    )
  }

  # Validate the requested frequency and find the nearest stored match =========
  if (!is.numeric(frequency) ||
      length(frequency) != 1 ||
      !is.finite(frequency)) {
    stop("'frequency' must be a single finite value in Hz.", call. = FALSE)
  }

  idx <- which.min(abs(available - frequency))
  tol <- max(1e-8 * max(1, abs(frequency)), 1e-6)
  if (abs(available[idx] - frequency) > tol) {
    warning(
      "Stored TMM results do not contain the requested frequency exactly. ",
      "Using the nearest stored frequency at ", available[idx], " Hz."
    )
  }

  # Return the resolved stored-frequency index =================================
  as.integer(idx)
}
