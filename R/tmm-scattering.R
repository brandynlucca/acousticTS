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

# Stored sharp-cylinder workflows now use retained axisymmetric blocks, so
# this hook simply stays silent while the public post-processing helpers reuse
# the standard spherical/spheroidal plotting and summary machinery.
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

# Validate the default 2D scattering-grid dimensions requested by the caller.
#' @noRd
.tmm_validate_scattering_grid_dims <- function(n_theta, n_phi) {
  # Require at least two samples along each angular axis =======================
  if (!is.numeric(n_theta) || length(n_theta) != 1 || !is.finite(n_theta) ||
    n_theta < 2 || n_theta %% 1 != 0) {
    stop("'n_theta' must be a single integer >= 2.", call. = FALSE)
  }
  if (!is.numeric(n_phi) || length(n_phi) != 1 || !is.finite(n_phi) ||
    n_phi < 2 || n_phi %% 1 != 0) {
    stop("'n_phi' must be a single integer >= 2.", call. = FALSE)
  }

  invisible(NULL)
}

# Resolve the incident and receive angles used by `tmm_scattering_grid()`.
#' @noRd
.tmm_scattering_grid_angles <- function(theta_body,
                                        phi_body,
                                        theta_scatter,
                                        phi_scatter,
                                        defaults,
                                        n_theta,
                                        n_phi) {
  # Resolve the incident direction from the stored defaults ====================
  theta_body <- .tmm_scalar_angle(theta_body, defaults$theta_body, "theta_body")
  phi_body <- .tmm_scalar_angle(phi_body, defaults$phi_body %||% pi, "phi_body")
  # Resolve the receive-angle grids ============================================
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

  list(
    theta_body = theta_body,
    phi_body = phi_body,
    theta_scatter = theta_scatter,
    phi_scatter = phi_scatter
  )
}

# Assemble the standard return payload for `tmm_scattering_grid()`.
#' @noRd
.tmm_scattering_grid_output <- function(acoustics,
                                        idx,
                                        theta_body,
                                        phi_body,
                                        theta_scatter,
                                        phi_scatter,
                                        f_scat) {
  # Convert the complex field onto differential scattering summaries ===========
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

# Evaluate the legacy cylindrical branch. This keeps the exact finite-cylinder
# family available for sharp-cylinder monostatic reuse without invoking the
# retained axisymmetric post-processing operator.
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

# Interpolate one stored hybrid-grid field value at the requested receive
# angles. This is used by the current axial elastic-shell prolate TMM rebuild
# branch, which stores a validated external hybrid grid rather than analytic
# retained modal coefficients.
#' @noRd
.tmm_scattering_hybrid_grid_point <- function(store_i,
                                              theta_scatter,
                                              phi_scatter) {
  theta_vals <- as.numeric(store_i$theta_scatter)
  phi_vals <- as.numeric(store_i$phi_scatter)
  field <- store_i$f_scat

  phi_target <- .tmm_wrap_angle_2pi(phi_scatter)
  if (abs(phi_target - 2 * pi) <= 1e-10) {
    phi_target <- 2 * pi
  }
  if (phi_target < min(phi_vals) - 1e-12 || phi_target > max(phi_vals) + 1e-12) {
    phi_target <- min(max(phi_target, min(phi_vals)), max(phi_vals))
  }
  theta_target <- min(max(theta_scatter, min(theta_vals)), max(theta_vals))

  theta_idx_hi <- which(theta_vals >= theta_target)[1]
  if (is.na(theta_idx_hi)) {
    theta_idx_hi <- length(theta_vals)
  }
  theta_idx_lo <- max(1L, theta_idx_hi - 1L)
  theta_idx_hi <- if (theta_vals[theta_idx_lo] == theta_target) theta_idx_lo else theta_idx_hi

  phi_idx_hi <- which(phi_vals >= phi_target)[1]
  if (is.na(phi_idx_hi)) {
    phi_idx_hi <- length(phi_vals)
  }
  phi_idx_lo <- max(1L, phi_idx_hi - 1L)
  phi_idx_hi <- if (phi_vals[phi_idx_lo] == phi_target) phi_idx_lo else phi_idx_hi

  if (theta_idx_lo == theta_idx_hi && phi_idx_lo == phi_idx_hi) {
    return(field[theta_idx_lo, phi_idx_lo])
  }

  theta_lo <- theta_vals[theta_idx_lo]
  theta_hi <- theta_vals[theta_idx_hi]
  phi_lo <- phi_vals[phi_idx_lo]
  phi_hi <- phi_vals[phi_idx_hi]

  theta_w <- if (theta_hi > theta_lo) {
    (theta_target - theta_lo) / (theta_hi - theta_lo)
  } else {
    0
  }
  phi_w <- if (phi_hi > phi_lo) {
    (phi_target - phi_lo) / (phi_hi - phi_lo)
  } else {
    0
  }

  f00 <- field[theta_idx_lo, phi_idx_lo]
  f01 <- field[theta_idx_lo, phi_idx_hi]
  f10 <- field[theta_idx_hi, phi_idx_lo]
  f11 <- field[theta_idx_hi, phi_idx_hi]

  (1 - theta_w) * (1 - phi_w) * f00 +
    (1 - theta_w) * phi_w * f01 +
    theta_w * (1 - phi_w) * f10 +
    theta_w * phi_w * f11
}

# Evaluate one stored hybrid-grid branch at arbitrary receive angles.
#' @noRd
.tmm_scattering_hybrid_grid <- function(t_store,
                                        acoustics,
                                        theta_body,
                                        phi_body,
                                        theta_scatter,
                                        phi_scatter) {
  vapply(
    seq_along(t_store),
    function(i) {
      store_i <- t_store[[i]]
      if (!isTRUE(all.equal(theta_body, store_i$theta_body, tolerance = 1e-8)) ||
        !isTRUE(all.equal(.tmm_wrap_angle_2pi(phi_body),
          .tmm_wrap_angle_2pi(store_i$phi_body),
          tolerance = 1e-8
        ))) {
        stop(
          "Stored ESPSMS hybrid-grid TMM currently supports only the stored axial incident geometry.",
          call. = FALSE
        )
      }
      .tmm_scattering_hybrid_grid_point(
        store_i = store_i,
        theta_scatter = theta_scatter,
        phi_scatter = phi_scatter
      )
    },
    complex(1)
  )
}

# Evaluate one stored hybrid-grid frequency over a requested theta-phi receive
# grid by interpolation from the retained external hybrid reference.
#' @noRd
.tmm_scattering_hybrid_grid_matrix <- function(model_params,
                                               frequency_idx,
                                               theta_body,
                                               phi_body,
                                               theta_scatter,
                                               phi_scatter) {
  store_i <- model_params$parameters$t_matrix[[frequency_idx]]
  if (!isTRUE(all.equal(theta_body, store_i$theta_body, tolerance = 1e-8)) ||
    !isTRUE(all.equal(.tmm_wrap_angle_2pi(phi_body),
      .tmm_wrap_angle_2pi(store_i$phi_body),
      tolerance = 1e-8
    ))) {
    stop(
      "Stored ESPSMS hybrid-grid TMM currently supports only the stored axial incident geometry.",
      call. = FALSE
    )
  }

  if (length(theta_scatter) == length(store_i$theta_scatter) &&
    length(phi_scatter) == length(store_i$phi_scatter) &&
    all(abs(theta_scatter - store_i$theta_scatter) <= 1e-12) &&
    all(abs(phi_scatter - store_i$phi_scatter) <= 1e-12)) {
    return(store_i$f_scat)
  }

  out <- matrix(0 + 0i, nrow = length(theta_scatter), ncol = length(phi_scatter))
  for (i in seq_along(theta_scatter)) {
    for (j in seq_along(phi_scatter)) {
      out[i, j] <- .tmm_scattering_hybrid_grid_point(
        store_i = store_i,
        theta_scatter = theta_scatter[i],
        phi_scatter = phi_scatter[j]
      )
    }
  }

  out
}

# Reject public cylinder TMM helpers that imply general-angle bistatic support.
#' @noRd
.tmm_stop_cylinder_bistatic_public <- function(helper) {
  stop(
    "Cylinder 'TMM' bistatic evaluation is outside the current public scope. ",
    "Cylinder TMM currently supports only exact monostatic backscatter and ",
    "orientation-averaged monostatic products. Helper: '", helper, "'.",
    call. = FALSE
  )
}

# Guard the public stored-cylinder helper scope to exact monostatic reuse only.
#' @noRd
.tmm_validate_cylinder_public_scattering <- function(shape_parameters,
                                                     theta_body,
                                                     phi_body,
                                                     theta_scatter,
                                                     phi_scatter,
                                                     helper = "tmm_scattering()") {
  if (!identical(.tmm_shape_name(shape_parameters), "Cylinder")) {
    return(invisible(NULL))
  }

  if (!.tmm_is_monostatic_direction(
    theta_body = theta_body,
    phi_body = phi_body,
    theta_scatter = theta_scatter,
    phi_scatter = phi_scatter
  )) {
    .tmm_stop_cylinder_bistatic_public(helper)
  }

  invisible(NULL)
}

# Evaluate the retained axisymmetric cylinder branch. Exact monostatic requests
# stay on the finite-cylinder family, while general-angle calls use the
# dedicated cylinder-native profile integral rather than the shared spherical
# retained operator.
#' @noRd
.tmm_cylinder_profile_quadrature <- function(shape_parameters,
                                             cylinder_endcap_fraction = NULL) {
  length_body <- as.numeric(shape_parameters$length)[1]
  half_length <- length_body / 2
  radius_profile <- as.numeric(shape_parameters$radius)
  if (!length(radius_profile)) {
    radius_profile <- rep(0, 2)
  }

  z_nodes <- seq(-half_length, half_length, length.out = length(radius_profile))
  taper_order <- if ("taper_order" %in% names(shape_parameters)) {
    as.numeric(shape_parameters$taper_order)[1]
  } else {
    NA_real_
  }

  if (!is.null(cylinder_endcap_fraction) &&
    cylinder_endcap_fraction > 0 &&
    !is.finite(taper_order)) {
    max_radius <- max(radius_profile, na.rm = TRUE)
    cap_length <- min(half_length * cylinder_endcap_fraction, half_length - 1e-9)
    z_cap_center <- half_length - cap_length
    radius_profile <- vapply(
      z_nodes,
      function(z_i) {
        abs_z <- abs(z_i)
        if (abs_z <= z_cap_center) {
          return(max_radius)
        }
        max_radius * sqrt(pmax(
          1 - ((abs_z - z_cap_center) / cap_length)^2,
          0
        ))
      },
      numeric(1)
    )
  }

  if (length(z_nodes) == 1) {
    return(list(
      z = z_nodes,
      radius = radius_profile,
      weight = length_body
    ))
  }

  dz <- diff(z_nodes)
  weights <- c(dz[1], dz[-length(dz)] + dz[-1], dz[length(dz)]) / 2

  list(
    z = z_nodes,
    radius = radius_profile,
    weight = weights
  )
}

#' @noRd
.tmm_cylinder_modal_coefficients <- function(boundary,
                                             k_sw,
                                             theta_body,
                                             theta_scatter,
                                             radius,
                                             n_max,
                                             g_body = NA_real_,
                                             h_body = NA_real_) {
  sin_eff <- sqrt((sin(theta_body)^2 + sin(theta_scatter)^2) / 2)
  k_eff <- pmax(abs(k_sw * sin_eff), 1e-10)
  m_seq <- 0:as.integer(n_max)
  nu <- neumann(m_seq)
  k1a <- k_eff * radius

  coeffs <- switch(boundary,
    liquid_filled = .fcms_bm_fluid(
      k1a = k1a,
      k2a = k_eff / h_body * radius,
      gh = g_body * h_body,
      nu = nu,
      m_limit = n_max
    ),
    gas_filled = .fcms_bm_fluid(
      k1a = k1a,
      k2a = k_eff / h_body * radius,
      gh = g_body * h_body,
      nu = nu,
      m_limit = n_max
    ),
    fixed_rigid = .fcms_bm_fixed_rigid(k1a = k1a, nu = nu, m_limit = n_max),
    pressure_release = .fcms_bm_pressure_release(
      k1a = k1a,
      nu = nu,
      m_limit = n_max
    ),
    stop("Unsupported boundary for cylinder-native TMM branch.", call. = FALSE)
  )

  if (is.matrix(coeffs)) {
    as.vector(coeffs[, 1])
  } else {
    as.vector(coeffs)
  }
}

#' @noRd
.tmm_cylinder_native_prefactor <- function(boundary) {
  if (boundary %in% c("liquid_filled", "gas_filled")) {
    return(-1 / pi)
  }

  1i / pi
}

#' @noRd
.tmm_scattering_cylinder_native <- function(model_params,
                                            shape_parameters,
                                            theta_body,
                                            phi_body,
                                            theta_scatter,
                                            phi_scatter,
                                            frequency_idx = NULL) {
  parameters <- model_params$parameters
  acoustics <- parameters$acoustics
  idx <- frequency_idx %||% seq_len(nrow(acoustics))
  profile <- .tmm_cylinder_profile_quadrature(
    shape_parameters = shape_parameters,
    cylinder_endcap_fraction = parameters$cylinder_endcap_fraction
  )
  boundary <- parameters$boundary
  dphi <- phi_body - phi_scatter

  if (.tmm_is_monostatic_direction(
    theta_body = theta_body,
    phi_body = phi_body,
    theta_scatter = theta_scatter,
    phi_scatter = phi_scatter
  )) {
    return(vapply(
      idx,
      function(i) {
        .tmm_cylindrical_monostatic_f_bs(
          acoustics_row = acoustics[i, , drop = FALSE],
          body_defaults = model_params$body,
          shape_parameters = shape_parameters,
          boundary = parameters$boundary,
          theta_body = theta_body
        )
      },
      complex(1)
    ))
  }

  vapply(
    idx,
    function(i) {
      acoustics_i <- acoustics[i, , drop = FALSE]
      q_axial <- acoustics_i$k_sw * (cos(theta_body) - cos(theta_scatter))
      axial_phase <- exp(1i * q_axial * profile$z)
      modal_integral <- vapply(
        seq_along(profile$z),
        function(j) {
          .tmm_cylinder_modal_coefficients(
            boundary = boundary,
            k_sw = acoustics_i$k_sw,
            theta_body = theta_body,
            theta_scatter = theta_scatter,
            radius = profile$radius[j],
            n_max = acoustics_i$n_max,
            g_body = model_params$body$g_body,
            h_body = model_params$body$h_body
          ) * profile$weight[j] * axial_phase[j]
        },
        complex(length = acoustics_i$n_max + 1L)
      )
      modal_sum <- rowSums(modal_integral)
      azimuth_term <- cos((0:acoustics_i$n_max) * dphi)
      .tmm_cylinder_native_prefactor(boundary) * sum(modal_sum * azimuth_term)
    },
    complex(1)
  )
}

# Backward-compatibility shim for any already-stored axisymmetric cylinder
# objects created before the dedicated cylinder-native backend landed.
#' @noRd
.tmm_scattering_axisymmetric_cylinder <- function(model_params,
                                                  shape_parameters,
                                                  theta_body,
                                                  phi_body,
                                                  theta_scatter,
                                                  phi_scatter,
                                                  frequency_idx = NULL) {
  .tmm_scattering_cylinder_native(
    model_params = model_params,
    shape_parameters = shape_parameters,
    theta_body = theta_body,
    phi_body = phi_body,
    theta_scatter = theta_scatter,
    phi_scatter = phi_scatter,
    frequency_idx = frequency_idx
  )
}

#' @noRd
.tmm_scattering_cylinder_native_grid <- function(model_params,
                                                 frequency_idx,
                                                 shape_parameters,
                                                 theta_body,
                                                 phi_body,
                                                 theta_scatter,
                                                 phi_scatter) {
  parameters <- model_params$parameters
  acoustics_i <- parameters$acoustics[frequency_idx, , drop = FALSE]
  profile <- .tmm_cylinder_profile_quadrature(
    shape_parameters = shape_parameters,
    cylinder_endcap_fraction = parameters$cylinder_endcap_fraction
  )
  boundary <- parameters$boundary
  prefactor <- .tmm_cylinder_native_prefactor(boundary)
  f_scat <- matrix(
    0 + 0i,
    nrow = length(theta_scatter),
    ncol = length(phi_scatter)
  )

  for (i in seq_along(theta_scatter)) {
    theta_i <- theta_scatter[i]
    q_axial <- acoustics_i$k_sw * (cos(theta_body) - cos(theta_i))
    axial_phase <- exp(1i * q_axial * profile$z)
    modal_integral <- vapply(
      seq_along(profile$z),
      function(j) {
        .tmm_cylinder_modal_coefficients(
          boundary = boundary,
          k_sw = acoustics_i$k_sw,
          theta_body = theta_body,
          theta_scatter = theta_i,
          radius = profile$radius[j],
          n_max = acoustics_i$n_max,
          g_body = model_params$body$g_body,
          h_body = model_params$body$h_body
        ) * profile$weight[j] * axial_phase[j]
      },
      complex(length = acoustics_i$n_max + 1L)
    )
    modal_sum <- rowSums(modal_integral)
    azimuth_term <- vapply(
      phi_scatter,
      function(phi_i) cos((0:acoustics_i$n_max) * (phi_body - phi_i)),
      numeric(acoustics_i$n_max + 1L)
    )
    f_scat[i, ] <- prefactor * drop(modal_sum %*% azimuth_term)
  }

  f_scat
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
  if (parameters$coordinate_system == "spheroidal") {
    if (!all(theta_body == theta_body[[1L]]) || !all(phi_body == phi_body[[1L]])) {
      stop(
        "Stored prolate TMM point evaluations currently require a single ",
        "incident direction per batch.",
        call. = FALSE
      )
    }
    incident_internal <- .tmm_public_to_spheroidal_angles(
      theta = theta_body[[1L]],
      phi = phi_body[[1L]]
    )
    scatter_internal <- .tmm_public_to_spheroidal_angles(
      theta = theta_scatter,
      phi = phi_scatter
    )
    return(
      prolate_spheroid_scattering_points_from_tmatrix_cpp(
        acoustics = acoustics[frequency_idx, , drop = FALSE],
        t_matrix = parameters$t_matrix[[frequency_idx]],
        theta_body = incident_internal$theta[[1L]],
        phi_body = incident_internal$phi[[1L]],
        theta_scatter = scatter_internal$theta,
        phi_scatter = scatter_internal$phi,
        precision = parameters$precision %||% "double"
      )
    )
  }

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
      } else if (parameters$coordinate_system == "sphere_modal") {
        .tmm_scattering_sphere_modal(
          t_store = parameters$t_matrix[frequency_idx],
          acoustics = acoustics[frequency_idx, , drop = FALSE],
          theta_body = theta_body[i],
          phi_body = phi_body[i],
          theta_scatter = theta_scatter[i],
          phi_scatter = phi_scatter[i]
        )[1]
      } else if (parameters$coordinate_system == "espsms_hybrid_grid") {
        .tmm_scattering_hybrid_grid(
          t_store = parameters$t_matrix[frequency_idx],
          acoustics = acoustics[frequency_idx, , drop = FALSE],
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
      } else if (parameters$coordinate_system == "axisymmetric") {
        .tmm_scattering_axisymmetric_cylinder(
          model_params = model_params,
          shape_parameters = shape_parameters,
          theta_body = theta_body[i],
          phi_body = phi_body[i],
          theta_scatter = theta_scatter[i],
          phi_scatter = phi_scatter[i],
          frequency_idx = frequency_idx
        )[1]
      } else if (parameters$coordinate_system == "cylinder_native") {
        .tmm_scattering_cylinder_native(
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
  } else if (parameters$coordinate_system %in% c("axisymmetric", "cylinder_native")) {
    f_scat <- .tmm_scattering_cylinder_native_grid(
      model_params = model_params,
      frequency_idx = frequency_idx,
      shape_parameters = shape_parameters,
      theta_body = theta_body,
      phi_body = phi_body,
      theta_scatter = theta_scatter,
      phi_scatter = phi_scatter
    )
  } else if (parameters$coordinate_system == "spheroidal") {
    incident_internal <- .tmm_public_to_spheroidal_angles(
      theta = theta_body,
      phi = phi_body
    )
    # The public-to-spheroidal transform depends on paired (theta, phi)
    # directions, so evaluate the full receive mesh point-by-point before
    # reshaping it back to the standard theta x phi grid.
    theta_mesh <- rep(theta_scatter, times = length(phi_scatter))
    phi_mesh <- rep(phi_scatter, each = length(theta_scatter))
    scatter_internal <- .tmm_public_to_spheroidal_angles(
      theta = theta_mesh,
      phi = phi_mesh
    )
    f_scat <- prolate_spheroid_scattering_points_from_tmatrix_cpp(
      acoustics = acoustics[frequency_idx, , drop = FALSE],
      t_matrix = t_store,
      theta_body = incident_internal$theta[[1L]],
      phi_body = incident_internal$phi[[1L]],
      theta_scatter = scatter_internal$theta,
      phi_scatter = scatter_internal$phi,
      precision = parameters$precision %||% "double"
    )
    f_scat <- matrix(
      f_scat,
      nrow = length(theta_scatter),
      ncol = length(phi_scatter)
    )
  } else if (parameters$coordinate_system == "sphere_modal") {
    f_scat <- .tmm_scattering_sphere_modal_grid(
      model_params = model_params,
      frequency_idx = frequency_idx,
      theta_body = theta_body,
      phi_body = phi_body,
      theta_scatter = theta_scatter,
      phi_scatter = phi_scatter
    )
  } else if (parameters$coordinate_system == "espsms_hybrid_grid") {
    f_scat <- .tmm_scattering_hybrid_grid_matrix(
      model_params = model_params,
      frequency_idx = frequency_idx,
      theta_body = theta_body,
      phi_body = phi_body,
      theta_scatter = theta_scatter,
      phi_scatter = phi_scatter
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

# Resolve cos(gamma) for the scattering angle between the incident and receive
# directions in the body-fixed spherical convention.
#' @noRd
.tmm_scattering_cosine <- function(theta_body,
                                   phi_body,
                                   theta_scatter,
                                   phi_scatter) {
  cos(theta_body) * cos(theta_scatter) +
    sin(theta_body) * sin(theta_scatter) * cos(phi_body - phi_scatter)
}

# Evaluate the exact stored sphere-modal coefficients at arbitrary receive
# directions. This is used by the shell-sphere TMM branch.
#' @noRd
.tmm_scattering_sphere_modal <- function(t_store,
                                         acoustics,
                                         theta_body,
                                         phi_body,
                                         theta_scatter,
                                         phi_scatter) {
  mu <- .tmm_scattering_cosine(
    theta_body = theta_body,
    phi_body = phi_body,
    theta_scatter = theta_scatter,
    phi_scatter = phi_scatter
  )

  vapply(
    seq_along(t_store),
    function(i) {
      store_i <- t_store[[i]]
      p_n <- drop(.tmm_assoc_legendre_table(0L, max(store_i$n_seq), mu))
      -1i / acoustics$k_sw[i] * sum(
        (2 * store_i$n_seq + 1) * store_i$A_n * p_n[seq_along(store_i$n_seq)],
        na.rm = TRUE
      )
    },
    complex(1)
  )
}

#' @noRd
.tmm_scattering_sphere_modal_grid <- function(model_params,
                                              frequency_idx,
                                              theta_body,
                                              phi_body,
                                              theta_scatter,
                                              phi_scatter) {
  store_i <- model_params$parameters$t_matrix[[frequency_idx]]
  acoustics_i <- model_params$parameters$acoustics[frequency_idx, , drop = FALSE]
  mu <- outer(
    theta_scatter,
    phi_scatter,
    Vectorize(function(theta_val, phi_val) {
      .tmm_scattering_cosine(
        theta_body = theta_body,
        phi_body = phi_body,
        theta_scatter = theta_val,
        phi_scatter = phi_val
      )
    })
  )
  p_mat <- array(0, dim = c(length(theta_scatter), length(phi_scatter), length(store_i$n_seq)))
  for (j in seq_along(store_i$n_seq)) {
    p_mat[, , j] <- matrix(
      .tmm_assoc_legendre_table(0L, max(store_i$n_seq), c(mu))[,
        store_i$n_seq[j] + 1],
      nrow = length(theta_scatter),
      ncol = length(phi_scatter)
    )
  }

  f_scat <- matrix(0 + 0i, nrow = length(theta_scatter), ncol = length(phi_scatter))
  for (j in seq_along(store_i$n_seq)) {
    f_scat <- f_scat + (2 * store_i$n_seq[j] + 1) * store_i$A_n[j] * p_mat[, , j]
  }

  -1i / acoustics_i$k_sw[1] * f_scat
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

# Evaluate the stored spherical-coordinate blocks at an arbitrary scattering
# geometry using the already-solved outgoing coefficients for the stored
# incident direction. This avoids one extra T %*% a_inc reconstruction when
# the caller reuses the exact stored geometry.
#' @noRd

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
  incident_internal <- .tmm_public_to_spheroidal_angles(
    theta = theta_body,
    phi = phi_body
  )
  scatter_internal <- .tmm_public_to_spheroidal_angles(
    theta = theta_scatter,
    phi = phi_scatter
  )
  # Delegate prolate general-angle evaluation to the compiled retained backend =
  prolate_spheroid_scattering_from_tmatrix_cpp(
    acoustics = acoustics,
    t_matrix = t_store,
    theta_body = incident_internal$theta,
    phi_body = incident_internal$phi,
    theta_scatter = scatter_internal$theta,
    phi_scatter = scatter_internal$phi,
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
#' receive-angle combinations without rebuilding the boundary solve. In the
#' current package build, stored cylinders are limited to exact monostatic
#' reuse through this helper; public general-angle cylinder bistatic evaluation
#' remains outside scope.
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
  .tmm_validate_cylinder_public_scattering(
    shape_parameters = shape_parameters,
    theta_body = theta_body,
    phi_body = phi_body,
    theta_scatter = theta_scatter,
    phi_scatter = phi_scatter,
    helper = "tmm_scattering()"
  )

  # Dispatch the retained evaluation through the active coordinate backend =====
  f_scat <- switch(parameters$coordinate_system,
    spherical = .tmm_scattering_spherical(
      t_store = parameters$t_matrix,
      acoustics = acoustics,
      theta_body = theta_body,
      phi_body = phi_body,
      theta_scatter = theta_scatter,
      phi_scatter = phi_scatter
    ),
    sphere_modal = .tmm_scattering_sphere_modal(
      t_store = parameters$t_matrix,
      acoustics = acoustics,
      theta_body = theta_body,
      phi_body = phi_body,
      theta_scatter = theta_scatter,
      phi_scatter = phi_scatter
    ),
    espsms_hybrid_grid = .tmm_scattering_hybrid_grid(
      t_store = parameters$t_matrix,
      acoustics = acoustics,
      theta_body = theta_body,
      phi_body = phi_body,
      theta_scatter = theta_scatter,
      phi_scatter = phi_scatter
    ),
    axisymmetric = .tmm_scattering_axisymmetric_cylinder(
      model_params = model_params,
      shape_parameters = shape_parameters,
      theta_body = theta_body,
      phi_body = phi_body,
      theta_scatter = theta_scatter,
      phi_scatter = phi_scatter
    ),
    cylinder_native = .tmm_scattering_cylinder_native(
      model_params = model_params,
      shape_parameters = shape_parameters,
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
#' stored branches. Stored cylinders are intentionally outside the public scope
#' of this helper and stop with an explicit error.
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
  .tmm_validate_scattering_grid_dims(n_theta, n_phi)

  # Resolve the incident direction and receive-angle grids =====================
  angles <- .tmm_scattering_grid_angles(
    theta_body = theta_body,
    phi_body = phi_body,
    theta_scatter = theta_scatter,
    phi_scatter = phi_scatter,
    defaults = defaults,
    n_theta = n_theta,
    n_phi = n_phi
  )
  theta_body <- angles$theta_body
  phi_body <- angles$phi_body
  theta_scatter <- angles$theta_scatter
  phi_scatter <- angles$phi_scatter
  if (identical(.tmm_shape_name(shape_parameters), "Cylinder")) {
    .tmm_stop_cylinder_bistatic_public("tmm_scattering_grid()")
  }

  # Evaluate the retained operator over the requested scattering grid ==========
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
  .tmm_scattering_grid_output(
    acoustics = acoustics,
    idx = idx,
    theta_body = theta_body,
    phi_body = phi_body,
    theta_scatter = theta_scatter,
    phi_scatter = phi_scatter,
    f_scat = f_scat
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
