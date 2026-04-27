################################################################################
# Transition matrix method (TMM) spherical-coordinate branch helpers
################################################################################

# Stable odd double-factorial used in the associated Legendre seed term.
#' @noRd
.tmm_double_factorial_odd <- function(m) {
  # Handle the seed case used by the associated Legendre recursion =============
  if (m <= 0) {
    return(1)
  }

  # Evaluate the odd double factorial through gamma functions ==================
  exp(lgamma(2 * m + 1) - m * log(2) - lgamma(m + 1))
}

# Build the associated Legendre table P_n^m(mu) for one azimuthal block.
#' @noRd
.tmm_assoc_legendre_table <- function(m, n_max, mu) {
  # Initialize the retained degree sequence and output table ===================
  mu <- as.numeric(mu)
  n_seq <- m:n_max
  n_terms <- length(n_seq)
  p_mat <- matrix(0, nrow = length(mu), ncol = n_terms)
  colnames(p_mat) <- n_seq

  # Seed the recursion with the diagonal associated Legendre term ==============
  if (m == 0) {
    p_mm <- rep(1, length(mu))
  } else {
    p_mm <- (-1)^m * .tmm_double_factorial_odd(m) * (1 - mu^2)^(m / 2)
  }
  p_mat[, 1] <- p_mm

  # Add the first super-diagonal term when it is retained ======================
  if (n_terms > 1) {
    p_m1m <- (2 * m + 1) * mu * p_mm
    p_mat[, 2] <- p_m1m
  }

  # March upward in degree using the three-term recursion ======================
  if (n_terms > 2) {
    for (n in (m + 2):n_max) {
      j <- n - m + 1
      p_mat[, j] <- (
        (2 * n - 1) * mu * p_mat[, j - 1] -
          (n + m - 1) * p_mat[, j - 2]
      ) / (n - m)
    }
  }

  # Return the full P_n^m(mu) block ============================================
  p_mat
}

# Differentiate the associated Legendre table with respect to theta for the
# curved-surface normal derivative.
#' @noRd
.tmm_assoc_legendre_theta_derivative <- function(m, n_seq, mu, p_mat) {
  # Preallocate the theta-derivative table for the retained block ==============
  sin_theta <- sqrt(pmax(1 - mu^2, 0))
  deriv <- matrix(0, nrow = nrow(p_mat), ncol = ncol(p_mat))

  # Apply the stable recurrence away from the polar singularities ==============
  if (all(sin_theta > 0)) {
    for (j in seq_along(n_seq)) {
      n <- n_seq[j]
      p_nm <- p_mat[, j]
      p_nm1 <- if (n == m) rep(0, length(mu)) else p_mat[, j - 1]
      deriv[, j] <- (n * mu * p_nm - (n + m) * p_nm1) / sin_theta
    }
  }

  # Return the derivative table used in the boundary operator ==================
  deriv
}

# Evaluate r(theta) and dr/dtheta for the geometries supported by the
# spherical-coordinate branch.
#' @noRd
.tmm_surface_radius <- function(shape_parameters, theta) {
  # Resolve the analytic radius law for a spherical target =====================
  if (shape_parameters$shape == "Sphere") {
    radius <- rep(as.numeric(shape_parameters$radius)[1], length(theta))
    return(list(radius = radius, radius_derivative = rep(0, length(theta))))
  }

  # Resolve the prolate spheroid radius and slope ==============================
  if (shape_parameters$shape == "ProlateSpheroid") {
    a <- as.numeric(shape_parameters$semimajor_length)[1]
    b <- as.numeric(shape_parameters$semiminor_length)[1]
    denom <- cos(theta)^2 / a^2 + sin(theta)^2 / b^2
    radius <- 1 / sqrt(denom)
    radius_derivative <- -sin(theta) * cos(theta) *
      (1 / b^2 - 1 / a^2) / (denom^(3 / 2))

    return(list(radius = radius, radius_derivative = radius_derivative))
  }

  # Resolve the oblate spheroid radius and slope ===============================
  if (shape_parameters$shape == "OblateSpheroid") {
    c_axial <- as.numeric(shape_parameters$semiminor_length)[1]
    a_equatorial <- as.numeric(shape_parameters$semimajor_length)[1]
    denom <- cos(theta)^2 / c_axial^2 + sin(theta)^2 / a_equatorial^2
    radius <- 1 / sqrt(denom)
    radius_derivative <- -sin(theta) * cos(theta) *
      (1 / a_equatorial^2 - 1 / c_axial^2) / (denom^(3 / 2))

    return(list(radius = radius, radius_derivative = radius_derivative))
  }

  # Resolve the piecewise finite-cylinder radius and slope =====================
  if (shape_parameters$shape == "Cylinder") {
    half_length <- as.numeric(shape_parameters$length)[1] / 2
    cyl_radius <- max(as.numeric(shape_parameters$radius), na.rm = TRUE)
    cos_theta <- cos(theta)
    sin_theta <- sin(theta)
    r_end <- half_length / pmax(abs(cos_theta), .Machine$double.eps)
    r_side <- cyl_radius / pmax(abs(sin_theta), .Machine$double.eps)
    radius <- pmin(r_end, r_side)

    radius_derivative <- numeric(length(theta))
    if (length(theta) > 1) {
      radius_derivative[1] <- (radius[2] - radius[1]) / (theta[2] - theta[1])
      radius_derivative[length(theta)] <- (
        radius[length(theta)] - radius[length(theta) - 1]
      ) / (theta[length(theta)] - theta[length(theta) - 1])
      if (length(theta) > 2) {
        radius_derivative[2:(length(theta) - 1)] <- (
          radius[3:length(theta)] - radius[1:(length(theta) - 2)]
        ) / (
          theta[3:length(theta)] - theta[1:(length(theta) - 2)]
        )
      }
    }

    return(list(radius = radius, radius_derivative = radius_derivative))
  }

  # Reject unsupported geometry labels =========================================
  stop("Unsupported TMM shape geometry.", call. = FALSE)
}

# Evaluate one spherical radial-function family across all retained degrees.
#' @noRd
.tmm_radial_matrix <- function(fun, n_seq, argument) {
  # Evaluate the requested radial family degree by degree ======================
  values <- lapply(n_seq, function(n) as.complex(fun(n, argument)))
  # Bind the retained modal vectors into one matrix ============================
  do.call(cbind, values)
}

# Apply the geometric normal derivative on a curved axisymmetric surface.
#' @noRd
.tmm_normal_derivative_matrix <- function(radial,
                                          radial_deriv,
                                          angular,
                                          angular_theta_deriv,
                                          k,
                                          radius,
                                          radius_derivative) {
  # Combine radial and angular derivatives into the surface normal operator ====
  k * radial_deriv * angular -
    sweep(
      radial * angular_theta_deriv,
      1,
      radius_derivative / (radius^2),
      `*`
    )
}

# Incident plane-wave coefficients for one azimuthal order in the spherical
# basis.
#' @noRd
.tmm_incident_plane_wave_coefficients <- function(m, n_seq, mu0, p_inc) {
  # Resolve the m-dependent normalization factors ==============================
  beta <- if (m == 0) {
    rep(1, length(n_seq))
  } else {
    2 * exp(lgamma(n_seq - m + 1) - lgamma(n_seq + m + 1))
  }

  # Return the incident regular-wave coefficients for this block ===============
  (1i)^n_seq * (2 * n_seq + 1) * beta * p_inc
}

# Reconstruct the monostatic far-field amplitude from the solved outgoing
# block coefficients.
#' @noRd
.tmm_backscatter_from_blocks <- function(blocks, k_sw, mu0) {
  # Initialize the accumulated monostatic far-field amplitude ==================
  f_bs <- 0 + 0i

  # Sum the outgoing coefficients from each retained azimuthal block ===========
  for (block in blocks) {
    n_seq <- block$n_seq
    p_inc <- as.numeric(block$p_inc)
    outgoing <- block$coefficients
    f_bs <- f_bs + sum(
      outgoing * ((-1)^n_seq * (-1i)^(n_seq + 1) * p_inc) / k_sw
    )
  }

  # Return the reconstructed monostatic amplitude ==============================
  f_bs
}

# Solve the projected block system with progressively more forgiving fallbacks.
#' @noRd
.tmm_solve_linear_system <- function(lhs, rhs) {
  # Use the direct solve when the projected system is square ===================
  if (nrow(lhs) == ncol(lhs)) {
    return(
      tryCatch(
        solve(lhs, rhs),
        error = function(...) {
          qr.solve(lhs, rhs)
        }
      )
    )
  }

  # Fall back to QR or normal equations for rectangular systems ================
  tryCatch(
    qr.solve(lhs, rhs),
    error = function(...) {
      lhs_h <- Conj(t(lhs))
      solve(lhs_h %*% lhs, lhs_h %*% rhs)
    }
  )
}

# Assemble the collocation system for one azimuthal block under the requested
# boundary condition.
#' @noRd
.tmm_boundary_block <- function(boundary,
                                n_seq,
                                p_mat,
                                dp_dtheta,
                                radius,
                                radius_derivative,
                                k_sw,
                                k_body = NULL,
                                rho_sw,
                                rho_body = NULL) {
  # Evaluate the exterior radial-function families on the target surface =======
  kr_sw <- k_sw * radius
  j_sw <- .tmm_radial_matrix(js, n_seq, kr_sw)
  dj_sw <- .tmm_radial_matrix(jsd, n_seq, kr_sw)
  h_sw <- .tmm_radial_matrix(hs, n_seq, kr_sw)
  dh_sw <- .tmm_radial_matrix(hsd, n_seq, kr_sw)

  reg_normal <- .tmm_normal_derivative_matrix(
    radial = j_sw,
    radial_deriv = dj_sw,
    angular = p_mat,
    angular_theta_deriv = dp_dtheta,
    k = k_sw,
    radius = radius,
    radius_derivative = radius_derivative
  )
  out_normal <- .tmm_normal_derivative_matrix(
    radial = h_sw,
    radial_deriv = dh_sw,
    angular = p_mat,
    angular_theta_deriv = dp_dtheta,
    k = k_sw,
    radius = radius,
    radius_derivative = radius_derivative
  )

  # Build the projected boundary system for rigid or soft targets ==============
  if (boundary == "pressure_release") {
    # Pressure-release targets enforce zero total pressure at the surface.
    return(list(lhs = h_sw * p_mat, rhs = -(j_sw * p_mat)))
  }

  if (boundary == "fixed_rigid") {
    # Rigid targets enforce zero normal velocity at the surface.
    return(list(lhs = out_normal, rhs = -reg_normal))
  }

  # Add the interior field block for penetrable targets ========================
  kr_body <- k_body * radius
  j_body <- .tmm_radial_matrix(js, n_seq, kr_body)
  dj_body <- .tmm_radial_matrix(jsd, n_seq, kr_body)
  in_normal <- .tmm_normal_derivative_matrix(
    radial = j_body,
    radial_deriv = dj_body,
    angular = p_mat,
    angular_theta_deriv = dp_dtheta,
    k = k_body,
    radius = radius,
    radius_derivative = radius_derivative
  )

  lhs <- rbind(
    cbind(h_sw * p_mat, -(j_body * p_mat)),
    cbind(out_normal / rho_sw, -in_normal / rho_body)
  )
  rhs <- -rbind(
    j_sw * p_mat,
    reg_normal / rho_sw
  )

  # Return the assembled collocation system ====================================
  list(lhs = lhs, rhs = rhs)
}

# Solve one frequency of the spherical-coordinate TMM branch and optionally
# retain the per-order projected blocks.
#' @noRd
.tmm_single_frequency_spherical <- function(k_sw,
                                            k_body,
                                            theta_body,
                                            boundary,
                                            shape_parameters,
                                            rho_sw,
                                            rho_body = NULL,
                                            n_max,
                                            store_t_matrix = FALSE) {
  # Initialize the per-order storage for one retained spherical solve ==========
  mu0 <- cos(theta_body)
  blocks <- vector("list", length = n_max + 1L)

  for (m in 0:n_max) {
    # Build the retained degree sequence and collocation grid ==================
    n_seq <- m:n_max
    n_terms <- length(n_seq)
    n_nodes <- .tmm_collocation_nodes(shape_parameters, boundary, n_terms)
    quad <- gauss_legendre(n_nodes, a = -1, b = 1)
    mu <- quad$nodes
    theta <- acos(mu)
    surface <- .tmm_surface_radius(shape_parameters, theta)
    p_mat <- .tmm_assoc_legendre_table(m, n_max, mu)
    dp_dtheta <- .tmm_assoc_legendre_theta_derivative(m, n_seq, mu, p_mat)

    # Assemble the collocation system for this azimuthal order =================
    system <- .tmm_boundary_block(
      boundary = boundary,
      n_seq = n_seq,
      p_mat = p_mat,
      dp_dtheta = dp_dtheta,
      radius = surface$radius,
      radius_derivative = surface$radius_derivative,
      k_sw = k_sw,
      k_body = k_body,
      rho_sw = rho_sw,
      rho_body = rho_body
    )

    # Project the collocation equations back onto the retained modal basis =====
    surface_weight <- surface$radius * sqrt(surface$radius^2 +
      surface$radius_derivative^2)
    weighted_test <- sweep(p_mat, 1, quad$weights * surface_weight, `*`)
    if (boundary %in% c("fixed_rigid", "pressure_release")) {
      lhs_proj <- Conj(t(weighted_test)) %*% system$lhs
      rhs_proj <- Conj(t(weighted_test)) %*% system$rhs
    } else {
      n_rows <- nrow(p_mat)
      projector <- rbind(
        cbind(Conj(t(weighted_test)), matrix(0, n_terms, n_rows)),
        cbind(matrix(0, n_terms, n_rows), Conj(t(weighted_test)))
      )
      lhs_proj <- projector %*% system$lhs
      rhs_proj <- projector %*% system$rhs
    }

    # Keep a lightweight conditioning indicator for diagnostics ================
    rcond_lhs <- tryCatch(
      rcond(lhs_proj),
      error = function(...) NA_real_
    )

    # Solve the projected system and recover the outgoing block ================
    solution <- .tmm_solve_linear_system(lhs_proj, rhs_proj)
    t_block <- if (boundary %in% c("fixed_rigid", "pressure_release")) {
      solution
    } else {
      solution[seq_len(n_terms), , drop = FALSE]
    }

    # Apply the block T-matrix to the incident coefficients ====================
    p_inc <- as.numeric(.tmm_assoc_legendre_table(m, n_max, mu0))
    a_inc <- .tmm_incident_plane_wave_coefficients(m, n_seq, mu0, p_inc)
    coeffs <- as.vector(t_block %*% a_inc)

    # Store the solved block and its diagnostic bookkeeping ====================
    blocks[[m + 1L]] <- list(
      m = m,
      n_seq = n_seq,
      p_inc = p_inc,
      coefficients = coeffs,
      rcond_lhs = rcond_lhs,
      T = if (isTRUE(store_t_matrix)) t_block else NULL
    )
  }

  # Return the monostatic amplitude and the retained block state ===============
  list(
    f_bs = .tmm_backscatter_from_blocks(blocks, k_sw = k_sw, mu0 = mu0),
    blocks = blocks
  )
}
