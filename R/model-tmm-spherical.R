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

# Build the collocation quadrature used by one spherical-coordinate TMM block.
# Cylinders are piecewise-smooth, so the retained axisymmetric cylinder backend splits the
# Gauss-Legendre nodes at the cap/sidewall transition instead of integrating
# across the derivative kink with one smooth panel.
#' @noRd
.tmm_collocation_quadrature <- function(shape_parameters,
                                        boundary,
                                        n_terms,
                                        collocation_multiplier = NULL) {
  n_nodes <- .tmm_collocation_nodes(shape_parameters, boundary, n_terms)
  if (!is.null(collocation_multiplier) &&
    is.finite(collocation_multiplier) &&
    collocation_multiplier > 0) {
    n_nodes <- max(n_nodes, as.integer(ceiling(collocation_multiplier * n_terms)))
  }
  if (!identical(.tmm_shape_name(shape_parameters), "Cylinder")) {
    return(gauss_legendre(n_nodes, a = -1, b = 1))
  }

  half_length <- as.numeric(shape_parameters$length)[1] / 2
  cyl_radius <- max(as.numeric(shape_parameters$radius), na.rm = TRUE)
  mu_switch <- half_length / sqrt(half_length^2 + cyl_radius^2)
  if (!is.finite(mu_switch) || mu_switch <= 0 || mu_switch >= 1) {
    return(gauss_legendre(n_nodes, a = -1, b = 1))
  }

  n_cap <- max(32L, ceiling(0.20 * n_nodes))
  n_side <- max(64L, n_nodes - 2L * n_cap)
  remaining <- n_nodes - n_side
  n_cap_lo <- as.integer(floor(remaining / 2L))
  n_cap_hi <- as.integer(ceiling(remaining / 2L))
  n_cap_lo <- max(24L, n_cap_lo)
  n_cap_hi <- max(24L, n_cap_hi)

  quad_cap_hi <- gauss_legendre(n_cap_hi, a = mu_switch, b = 1)
  quad_side <- gauss_legendre(n_side, a = -mu_switch, b = mu_switch)
  quad_cap_lo <- gauss_legendre(n_cap_lo, a = -1, b = -mu_switch)

  list(
    nodes = c(
      quad_cap_hi$nodes,
      quad_side$nodes,
      quad_cap_lo$nodes
    ),
    weights = c(
      quad_cap_hi$weights,
      quad_side$weights,
      quad_cap_lo$weights
    )
  )
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
.tmm_cylinder_endcap_fraction <- function(shape_parameters,
                                          cylinder_endcap_fraction = NULL) {
  taper_order <- if ("taper_order" %in% names(shape_parameters)) {
    as.numeric(shape_parameters$taper_order)[1]
  } else {
    NA_real_
  }

  if (!is.null(cylinder_endcap_fraction)) {
    return(.tmm_resolve_cylinder_endcap_fraction(cylinder_endcap_fraction))
  }

  if (is.finite(taper_order)) {
    return(0)
  }

  0
}

# Evaluate a short spheroidal-endcap cylinder radius law used to regularize the
# sharp corner in the retained axisymmetric cylinder branch.
#' @noRd
.tmm_cylinder_spheroidal_radius <- function(theta,
                                            half_length,
                                            cyl_radius,
                                            endcap_fraction) {
  cap_length <- min(half_length * endcap_fraction, half_length - 1e-9)
  if (!is.finite(cap_length) || cap_length <= 0) {
    return(NULL)
  }

  z_cap_center <- half_length - cap_length
  theta_reduced <- pmin(theta, pi - theta)
  radius <- vapply(
    theta_reduced,
    function(th) {
      if (!is.finite(th) || th <= sqrt(.Machine$double.eps)) {
        return(half_length)
      }
      if (abs(th - pi / 2) <= sqrt(.Machine$double.eps)) {
        return(cyl_radius)
      }

      cos_theta <- cos(th)
      sin_theta <- sin(th)
      side_radius <- cyl_radius / sin_theta
      side_height <- side_radius * cos_theta
      if (side_height <= z_cap_center + 1e-10) {
        return(side_radius)
      }

      quad_a <- (cos_theta^2) / (cap_length^2) +
        (sin_theta^2) / (cyl_radius^2)
      quad_b <- -2 * z_cap_center * cos_theta / (cap_length^2)
      quad_c <- (z_cap_center^2) / (cap_length^2) - 1
      disc <- pmax(quad_b^2 - 4 * quad_a * quad_c, 0)
      root_pos <- (-quad_b + sqrt(disc)) / (2 * quad_a)
      root_neg <- (-quad_b - sqrt(disc)) / (2 * quad_a)
      max(root_pos, root_neg)
    },
    numeric(1)
  )

  radius_derivative <- numeric(length(theta))
  if (length(theta) > 1) {
    radius_derivative[1] <- (radius[2] - radius[1]) /
      (theta[2] - theta[1])
    radius_derivative[length(theta)] <- (
      radius[length(theta)] - radius[length(theta) - 1]
    ) / (theta[length(theta)] - theta[length(theta) - 1])
  }
  if (length(theta) > 2) {
    radius_derivative[2:(length(theta) - 1)] <- (
      radius[3:length(theta)] - radius[1:(length(theta) - 2)]
    ) / (
      theta[3:length(theta)] - theta[1:(length(theta) - 2)]
    )
  }

  list(radius = radius, radius_derivative = radius_derivative)
}

# Evaluate r(theta) and dr/dtheta for the geometries supported by the
# spherical-coordinate branch.
#' @noRd
.tmm_surface_radius <- function(shape_parameters,
                                theta,
                                cylinder_endcap_fraction = NULL) {
  shape_source <- shape_parameters$shell %||% shape_parameters
  shape_name <- .tmm_shape_name(shape_parameters)
  # Resolve the analytic radius law for a spherical target =====================
  if (identical(shape_name, "Sphere")) {
    radius <- rep(as.numeric(shape_source$radius)[1], length(theta))
    return(list(radius = radius, radius_derivative = rep(0, length(theta))))
  }

  # Resolve the prolate spheroid radius and slope ==============================
  if (identical(shape_name, "ProlateSpheroid")) {
    a <- as.numeric(shape_source$semimajor_length)[1]
    b <- as.numeric(shape_source$semiminor_length)[1]
    denom <- cos(theta)^2 / a^2 + sin(theta)^2 / b^2
    radius <- 1 / sqrt(denom)
    radius_derivative <- -sin(theta) * cos(theta) *
      (1 / b^2 - 1 / a^2) / (denom^(3 / 2))

    return(list(radius = radius, radius_derivative = radius_derivative))
  }

  # Resolve the oblate spheroid radius and slope ===============================
  if (identical(shape_name, "OblateSpheroid")) {
    c_axial <- as.numeric(shape_source$semiminor_length)[1]
    a_equatorial <- as.numeric(shape_source$semimajor_length)[1]
    denom <- cos(theta)^2 / c_axial^2 + sin(theta)^2 / a_equatorial^2
    radius <- 1 / sqrt(denom)
    radius_derivative <- -sin(theta) * cos(theta) *
      (1 / a_equatorial^2 - 1 / c_axial^2) / (denom^(3 / 2))

    return(list(radius = radius, radius_derivative = radius_derivative))
  }

  # Resolve the piecewise finite-cylinder radius and slope =====================
  if (identical(shape_name, "Cylinder")) {
    half_length <- as.numeric(shape_source$length)[1] / 2
    cyl_radius <- max(as.numeric(shape_source$radius), na.rm = TRUE)
    taper_order <- if ("taper_order" %in% names(shape_parameters)) {
      as.numeric(shape_parameters$taper_order)[1]
    } else {
      NA_real_
    }

    endcap_fraction <- .tmm_cylinder_endcap_fraction(
      shape_parameters = shape_parameters,
      cylinder_endcap_fraction = cylinder_endcap_fraction
    )
    if (endcap_fraction > 0) {
      endcap_profile <- .tmm_cylinder_spheroidal_radius(
        theta = theta,
        half_length = half_length,
        cyl_radius = cyl_radius,
        endcap_fraction = endcap_fraction
      )
      if (!is.null(endcap_profile)) {
        return(endcap_profile)
      }
    }

    if (is.finite(taper_order)) {
      # A tapered cylinder stays axisymmetric but smooths the cap/sidewall
      # corner. Solve the ray/profile intersection numerically and recover the
      # theta-derivative by finite differences on the collocation grid.
      taper_order <- max(2, abs(taper_order))
      theta_reduced <- pmin(theta, pi - theta)
      x_hit <- vapply(
        theta_reduced,
        function(th) {
          if (!is.finite(th) || th <= sqrt(.Machine$double.eps)) {
            return(half_length)
          }
          if (abs(th - pi / 2) <= sqrt(.Machine$double.eps)) {
            return(0)
          }

          root_fun <- function(x) {
            profile_ratio <- pmax(1 - (abs(x) / half_length)^taper_order, 0)
            x * tan(th) - cyl_radius * sqrt(profile_ratio)
          }

          uniroot(
            root_fun,
            lower = 0,
            upper = half_length,
            tol = 1e-10
          )$root
        },
        numeric(1)
      )
      radius <- x_hit / pmax(abs(cos(theta_reduced)), .Machine$double.eps)
      radius[abs(theta_reduced - pi / 2) <= sqrt(.Machine$double.eps)] <- cyl_radius

      radius_derivative <- numeric(length(theta))
      if (length(theta) > 1) {
        radius_derivative[1] <- (radius[2] - radius[1]) /
          (theta[2] - theta[1])
        radius_derivative[length(theta)] <- (
          radius[length(theta)] - radius[length(theta) - 1]
        ) / (theta[length(theta)] - theta[length(theta) - 1])
      }
      if (length(theta) > 2) {
        radius_derivative[2:(length(theta) - 1)] <- (
          radius[3:length(theta)] - radius[1:(length(theta) - 2)]
        ) / (
          theta[3:length(theta)] - theta[1:(length(theta) - 2)]
        )
      }

      return(list(radius = radius, radius_derivative = radius_derivative))
    }

    theta_switch <- atan2(cyl_radius, half_length)
    radius <- numeric(length(theta))
    radius_derivative <- numeric(length(theta))

    cap_lo <- theta <= theta_switch + 1e-12
    cap_hi <- theta >= (pi - theta_switch - 1e-12)
    side <- !(cap_lo | cap_hi)

    if (any(cap_lo)) {
      th <- theta[cap_lo]
      radius[cap_lo] <- half_length / pmax(cos(th), .Machine$double.eps)
      radius_derivative[cap_lo] <- half_length * sin(th) /
        pmax(cos(th)^2, .Machine$double.eps)
    }
    if (any(cap_hi)) {
      th <- theta[cap_hi]
      radius[cap_hi] <- -half_length / pmin(cos(th), -.Machine$double.eps)
      radius_derivative[cap_hi] <- -half_length * sin(th) /
        pmax(cos(th)^2, .Machine$double.eps)
    }
    if (any(side)) {
      th <- theta[side]
      radius[side] <- cyl_radius / pmax(sin(th), .Machine$double.eps)
      radius_derivative[side] <- -cyl_radius * cos(th) /
        pmax(sin(th)^2, .Machine$double.eps)
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

# Assemble piecewise cap/side quadrature for a sharp finite cylinder. Using the
# true generator segments keeps the cap/side junction out of the surface
# parameterization and avoids forcing the sharp edge through a smooth r(theta)
# representation.
#' @noRd
.tmm_cylinder_piecewise_surface <- function(shape_parameters, n_terms) {
  half_length <- as.numeric(shape_parameters$length)[1] / 2
  cyl_radius <- max(as.numeric(shape_parameters$radius), na.rm = TRUE)

  n_side <- max(128L, 10L * n_terms)
  n_cap <- max(96L, 8L * n_terms)

  quad_side <- gauss_legendre(n_side, a = -half_length, b = half_length)
  quad_cap <- gauss_legendre(n_cap, a = 0, b = cyl_radius)

  z_side <- quad_side$nodes
  rho_cap <- quad_cap$nodes

  r_side <- sqrt(cyl_radius^2 + z_side^2)
  mu_side <- z_side / r_side
  theta_side <- acos(mu_side)
  weight_side <- cyl_radius * quad_side$weights
  normal_r_side <- sin(theta_side)
  normal_theta_side <- cos(theta_side)

  r_top <- sqrt(rho_cap^2 + half_length^2)
  mu_top <- half_length / r_top
  theta_top <- acos(mu_top)
  weight_top <- rho_cap * quad_cap$weights
  normal_r_top <- cos(theta_top)
  normal_theta_top <- -sin(theta_top)

  r_bottom <- r_top
  mu_bottom <- -mu_top
  theta_bottom <- acos(mu_bottom)
  weight_bottom <- weight_top
  normal_r_bottom <- -cos(theta_bottom)
  normal_theta_bottom <- sin(theta_bottom)

  list(
    mu = c(mu_top, mu_side, mu_bottom),
    theta = c(theta_top, theta_side, theta_bottom),
    radius = c(r_top, r_side, r_bottom),
    area_weight = c(weight_top, weight_side, weight_bottom),
    normal_r = c(normal_r_top, normal_r_side, normal_r_bottom),
    normal_theta = c(normal_theta_top, normal_theta_side, normal_theta_bottom)
  )
}

# Apply the geometric normal derivative on one explicitly parameterized
# axisymmetric surface. This is used for sharp finite cylinders where the cap
# and sidewall normals are known analytically on each patch.
#' @noRd
.tmm_normal_derivative_explicit <- function(radial,
                                            radial_deriv,
                                            angular,
                                            angular_theta_deriv,
                                            k,
                                            radius,
                                            normal_r,
                                            normal_theta) {
  sweep(k * radial_deriv * angular, 1, normal_r, `*`) +
    sweep((radial * angular_theta_deriv) / radius, 1, normal_theta, `*`)
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

# Solve one projected cylinder block after explicit row/column equilibration.
# The sharp-cylinder retained systems can become catastrophically ill-
# conditioned once the projected Legendre basis grows, so a scaled SVD solve is
# materially more stable than a raw dense solve.
#' @noRd
.tmm_solve_linear_system_stabilized <- function(lhs,
                                                rhs,
                                                svd_rel_tol = 1e-10) {
  row_scale <- sqrt(rowSums(Mod(lhs)^2))
  row_scale[!is.finite(row_scale) | row_scale <= 0] <- 1
  lhs_row <- lhs / row_scale
  rhs_row <- rhs / row_scale

  col_scale <- sqrt(colSums(Mod(lhs_row)^2))
  col_scale[!is.finite(col_scale) | col_scale <= 0] <- 1
  lhs_scaled <- sweep(lhs_row, 2, col_scale, `/`)
  if (any(!is.finite(lhs_scaled)) || any(!is.finite(rhs_row))) {
    return(.tmm_solve_linear_system(lhs, rhs))
  }

  svd_fit <- svd(lhs_scaled)
  if (!length(svd_fit$d) || !all(is.finite(svd_fit$d))) {
    return(.tmm_solve_linear_system(lhs, rhs))
  }

  max_sv <- max(svd_fit$d)
  keep <- svd_fit$d > max_sv * svd_rel_tol
  if (!any(keep)) {
    keep[[1]] <- TRUE
  }

  u_keep <- svd_fit$u[, keep, drop = FALSE]
  v_keep <- svd_fit$v[, keep, drop = FALSE]
  d_keep <- svd_fit$d[keep]

  coeff_scaled <- v_keep %*% (
    diag(1 / d_keep, nrow = length(d_keep)) %*% (
      Conj(t(u_keep)) %*% rhs_row
    )
  )

  sweep(coeff_scaled, 1, col_scale, `/`)
}

# Compute the principal curvatures of a prolate spheroid of revolution at the
# supplied polar angles. The retained elastic-shell prolate branch uses these
# curvatures in a thin-shell surface-impedance approximation.
#' @noRd
.tmm_prolate_principal_curvatures <- function(shape_parameters, theta) {
  shape_source <- shape_parameters$shell %||% shape_parameters
  a <- as.numeric(shape_source$semimajor_length)[1]
  b <- as.numeric(shape_source$semiminor_length)[1]
  denom <- a^2 * sin(theta)^2 + b^2 * cos(theta)^2

  list(
    meridional = a * b / (denom^(3 / 2)),
    circumferential = a / (b * sqrt(denom))
  )
}

# Build the thin-shell surface-density jump coefficient used by the elastic
# shell projected TMM branch. The coefficient multiplies the normal-derivative
# continuity variable, so it has units of surface density.
#' @noRd
.tmm_elastic_shell_alpha <- function(shape_parameters,
                                     shell_body,
                                     theta,
                                     frequency_hz) {
  if (!identical(.tmm_shape_name(shape_parameters), "ProlateSpheroid")) {
    stop(
      "Elastic-shelled projected TMM currently supports prolate spheroids only.",
      call. = FALSE
    )
  }

  if (is.null(shell_body)) {
    stop(
      "Elastic-shelled projected TMM requires explicit shell-body parameters.",
      call. = FALSE
    )
  }

  shell_thickness <- as.numeric(shell_body$shell_thickness)[1]
  shell_density <- as.numeric(shell_body$shell_density)[1]
  shell_E <- as.numeric(shell_body$shell_E)[1]
  shell_nu <- as.numeric(shell_body$shell_nu)[1]
  if (!all(is.finite(c(shell_thickness, shell_density, shell_E, shell_nu)))) {
    stop(
      "Elastic-shelled projected TMM requires finite shell thickness, density, ",
      "Young's modulus, and Poisson ratio.",
      call. = FALSE
    )
  }

  omega <- 2 * pi * as.numeric(frequency_hz)[1]
  if (!is.finite(omega) || omega <= 0) {
    stop(
      "'frequency_hz' must be one positive finite value for the elastic shell branch.",
      call. = FALSE
    )
  }

  curvatures <- .tmm_prolate_principal_curvatures(shape_parameters, theta)
  kappa_m <- curvatures$meridional
  kappa_c <- curvatures$circumferential
  shell_extensional <- shell_E * shell_thickness / (1 - shell_nu^2)
  shell_bending <- shell_E * shell_thickness^3 / (12 * (1 - shell_nu^2))
  membrane_stiffness <- shell_extensional * (
    kappa_m^2 + kappa_c^2 + 2 * shell_nu * kappa_m * kappa_c
  )
  bending_stiffness <- shell_bending * (kappa_m^2 + kappa_c^2)^2

  (membrane_stiffness + bending_stiffness) / (omega^2) -
    shell_density * shell_thickness
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
                                rho_body = NULL,
                                shell_alpha = NULL) {
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

  if (identical(boundary, "elastic_shelled")) {
    lhs <- rbind(
      cbind(
        h_sw * p_mat - sweep(out_normal / rho_sw, 1, shell_alpha, `*`),
        -(j_body * p_mat)
      ),
      cbind(out_normal / rho_sw, -in_normal / rho_body)
    )
    rhs <- -rbind(
      j_sw * p_mat - sweep(reg_normal / rho_sw, 1, shell_alpha, `*`),
      reg_normal / rho_sw
    )

    return(list(lhs = lhs, rhs = rhs))
  }

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

# Assemble the piecewise sharp-cylinder system for one azimuthal block. This
# uses the explicit cap and sidewall normals instead of the smooth star-shaped
# derivative identity used by the generic spherical branch.
#' @noRd
.tmm_boundary_block_cylinder_piecewise <- function(boundary,
                                                   n_seq,
                                                   p_mat,
                                                   dp_dtheta,
                                                   radius,
                                                   normal_r,
                                                   normal_theta,
                                                   k_sw,
                                                   k_body = NULL,
                                                   rho_sw,
                                                   rho_body = NULL) {
  kr_sw <- k_sw * radius
  j_sw <- .tmm_radial_matrix(js, n_seq, kr_sw)
  dj_sw <- .tmm_radial_matrix(jsd, n_seq, kr_sw)
  h_sw <- .tmm_radial_matrix(hs, n_seq, kr_sw)
  dh_sw <- .tmm_radial_matrix(hsd, n_seq, kr_sw)

  reg_normal <- .tmm_normal_derivative_explicit(
    radial = j_sw,
    radial_deriv = dj_sw,
    angular = p_mat,
    angular_theta_deriv = dp_dtheta,
    k = k_sw,
    radius = radius,
    normal_r = normal_r,
    normal_theta = normal_theta
  )
  out_normal <- .tmm_normal_derivative_explicit(
    radial = h_sw,
    radial_deriv = dh_sw,
    angular = p_mat,
    angular_theta_deriv = dp_dtheta,
    k = k_sw,
    radius = radius,
    normal_r = normal_r,
    normal_theta = normal_theta
  )

  if (boundary == "pressure_release") {
    return(list(lhs = h_sw * p_mat, rhs = -(j_sw * p_mat)))
  }

  if (boundary == "fixed_rigid") {
    return(list(lhs = out_normal, rhs = -reg_normal))
  }

  kr_body <- k_body * radius
  j_body <- .tmm_radial_matrix(js, n_seq, kr_body)
  dj_body <- .tmm_radial_matrix(jsd, n_seq, kr_body)
  in_normal <- .tmm_normal_derivative_explicit(
    radial = j_body,
    radial_deriv = dj_body,
    angular = p_mat,
    angular_theta_deriv = dp_dtheta,
    k = k_body,
    radius = radius,
    normal_r = normal_r,
    normal_theta = normal_theta
  )

  lhs <- rbind(
    cbind(h_sw * p_mat, -(j_body * p_mat)),
    cbind(out_normal / rho_sw, -in_normal / rho_body)
  )
  rhs <- -rbind(
    j_sw * p_mat,
    reg_normal / rho_sw
  )

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
                                            cylinder_endcap_fraction = NULL,
                                            store_t_matrix = FALSE,
                                            shell_body = NULL,
                                            frequency_hz = NULL,
                                            collocation_multiplier = NULL,
                                            svd_rel_tol = NULL) {
  # Initialize the per-order storage for one retained spherical solve ==========
  mu0 <- cos(theta_body)
  blocks <- vector("list", length = n_max + 1L)
  use_piecewise_cylinder <- identical(.tmm_shape_name(shape_parameters), "Cylinder") &&
    .tmm_cylinder_endcap_fraction(
      shape_parameters = shape_parameters,
      cylinder_endcap_fraction = cylinder_endcap_fraction
    ) <= 0 &&
    !("taper_order" %in% names(shape_parameters) &&
      is.finite(as.numeric(shape_parameters$taper_order)[1]))

  for (m in 0:n_max) {
    # Build the retained degree sequence and collocation grid ==================
    n_seq <- m:n_max
    n_terms <- length(n_seq)
    if (use_piecewise_cylinder) {
      surface <- .tmm_cylinder_piecewise_surface(shape_parameters, n_terms)
      mu <- surface$mu
      p_mat <- .tmm_assoc_legendre_table(m, n_max, mu)
      dp_dtheta <- .tmm_assoc_legendre_theta_derivative(m, n_seq, mu, p_mat)
      system <- .tmm_boundary_block_cylinder_piecewise(
        boundary = boundary,
        n_seq = n_seq,
        p_mat = p_mat,
        dp_dtheta = dp_dtheta,
        radius = surface$radius,
        normal_r = surface$normal_r,
        normal_theta = surface$normal_theta,
        k_sw = k_sw,
        k_body = k_body,
        rho_sw = rho_sw,
        rho_body = rho_body
      )
      weighted_test <- sweep(p_mat, 1, surface$area_weight, `*`)
    } else {
      quad <- .tmm_collocation_quadrature(
        shape_parameters = shape_parameters,
        boundary = boundary,
        n_terms = n_terms,
        collocation_multiplier = collocation_multiplier
      )
      mu <- quad$nodes
      theta <- acos(mu)
      surface <- .tmm_surface_radius(
        shape_parameters = shape_parameters,
        theta = theta,
        cylinder_endcap_fraction = cylinder_endcap_fraction
      )
      p_mat <- .tmm_assoc_legendre_table(m, n_max, mu)
      dp_dtheta <- .tmm_assoc_legendre_theta_derivative(m, n_seq, mu, p_mat)

      # Assemble the collocation system for this azimuthal order ===============
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
        rho_body = rho_body,
        shell_alpha = if (identical(boundary, "elastic_shelled")) {
          .tmm_elastic_shell_alpha(
            shape_parameters = shape_parameters,
            shell_body = shell_body,
            theta = theta,
            frequency_hz = frequency_hz
          )
        } else {
          NULL
        }
      )

      # Project the collocation equations back onto the retained modal basis ===
      surface_weight <- surface$radius * sqrt(surface$radius^2 +
        surface$radius_derivative^2)
      weighted_test <- sweep(p_mat, 1, quad$weights * surface_weight, `*`)
    }

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

    p_inc <- as.numeric(.tmm_assoc_legendre_table(m, n_max, mu0))
    a_inc <- .tmm_incident_plane_wave_coefficients(m, n_seq, mu0, p_inc)
    # Keep a lightweight conditioning indicator for diagnostics ================
    rcond_lhs <- tryCatch(
      rcond(lhs_proj),
      error = function(...) NA_real_
    )

    # Solve the projected system and recover the outgoing block ================
    use_stabilized_solver <- use_piecewise_cylinder ||
      identical(boundary, "elastic_shelled")
    if (isTRUE(store_t_matrix)) {
      solution <- if (use_stabilized_solver) {
        .tmm_solve_linear_system_stabilized(
          lhs_proj,
          rhs_proj,
          svd_rel_tol = svd_rel_tol %||% 1e-10
        )
      } else {
        .tmm_solve_linear_system(lhs_proj, rhs_proj)
      }
      t_block <- if (boundary %in% c("fixed_rigid", "pressure_release")) {
        solution
      } else {
        solution[seq_len(n_terms), , drop = FALSE]
      }
      coeffs <- as.vector(t_block %*% a_inc)
    } else {
      rhs_excited <- as.vector(rhs_proj %*% a_inc)
      solution <- if (use_stabilized_solver) {
        .tmm_solve_linear_system_stabilized(
          lhs_proj,
          rhs_excited,
          svd_rel_tol = svd_rel_tol %||% 1e-10
        )
      } else {
        .tmm_solve_linear_system(lhs_proj, rhs_excited)
      }
      coeffs <- if (boundary %in% c("fixed_rigid", "pressure_release")) {
        as.vector(solution)
      } else {
        as.vector(solution[seq_len(n_terms)])
      }
      t_block <- NULL
    }

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
