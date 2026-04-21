################################################################################
# Transition matrix method (TMM) orientation helpers
################################################################################

# Normalize one user-facing orientation distribution into the shared
# theta/phi/weight data frame used by the TMM post-processing helpers.
#' @noRd
.tmm_validate_orientation_distribution <- function(distribution) {
  # Confirm the helper output carries the expected distribution class ==========
  if (!inherits(distribution, "TMMOrientationDistribution")) {
    stop(
      "'distribution' must be created by 'tmm_orientation_distribution()'.",
      call. = FALSE
    )
  }

  # Check that the required angle and weight columns are present ===============
  required_cols <- c("theta_body", "phi_body", "weights", "distribution")
  missing_cols <- setdiff(required_cols, names(distribution))
  if (length(missing_cols)) {
    stop(
      "'distribution' is missing required column(s): ",
      paste(missing_cols, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  # Validate the numeric angle and weight columns ==============================
  if (!is.numeric(distribution$theta_body) ||
    any(!is.finite(distribution$theta_body))) {
    stop(
      "'distribution$theta_body' must be finite numeric angles in radians.",
      call. = FALSE
    )
  }
  if (!is.numeric(distribution$phi_body) ||
    any(!is.finite(distribution$phi_body))) {
    stop(
      "'distribution$phi_body' must be finite numeric angles in radians.",
      call. = FALSE
    )
  }
  if (!is.numeric(distribution$weights) ||
    any(!is.finite(distribution$weights)) ||
    any(distribution$weights < 0) || sum(distribution$weights) <= 0) {
    stop(
      "'distribution$weights' must be finite, non-negative, and sum to a ",
      "positive value.",
      call. = FALSE
    )
  }

  # Normalize the weights onto a proper averaging distribution =================
  distribution$weights <- distribution$weights / sum(distribution$weights)
  distribution
}

# Convert one theta-grid density into normalized quadrature weights.
#' @noRd
.tmm_distribution_weights <- function(theta_body,
                                      density_values,
                                      lower = NULL,
                                      upper = NULL,
                                      name = "theta_body") {
  # Validate the supplied density values against the theta grid ================
  if (!is.numeric(density_values) ||
    length(density_values) != length(theta_body) ||
    any(!is.finite(density_values)) || any(density_values < 0)) {
    stop(
      "'density_values' must be a non-negative numeric vector with the ",
      "same length as '",
      name,
      "'.",
      call. = FALSE
    )
  }

  # Convert the density into cell-integrated quadrature weights ================
  interval_weights <- .tmm_interval_weights(
    theta_body,
    lower = lower,
    upper = upper,
    name = name
  )
  weights <- density_values * interval_weights
  if (sum(weights) <= 0) {
    stop("Orientation weights must sum to a positive value.", call. = FALSE)
  }

  # Normalize the weights before returning them ================================
  weights / sum(weights)
}
#' Build an orientation distribution for stored TMM post-processing
#'
#' @description
#' Creates a validated set of incident angles and normalized weights that can be
#' reused by \code{\link{tmm_average_orientation}} or
#' \code{\link{tmm_products}}. The distributions
#' defined here are distributions in `theta_body` itself rather than isotropic
#' solid-angle distributions.
#'
#' @param distribution Orientation-distribution type. One of `"uniform"`,
#'   `"normal"`, `"truncated_normal"`, `"quadrature"`, or `"pdf"`.
#' @param theta_body Optional numeric vector of incident polar angles (radians).
#'   Required for the `"quadrature"` and `"pdf"` pathways.
#' @param weights Optional numeric quadrature weights paired with `theta_body`
#'   for `distribution = "quadrature"`.
#' @param pdf Optional user-supplied density over `theta_body` for
#'   `distribution = "pdf"`. This can be either a numeric vector the same
#'   length as `theta_body` or a function evaluated at `theta_body`.
#' @param phi_body Incident azimuth angle(s) (radians). Either scalar or the
#'   same length as the resolved `theta_body` grid.
#' @param mean_theta Mean angle (radians) for the normal-family distributions.
#' @param sd_theta Standard deviation (radians) for the normal-family
#'   distributions.
#' @param lower Lower bound (radians) for the uniform and truncated-normal
#'   distributions.
#' @param upper Upper bound (radians) for the uniform and truncated-normal
#'   distributions.
#' @param n_theta Number of grid points for the analytic distributions.
#'
#' @return A data frame with normalized orientation weights and class
#'   `"TMMOrientationDistribution"`.
#'
#' @seealso \code{\link{tmm_average_orientation}},
#'   \code{\link{tmm_products}}
#' @export
tmm_orientation_distribution <- function(distribution = c(
                                           "uniform",
                                           "normal",
                                           "truncated_normal",
                                           "quadrature",
                                           "pdf"
                                         ),
                                         theta_body = NULL,
                                         weights = NULL,
                                         pdf = NULL,
                                         phi_body = pi,
                                         mean_theta = pi / 2,
                                         sd_theta = pi / 12,
                                         lower = 0,
                                         upper = pi,
                                         n_theta = 91) {
  # Resolve and validate the requested distribution family =====================
  distribution <- match.arg(distribution)
  n_theta <- .tmm_validate_orientation_n_theta(n_theta)

  # Build the requested theta grid and associated normalized weights ===========
  orientation <- switch(distribution,
    quadrature = {
      if (is.null(theta_body)) {
        stop(
          "'theta_body' must be supplied for ",
          "'distribution = \"quadrature\"'.",
          call. = FALSE
        )
      }
      .tmm_orientation_quadrature(theta_body, weights)
    },
    pdf = {
      if (is.null(theta_body)) {
        stop(
          "'theta_body' must be supplied for ",
          "'distribution = \"pdf\"'.",
          call. = FALSE
        )
      }
      .tmm_orientation_pdf(theta_body, pdf)
    },
    uniform = .tmm_orientation_uniform(lower, upper, n_theta),
    normal = .tmm_orientation_normal_density(mean_theta, sd_theta, n_theta),
    truncated_normal = .tmm_orientation_truncated_normal(
      mean_theta = mean_theta,
      sd_theta = sd_theta,
      lower = lower,
      upper = upper,
      n_theta = n_theta
    )
  )

  theta_body <- orientation$theta_body
  weights <- orientation$weights

  # Expand the azimuth input to the resolved theta grid ========================
  phi_body <- .tmm_resolve_orientation_phi(phi_body, length(theta_body))

  # Return the standardized orientation-distribution data frame ================
  distribution_df <- data.frame(
    theta_body = theta_body,
    phi_body = phi_body,
    weights = weights,
    distribution = distribution
  )
  class(distribution_df) <- c(
    "TMMOrientationDistribution",
    class(distribution_df)
  )
  distribution_df
}

# Resolve either a pre-built orientation distribution or the direct averaging
# inputs supplied to `tmm_average_orientation()`.
#' @noRd
.tmm_average_orientation_inputs <- function(distribution,
                                            theta_body,
                                            weights,
                                            phi_body,
                                            default_phi_body = pi) {
  # Reuse the shared distribution validator when a pre-built grid is supplied ==
  if (!is.null(distribution)) {
    distribution <- .tmm_validate_orientation_distribution(distribution)
    return(list(
      theta_body = distribution$theta_body,
      phi_body = distribution$phi_body,
      weights = distribution$weights
    ))
  }

  # Otherwise validate the direct theta/phi/weight inputs ======================
  .tmm_average_orientation_direct_inputs(
    theta_body = theta_body,
    weights = weights,
    phi_body = phi_body,
    default_phi_body = default_phi_body
  )
}

# Validate direct theta/phi/weight inputs for `tmm_average_orientation()`.
#' @noRd
.tmm_average_orientation_direct_inputs <- function(theta_body,
                                                   weights,
                                                   phi_body,
                                                   default_phi_body = pi) {
  # Require a concrete angle vector when no stored distribution is supplied ====
  if (!is.numeric(theta_body) ||
    !length(theta_body) ||
    any(!is.finite(theta_body))) {
    stop(
      "'theta_body' must be a non-empty numeric vector of angles in radians.",
      call. = FALSE
    )
  }

  n_angles <- length(theta_body)
  if (is.null(phi_body)) {
    phi_body <- default_phi_body
  }
  phi_body <- rep_len(phi_body, n_angles)
  if (any(!is.finite(phi_body))) {
    stop("'phi_body' must be finite.", call. = FALSE)
  }

  # Resolve equal weights or normalize the user-supplied weight vector =========
  weights <- .tmm_average_orientation_weights(weights, n_angles, theta_body)

  list(
    theta_body = theta_body,
    phi_body = phi_body,
    weights = weights
  )
}

# Normalize the averaging weights for `tmm_average_orientation()`.
#' @noRd
.tmm_average_orientation_weights <- function(weights, n_angles, theta_body) {
  # Default to equal weighting across the supplied angle grid ==================
  if (is.null(weights)) {
    return(rep(1 / n_angles, n_angles))
  }

  # Validate and normalize explicit weights ====================================
  if (!is.numeric(weights) ||
    length(weights) != n_angles ||
    any(!is.finite(weights)) ||
    any(weights < 0) ||
    sum(weights) <= 0) {
    stop(
      "'weights' must be a non-negative numeric vector with the same length ",
      "as 'theta_body'.",
      call. = FALSE
    )
  }

  weights / sum(weights)
}

# Resolve the receive-angle vectors for one orientation-averaged TMM solve.
#' @noRd
.tmm_average_scatter_angles <- function(theta_body,
                                        phi_body,
                                        theta_scatter,
                                        phi_scatter) {
  # Default to the exact monostatic receive direction when omitted ============
  n_angles <- length(theta_body)
  theta_scatter <- if (is.null(theta_scatter)) {
    pi - theta_body
  } else {
    rep_len(theta_scatter, n_angles)
  }
  phi_scatter <- if (is.null(phi_scatter)) {
    phi_body + pi
  } else {
    rep_len(phi_scatter, n_angles)
  }

  list(theta_scatter = theta_scatter, phi_scatter = phi_scatter)
}

#' Orientation-average scattering from a stored TMM object
#'
#' @description
#' Reuses the stored TMM blocks to average the differential scattering cross
#' section over a user-supplied set of incident orientations. By default the
#' receive direction is taken to be the exact monostatic direction for each
#' supplied orientation.
#'
#' @param object Scatterer-object previously evaluated with
#'   `target_strength(..., model = "TMM", store_t_matrix = TRUE)`.
#' @param theta_body Numeric vector of incident polar angles (radians).
#' @param weights Optional numeric vector of averaging weights. When omitted,
#'   equal weights are used.
#' @param phi_body Incident azimuth angle(s) (radians). Either scalar or the
#'   same length as `theta_body`. Defaults to `pi`.
#' @param theta_scatter Receive polar angle(s) (radians). Either scalar or the
#'   same length as `theta_body`. Defaults to the monostatic direction.
#' @param phi_scatter Receive azimuth angle(s) (radians). Either scalar or the
#'   same length as `theta_body`. Defaults to the monostatic direction.
#' @param distribution Optional orientation distribution created by
#'   \code{\link{tmm_orientation_distribution}}. When supplied, it overrides
#'   the direct `theta_body`, `weights`, and `phi_body` inputs.
#'
#' @return A data frame containing the frequency, the orientation-averaged
#'   differential backscattering cross section and the corresponding
#'   orientation-averaged target strength.
#'
#' @seealso \code{\link{tmm_scattering}}
#' @export
tmm_average_orientation <- function(object,
                                    theta_body = NULL,
                                    weights = NULL,
                                    phi_body = NULL,
                                    theta_scatter = NULL,
                                    phi_scatter = NULL,
                                    distribution = NULL) {
  # Recover the stored TMM state and any pre-built distribution ================
  model_params <- .tmm_require_stored_blocks(object)
  .tmm_warn_exploratory_cylinder_blocks(object, model_params)
  orientation <- .tmm_average_orientation_inputs(
    distribution = distribution,
    theta_body = theta_body,
    weights = weights,
    phi_body = phi_body,
    default_phi_body = model_params$body$phi_body %||% pi
  )
  theta_body <- orientation$theta_body
  phi_body <- orientation$phi_body
  weights <- orientation$weights

  # Resolve the receive direction, defaulting to exact monostatic geometry =====
  scatter_angles <- .tmm_average_scatter_angles(
    theta_body = theta_body,
    phi_body = phi_body,
    theta_scatter = theta_scatter,
    phi_scatter = phi_scatter
  )
  theta_scatter <- scatter_angles$theta_scatter
  phi_scatter <- scatter_angles$phi_scatter

  # Evaluate the stored scattering operator at each supplied orientation =======
  n_angles <- length(theta_body)
  sigma_mat <- vapply(
    seq_len(n_angles),
    function(i) {
      suppressWarnings(tmm_scattering(
        object = object,
        theta_body = theta_body[i],
        phi_body = phi_body[i],
        theta_scatter = theta_scatter[i],
        phi_scatter = phi_scatter[i]
      ))$sigma_scat
    },
    numeric(length(model_params$parameters$acoustics$frequency))
  )
  if (is.null(dim(sigma_mat))) {
    sigma_mat <- matrix(
      sigma_mat,
      nrow = length(model_params$parameters$acoustics$frequency),
      ncol = n_angles
    )
  }
  sigma_avg <- as.vector(sigma_mat %*% weights)

  # Return the orientation-averaged backscatter summaries ======================
  data.frame(
    frequency = model_params$parameters$acoustics$frequency,
    sigma_bs = sigma_avg,
    TS = db(sigma_avg)
  )
}
