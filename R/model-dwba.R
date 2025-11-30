#' Calculates the theoretical TS of a fluid-like scatterer at a given frequency
#' using the distorted Born wave approximation (DWBA) model.
#'
#' @param object FLS-class scatterer.
#' @references
#' Stanton, T.K., Chu, D., and Wiebe, P.H. 1998. Sound scattering by several
#' zooplankton groups. II. Scattering models. J. Acoust. Soc. Am., 103, 236-253.
#'
#' Demer, D.A., and Conti, S.G. 2003. Reconciling theoretical versus empirical
#' target strengths of krill: effects of phase variability on the distorted-wave
#' Born approximation. ICES J. Mar. Sci., 60: 429-434.
#' @export
DWBA <- function(object) {
  # Extract model parameters/inputs ============================================
  model <- acousticTS::extract(object, "model_parameters")$DWBA
  body <- acousticTS::extract(object, "body")
  theta <- body$theta
  r0 <- body$rpos[1:3, ]
  # Material properties calculation ============================================
  g <- body$g
  h <- body$h
  R <- 1 / (g * h * h) + 1 / g - 2
  # Calculate rotation matrix and update wavenumber matrix =====================
  rotation_matrix <- matrix(
    c(
      cos(theta),
      0.0,
      sin(theta)
    ),
    1
  )
  k_sw_rot <- model$parameters$acoustics$k_sw %*% rotation_matrix
  # Calculate Euclidean norms ==================================================
  k_sw_norm <- acousticTS::vecnorm(k_sw_rot)
  # Update position matrices  ==================================================
  rpos <- rbind(r0, a = model$body$radius)
  # Calculate position matrix lags  ============================================
  rpos_diff <- t(diff(t(rpos)))
  # Multiply wavenumber and body matrices ======================================
  rpos_diff_k <- t(
    vapply(
      X = seq_along(k_sw_norm),
      FUN = function(x) {
        colSums(rpos_diff[1:3, ] * k_sw_rot[x, ])
      },
      FUN.VALUE = numeric(ncol(rpos_diff))
    )
  )
  # Calculate Euclidean norms ==================================================
  rpos_diff_norm <- sqrt(colSums(rpos_diff[1:3, ] * rpos_diff[1:3, ]))
  # Estimate angles between body cylinders =====================================
  alpha <- acos(rpos_diff_k / (k_sw_norm %*% t(rpos_diff_norm)))
  beta <- abs(alpha - pi / 2)
  # Define integrand ===========================================================
  integrand <- function(s, x) {
    rint_mat <- s * rpos_diff + rpos[, seq_len(ncol(rpos_diff))]
    rint_k1_h_mat <- k_sw_rot[x, ] %*% rint_mat[1:3, ] / h
    bessel <- jc(1, 2 * (k_sw_norm[x] * rint_mat[4, ] / h *
                           cos(beta[x, ]))) / cos(beta[x, ])
    fb_a <- k_sw_norm[x] / 4 * R * rint_mat[4, ] *
      exp(2i * rint_k1_h_mat) * bessel * rpos_diff_norm
    sum(fb_a, na.rm = TRUE)
  }
  # Vectorize integrand function ===============================================
  integrand_vec <- Vectorize(integrand)
  # Calculate linear scatter response ==========================================
  f_bs <- rep(NA, length(k_sw_norm))
  f_bs <- vapply(
    seq_along(k_sw_norm),
    FUN = function(x) {
      # Real ===================================================================
      Ri <- stats::integrate(function(s) {
        Re(integrand_vec(s, x))
      }, 0, 1)$value
      # Real ===================================================================
      Ii <- stats::integrate(function(s) {
        Im(integrand_vec(s, x))
      }, 0, 1)$value
      # Return =================================================================
      sqrt(Ri^2 + Ii^2)
    },
    FUN.VALUE = numeric(1)
  )
  # Update scatterer object ====================================================
  methods::slot(object, "model")$DWBA <- data.frame(
    frequency = model$parameters$acoustics$frequency,
    ka = model$parameters$acoustics$k_sw * stats::median(model$body$radius),
    f_bs = f_bs,
    sigma_bs = abs(f_bs) * abs(f_bs),
    TS = 10 * log10(abs(f_bs) * abs(f_bs))
  )
  # Return object ==============================================================
  object
}
#' Calculates the theoretical TS of a fluid-like scatterer at a given frequency
#' using the distorted Born wave approximation (DWBA) model.
#'
#' @param object FLS-class scatterer.
#' @references
#' Stanton, T.K., Chu, D., and Wiebe, P.H. 1998. Sound scattering by several
#' zooplankton groups. II. Scattering models. J. Acoust. Soc. Am., 103, 236-253.
#'
#' Demer, D.A., and Conti, S.G. 2003. Reconciling theoretical versus empirical
#' target strengths of krill: effects of phase variability on the distorted-wave
#' Born approximation. ICES J. Mar. Sci., 60: 429-434.
#' @export
DWBA_curved <- function(object) {
  # Extract model parameters/inputs ============================================
  model <- acousticTS::extract(object, "model_parameters")$DWBA_curved
  body <- acousticTS::extract(object, "body")
  # Parse position matrices ====================================================
  r0 <- model$body$rpos[1:3, ]
  theta <- model$body$theta
  # Material properties calculation ============================================
  g <- body$g
  h <- body$h
  R <- 1 / (g * h * h) + 1 / g - 2
  # Calculate rotation matrix and update wavenumber matrix =====================
  rotation_matrix <- matrix(
    c(
      cos(theta),
      0.0,
      sin(theta)
    ),
    1
  )
  k_sw_rot <- model$parameters$acoustics$k_sw %*% rotation_matrix
  # Calculate Euclidean norms ==================================================
  k_sw_norm <- acousticTS::vecnorm(k_sw_rot)
  # Update position matrices  ==================================================
  rpos <- rbind(r0, a = model$body$radius)
  # Calculate position matrix lags  ============================================
  rpos_diff <- t(diff(t(rpos)))
  # Multiply wavenumber and body matrices ======================================
  rpos_diff_k <- t(vapply(seq_along(k_sw_norm),
                          FUN = function(x) {
                            colSums(rpos_diff[1:3, ] * k_sw_rot[x, ])
                          },
                          FUN.VALUE = numeric(ncol(rpos_diff))
  ))
  # Calculate Euclidean norms ==================================================
  rpos_diff_norm <- sqrt(colSums(rpos_diff[1:3, ] * rpos_diff[1:3, ]))
  # Estimate angles between body cylinders =====================================
  alpha <- acos(rpos_diff_k / (k_sw_norm %*% t(rpos_diff_norm)))
  beta <- abs(alpha - pi / 2)
  # Define integrand ===========================================================
  integrand <- function(s, x) {
    rint_mat <- s * rpos_diff + rpos[, seq_len(ncol(rpos_diff))]
    rint_k1_h_mat <- k_sw_rot[x, ] %*% rint_mat[1:3, ] / h
    bessel <- jc(1, 2 * (k_sw_norm[x] * rint_mat[4, ] / h *
                           cos(beta[x, ]))) / cos(beta[x, ])
    fb_a <- k_sw_norm[x] / 4 * R * rint_mat[4, ] *
      exp(2i * rint_k1_h_mat) * bessel * rpos_diff_norm
    sum(fb_a, na.rm = TRUE)
  }
  # Vectorize integrand function ===============================================
  integrand_vec <- Vectorize(integrand)
  # Calculate linear scatter response ==========================================
  f_bs <- rep(NA, length(k_sw_norm))
  f_bs <- vapply(seq_along(k_sw_norm),
                 FUN = function(x) {
                   # Real ===============================================
                   Ri <- stats::integrate(function(s) {
                     Re(integrand_vec(s, x))
                   }, 0, 1)$value
                   # Real ===============================================
                   Ii <- stats::integrate(function(s) {
                     Im(integrand_vec(s, x))
                   }, 0, 1)$value
                   # Return =============================================
                   sqrt(Ri^2 + Ii^2)
                 },
                 FUN.VALUE = numeric(1)
  )
  # Update scatterer object ====================================================
  methods::slot(object, "model")$DWBA_curved <- data.frame(
    frequency = model$parameters$acoustics$frequency,
    ka = model$parameters$acoustics$k_sw * stats::median(model$body$radius),
    f_bs = f_bs,
    sigma_bs = abs(f_bs) * abs(f_bs),
    TS = 10 * log10(abs(f_bs) * abs(f_bs))
  )
  # Return object ==============================================================
  object
}
