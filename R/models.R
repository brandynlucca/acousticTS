#' Model registry comprising target strength models currently supported by the
#' acousticTS R-package.
#' @keywords internal
#' @noRd
.get_models <- function() {
  list(
    DWBA = DWBA,
    DWBA_curved = DWBA_curved,
    SDWBA = SDWBA,
    SDWBA_curved = SDWBA_curved,
    calibration = calibration,
    SOEMS = calibration,
    ESSMS = ESSMS,
    SPHMS = SPHMS,
    KRM = KRM,
    HPA = HPA,
    SPHMS = SPHMS,
    ESSMS = ESSMS,
    PSMS = PSMS,
    FCMS = FCMS
  )
}


#' ################################################################################
#' # Primary scattering models for fluid-like scatterers (FLS)
#' ################################################################################
#' #' Calculates the theoretical TS of a target using the deformed cylinder
#' #' model (DCM).
#' #'
#' #' @param object FLS-class scatterer.
#' #' @usage
#' #' DCM( object )
#' #' @return
#' #' Target strength (TS, dB re. 1 \eqn{m^2}) of a FLS-class object.
#' #' #' @references
#' #' Stanton, T.K., Chu, D., and Wiebe, P.H. 1998. Sound scattering by several
#' #' zooplankton groups. II. Scattering models. Journal of the Acoustical Society
#' #' of America, 103(1), 236-253.
#' #' @export
#' DCM <- function(object) {
#'   # Extract model parameters/inputs ============================================
#'   model <- acousticTS::extract(object, "model_parameters")$DCM
#'   body <- acousticTS::extract(object, "body")
#'   body$density <- model$medium$density * body$g
#'   body$sound_speed <- model$medium$sound_speed * body$h
#'   # Multiply acoustic wavenumber by body radius ================================
#'   k1a <- model$parameters$acoustics$k_sw * model$body$radius
#'   k2a <- model$parameters$acoustics$k_b * model$ body$radius
#'   # Calculate Reflection Coefficient ===========================================
#'   R12 <- acousticTS::reflection_coefficient(model$medium, body)
#'   # Calculate Transmission Coefficients ========================================
#'   T12 <- 1 - R12
#'   T21 <- 1 + R12
#'   # Calculate phase shift term for the ray model ===============================
#'   mu <- (-pi / 2 * k1a) / (k1a + 0.4)
#'   # Calculate Interference Term between echoes from front/back interfaces ======
#'   I0 <- 1 - T12 * T21 * exp(4i * k2a) * exp(1i * mu * k1a)
#'   # Calculate half-width of directivity pattern ================================
#'   w12 <- sqrt(model$body$radius_curvature * model$body$radius)
#'   # Apply shift to orientation required for this model =========================
#'   theta_shift <- body$theta - pi / 2
#'   # Calculate Taylor expansion of directivity function =========================
#'   D_theta <- exp(-model$body$alpha_B *
#'                          (2 * theta_shift * model$body$radius_curvature /
#'                             model$body$length)^2)
#'   # Calculate linear backscatter function, f_bs ================================
#'   f_bs <- 0.5 * w12 * R12 * exp(-2i * k1a) * I0 * D_theta
#'   # Update scatterer object ====================================================
#'   methods::slot(object, "model")$DCM <- data.frame(
#'     frequency = model$parameters$acoustics$frequency,
#'     ka = k1a,
#'     f_bs = f_bs,
#'     sigma_bs = abs(f_bs)^2,
#'     TS = 10 * log10(abs(f_bs)^2)
#'   )
#'   # Return object ==============================================================
#'   object
#' }
#' #' Calculates the theoretical TS of a fluid-like scatterer at a given frequency
#' #' using the distorted Born wave approximation (DWBA) model.
#' #'
#' #' @param object FLS-class scatterer.
#' #' @references
#' #' Stanton, T.K., Chu, D., and Wiebe, P.H. 1998. Sound scattering by several
#' #' zooplankton groups. II. Scattering models. J. Acoust. Soc. Am., 103, 236-253.
#' #'
#' #' Demer, D.A., and Conti, S.G. 2003. Reconciling theoretical versus empirical
#' #' target strengths of krill: effects of phase variability on the distorted-wave
#' #' Born approximation. ICES J. Mar. Sci., 60: 429-434.
#' #' @export
#' DWBA <- function(object) {
#'   # Extract model parameters/inputs ============================================
#'   model <- acousticTS::extract(object, "model_parameters")$DWBA
#'   body <- acousticTS::extract(object, "body")
#'   theta <- body$theta
#'   r0 <- body$rpos[1:3, ]
#'   # Material properties calculation ============================================
#'   g <- body$g
#'   h <- body$h
#'   R <- 1 / (g * h * h) + 1 / g - 2
#'   # Calculate rotation matrix and update wavenumber matrix =====================
#'   rotation_matrix <- matrix(
#'     c(
#'       cos(theta),
#'       0.0,
#'       sin(theta)
#'     ),
#'     1
#'   )
#'   k_sw_rot <- model$parameters$acoustics$k_sw %*% rotation_matrix
#'   # Calculate Euclidean norms ==================================================
#'   k_sw_norm <- acousticTS::vecnorm(k_sw_rot)
#'   # Update position matrices  ==================================================
#'   rpos <- rbind(r0, a = model$body$radius)
#'   # Calculate position matrix lags  ============================================
#'   rpos_diff <- t(diff(t(rpos)))
#'   # Multiply wavenumber and body matrices ======================================
#'   rpos_diff_k <- t(
#'     vapply(
#'       X = seq_along(k_sw_norm),
#'       FUN = function(x) {
#'         colSums(rpos_diff[1:3, ] * k_sw_rot[x, ])
#'       },
#'       FUN.VALUE = numeric(ncol(rpos_diff))
#'     )
#'   )
#'   # Calculate Euclidean norms ==================================================
#'   rpos_diff_norm <- sqrt(colSums(rpos_diff[1:3, ] * rpos_diff[1:3, ]))
#'   # Estimate angles between body cylinders =====================================
#'   alpha <- acos(rpos_diff_k / (k_sw_norm %*% t(rpos_diff_norm)))
#'   beta <- abs(alpha - pi / 2)
#'   # Define integrand ===========================================================
#'   integrand <- function(s, x) {
#'     rint_mat <- s * rpos_diff + rpos[, seq_len(ncol(rpos_diff))]
#'     rint_k1_h_mat <- k_sw_rot[x, ] %*% rint_mat[1:3, ] / h
#'     bessel <- jc(1, 2 * (k_sw_norm[x] * rint_mat[4, ] / h *
#'                            cos(beta[x, ]))) / cos(beta[x, ])
#'     fb_a <- k_sw_norm[x] / 4 * R * rint_mat[4, ] *
#'       exp(2i * rint_k1_h_mat) * bessel * rpos_diff_norm
#'     sum(fb_a, na.rm = TRUE)
#'   }
#'   # Vectorize integrand function ===============================================
#'   integrand_vec <- Vectorize(integrand)
#'   # Calculate linear scatter response ==========================================
#'   f_bs <- rep(NA, length(k_sw_norm))
#'   f_bs <- vapply(
#'     seq_along(k_sw_norm),
#'     FUN = function(x) {
#'       # Real ===================================================================
#'       Ri <- stats::integrate(function(s) {
#'         Re(integrand_vec(s, x))
#'       }, 0, 1)$value
#'       # Real ===================================================================
#'       Ii <- stats::integrate(function(s) {
#'         Im(integrand_vec(s, x))
#'       }, 0, 1)$value
#'       # Return =================================================================
#'       sqrt(Ri^2 + Ii^2)
#'     },
#'     FUN.VALUE = numeric(1)
#'   )
#'   # Update scatterer object ====================================================
#'   methods::slot(object, "model")$DWBA <- data.frame(
#'     frequency = model$parameters$acoustics$frequency,
#'     ka = model$parameters$acoustics$k_sw * stats::median(model$body$radius),
#'     f_bs = f_bs,
#'     sigma_bs = abs(f_bs) * abs(f_bs),
#'     TS = 10 * log10(abs(f_bs) * abs(f_bs))
#'   )
#'   # Return object ==============================================================
#'   object
#' }
#' #' Calculates the theoretical TS of a fluid-like scatterer at a given frequency
#' #' using the distorted Born wave approximation (DWBA) model.
#' #'
#' #' @param object FLS-class scatterer.
#' #' @references
#' #' Stanton, T.K., Chu, D., and Wiebe, P.H. 1998. Sound scattering by several
#' #' zooplankton groups. II. Scattering models. J. Acoust. Soc. Am., 103, 236-253.
#' #'
#' #' Demer, D.A., and Conti, S.G. 2003. Reconciling theoretical versus empirical
#' #' target strengths of krill: effects of phase variability on the distorted-wave
#' #' Born approximation. ICES J. Mar. Sci., 60: 429-434.
#' #' @export
#' DWBA_curved <- function(object) {
#'   # Extract model parameters/inputs ============================================
#'   model <- acousticTS::extract(object, "model_parameters")$DWBA_curved
#'   body <- acousticTS::extract(object, "body")
#'   # Parse position matrices ====================================================
#'   r0 <- model$body$rpos[1:3, ]
#'   theta <- model$body$theta
#'   # Material properties calculation ============================================
#'   g <- body$g
#'   h <- body$h
#'   R <- 1 / (g * h * h) + 1 / g - 2
#'   # Calculate rotation matrix and update wavenumber matrix =====================
#'   rotation_matrix <- matrix(
#'     c(
#'       cos(theta),
#'       0.0,
#'       sin(theta)
#'     ),
#'     1
#'   )
#'   k_sw_rot <- model$parameters$acoustics$k_sw %*% rotation_matrix
#'   # Calculate Euclidean norms ==================================================
#'   k_sw_norm <- acousticTS::vecnorm(k_sw_rot)
#'   # Update position matrices  ==================================================
#'   rpos <- rbind(r0, a = model$body$radius)
#'   # Calculate position matrix lags  ============================================
#'   rpos_diff <- t(diff(t(rpos)))
#'   # Multiply wavenumber and body matrices ======================================
#'   rpos_diff_k <- t(vapply(seq_along(k_sw_norm),
#'     FUN = function(x) {
#'       colSums(rpos_diff[1:3, ] * k_sw_rot[x, ])
#'     },
#'     FUN.VALUE = numeric(ncol(rpos_diff))
#'   ))
#'   # Calculate Euclidean norms ==================================================
#'   rpos_diff_norm <- sqrt(colSums(rpos_diff[1:3, ] * rpos_diff[1:3, ]))
#'   # Estimate angles between body cylinders =====================================
#'   alpha <- acos(rpos_diff_k / (k_sw_norm %*% t(rpos_diff_norm)))
#'   beta <- abs(alpha - pi / 2)
#'   # Define integrand ===========================================================
#'   integrand <- function(s, x) {
#'     rint_mat <- s * rpos_diff + rpos[, seq_len(ncol(rpos_diff))]
#'     rint_k1_h_mat <- k_sw_rot[x, ] %*% rint_mat[1:3, ] / h
#'     bessel <- jc(1, 2 * (k_sw_norm[x] * rint_mat[4, ] / h *
#'                            cos(beta[x, ]))) / cos(beta[x, ])
#'     fb_a <- k_sw_norm[x] / 4 * R * rint_mat[4, ] *
#'       exp(2i * rint_k1_h_mat) * bessel * rpos_diff_norm
#'     sum(fb_a, na.rm = TRUE)
#'   }
#'   # Vectorize integrand function ===============================================
#'   integrand_vec <- Vectorize(integrand)
#'   # Calculate linear scatter response ==========================================
#'   f_bs <- rep(NA, length(k_sw_norm))
#'   f_bs <- vapply(seq_along(k_sw_norm),
#'     FUN = function(x) {
#'       # Real ===============================================
#'       Ri <- stats::integrate(function(s) {
#'         Re(integrand_vec(s, x))
#'       }, 0, 1)$value
#'       # Real ===============================================
#'       Ii <- stats::integrate(function(s) {
#'         Im(integrand_vec(s, x))
#'       }, 0, 1)$value
#'       # Return =============================================
#'       sqrt(Ri^2 + Ii^2)
#'     },
#'     FUN.VALUE = numeric(1)
#'   )
#'   # Update scatterer object ====================================================
#'   methods::slot(object, "model")$DWBA_curved <- data.frame(
#'     frequency = model$parameters$acoustics$frequency,
#'     ka = model$parameters$acoustics$k_sw * stats::median(model$body$radius),
#'     f_bs = f_bs,
#'     sigma_bs = abs(f_bs) * abs(f_bs),
#'     TS = 10 * log10(abs(f_bs) * abs(f_bs))
#'   )
#'   # Return object ==============================================================
#'   object
#' }
#' #' Calculates the theoretical TS of a fluid-like scatterer at a given frequency
#' #' using the stochastic distorted Born wave approximation (DWBA) model.
#' #'
#' #' @param object FLS-class scatterer.
#' #' @references
#' #' Demer, D.A., and Conti, S.G. 2003. Reconciling theoretical versus empirical
#' #' target strengths of krill: effects of phase variability on the distorted-wave
#' #' Born approximation. ICES J. Mar. Sci., 60: 429-434.
#' #' @export
#' SDWBA <- function(object) {
#'   # Extract model parameters/inputs ============================================
#'   model <- acousticTS::extract(object, "model_parameters")$SDWBA
#'   body <- acousticTS::extract(object, "body")
#'   theta <- body$theta
#'   # Material properties calculation ============================================
#'   g <- body$g
#'   h <- body$h
#'   R <- 1 / (g * h * h) + 1 / g - 2
#'   # Calculate rotation matrix and update wavenumber matrix =====================
#'   rotation_matrix <- matrix(
#'     c(
#'       cos(theta),
#'       0.0,
#'       sin(theta)
#'     ),
#'     1
#'   )
#'   SDWBA_resampled <- function(i) {
#'     # Parse phase-specific parameters ==========================================
#'     sub_params <- model$parameters[[i]]
#'     # Calculate rotation matrix and update wavenumber matrix ===================
#'     k_sw_rot <- sub_params$acoustics$k_sw %*% rotation_matrix
#'     # Calculate Euclidean norms ================================================
#'     k_sw_norm <- vecnorm(k_sw_rot)
#'     # Update position matrices  ================================================
#'     r0 <- rbind(
#'       sub_params$body_params$rpos[1:3, ],
#'       sub_params$body_params$radius
#'     )
#'     # Calculate position matrix lags  ==========================================
#'     r0_diff <- t(diff(t(r0)))
#'     # Multiply wavenumber and body matrices ====================================
#'     r0_diff_k <- t(vapply(seq_along(k_sw_norm),
#'       FUN = function(x) {
#'         colSums(r0_diff[1:3, ] *
#'                   k_sw_rot[x, ])
#'       },
#'       FUN.VALUE = numeric(ncol(r0_diff))
#'     ))
#'     # Calculate Euclidean norms ================================================
#'     r0_diff_norm <- sqrt(colSums(r0_diff[1:3, ] * r0_diff[1:3, ]))
#'     # Estimate angles between body cylinders ===================================
#'     alpha <- acos(r0_diff_k / (k_sw_norm %*% t(r0_diff_norm)))
#'     beta <- abs(alpha - pi / 2)
#'     # Call in metrics ==========================================================
#'     phase_sd <- sub_params$meta_params$phase_sd
#'     # r0_diff_h <- r0_diff / h
#'     # r0_h <- r0 / h
#'     # Define integrand =========================================================
#'     integrand <- function(s, x, y) {
#'       # integrand <- function( s , x ) {
#'       # rint_mat <- s * r0_diff_h[ , y ] + r0_h[ , y ]
#'       rint_mat <- s * r0_diff[, y] + r0[, y]
#'       # rint_k1_h_mat <- k_sw_rot[ x , ] %*% rint_mat[ 1 : 3 ]
#'       rint_k1_h_mat <- k_sw_rot[x, ] %*% rint_mat[1:3] / h
#'       bessel <- jc(1, 2 * (k_sw_norm[x] * rint_mat[4] / h *
#'                              cos(beta[x, y]))) / cos(beta[x, y])
#'       # bessel <- jc( 1 , 2 * ( k_sw_norm[ x ] * rint_mat[ 4 ] *
#'       #                           cos( beta[ x , y ] ) ) ) / cos( beta[ x , y] )
#'       # fb_a <- k_sw_norm[ x ] / 4 * R * rint_mat[ 4 ] * h *
#'       #   exp( 2i * rint_k1_h_mat ) * bessel[y] * r0_diff_norm[ y ]
#'       fb_a <- k_sw_norm[x] / 4 * R * rint_mat[4] *
#'         exp(2i * rint_k1_h_mat) * bessel * r0_diff_norm[y]
#'       sum(fb_a, na.rm = TRUE)
#'     }
#'     # Vectorize integrand function =============================================
#'     integrand_vectorized <- Vectorize(integrand)
#'     stochastic_TS <- function(n_k, n_segments, n_iterations) {
#'       phase_cyl <- vapply(seq_len(n_k),
#'         FUN = function(x) {
#'           cyl_phase <- array(
#'             vapply(seq_len(n_segments),
#'               FUN = function(y) {
#'                 phase_integrate(
#'                   x, y,
#'                   n_iterations,
#'                   integrand_vectorized,
#'                   phase_sd
#'                 )
#'               },
#'               FUN.VALUE = complex(n_iterations)
#'             ),
#'             dim = c(n_iterations, n_segments)
#'           )
#'           cyl_sum_phase <- rowSums(cyl_phase, na.rm = TRUE)
#'           cyl_sum_phase
#'         },
#'         FUN.VALUE = complex(n_iterations)
#'       )
#'       phase_cyl <- array(phase_cyl, dim = c(n_iterations, n_k))
#'       data.frame(
#'         f_bs = colMeans(phase_cyl),
#'         sigma_bs = colMeans(.sigma_bs(phase_cyl)),
#'         TS_mean = db(colMeans(.sigma_bs(phase_cyl))),
#'         TS_sd = db(apply(.sigma_bs(phase_cyl), 2, stats::sd))
#'       )
#'     }
#'     # Calculate linear scatter response ========================================
#'     backscatter_df <- stochastic_TS(
#'       n_k = length(k_sw_norm),
#'       n_segments = sub_params$n_segments - 1,
#'       n_iterations = sub_params$meta_params$n_iterations
#'     )
#'     backscatter_df
#'   }
#'   # Generate results dataframe by collating resampled results ==================
#'   results <- do.call(
#'     "rbind",
#'     lapply(seq_along(model$parameters),
#'       FUN = function(i) SDWBA_resampled(i)
#'     )
#'   )
#'   # Update scatterer object ====================================================
#'   methods::slot(object, "model")$SDWBA$f_bs <- results$f_bs
#'   methods::slot(object, "model")$SDWBA$sigma_bs <- results$sigma_bs
#'   methods::slot(object, "model")$SDWBA$TS <- results$TS_mean
#'   methods::slot(object, "model")$SDWBA$TS_sd <- results$TS_sd
#'   # Return object ==============================================================
#'   object
#' }
#' #' Calculates the theoretical TS of a fluid-like scatterer at a given frequency
#' #' using the stochastic distorted Born wave approximation (DWBA) model using
#' #' Eq. (6) from Stanton et al. (1998).
#' #' @param object FLS-class scatterer.
#' #' @references
#' #' Stanton, T.K., Chu, D., and Wiebe, P.H. 1998. Sound scattering by several
#' #' zooplankton groups. II. Scattering models. J. Acoust. Soc. Am., 103, 236-253.
#' #' @rdname SDWBA_curved
#' #' @export
#' SDWBA_curved <- function(object) {
#'   # Extract model parameters/inputs ============================================
#'   model <- acousticTS::extract(object, "model_parameters")$SDWBA_curved
#'   body <- acousticTS::extract(object, "body")
#'   # Parse position matrices ====================================================
#'   theta <- model$body$theta
#'   # Material properties calculation ============================================
#'   g <- body$g
#'   h <- body$h
#'   R <- 1 / (g * h * h) + 1 / g - 2
#'   # Calculate rotation matrix and update wavenumber matrix =====================
#'   rotation_matrix <- matrix(
#'     c(
#'       cos(theta),
#'       0.0,
#'       sin(theta)
#'     ),
#'     1
#'   )
#'   SDWBA_resampled_c <- function(i) {
#'     # Parse phase-specific parameters ==========================================
#'     sub_params <- model$parameters[[i]]
#'     # Calculate rotation matrix and update wavenumber matrix ===================
#'     k_sw_rot <- sub_params$acoustics$k_sw %*% rotation_matrix
#'     # Calculate Euclidean norms ================================================
#'     k_sw_norm <- vecnorm(k_sw_rot)
#'     # Update position matrices  ================================================
#'     r0 <- rbind(
#'       sub_params$body_params$rpos[1:3, ],
#'       sub_params$body_params$radius
#'     )
#'     # Calculate position matrix lags  ==========================================
#'     r0_diff <- t(diff(t(r0)))
#'     # Multiply wavenumber and body matrices ====================================
#'     r0_diff_k <- t(vapply(seq_along(k_sw_norm),
#'       FUN = function(x) {
#'         colSums(r0_diff[1:3, ] *
#'                   k_sw_rot[x, ])
#'       },
#'       FUN.VALUE = numeric(ncol(r0_diff))
#'     ))
#'     # Calculate Euclidean norms ================================================
#'     r0_diff_norm <- sqrt(colSums(r0_diff[1:3, ] * r0_diff[1:3, ]))
#'     # Estimate angles between body cylinders ===================================
#'     alpha <- acos(r0_diff_k / (k_sw_norm %*% t(r0_diff_norm)))
#'     beta <- abs(alpha - pi / 2)
#'     # Call in metrics ==========================================================
#'     phase_sd <- sub_params$meta_params$phase_sd
#'     r0_diff_h <- r0_diff / h
#'     r0_h <- r0 / h
#'     # Define integrand =========================================================
#'     integrand_c <- function(s, x, y) {
#'       rint_mat <- s * r0_diff_h[, y] + r0_h[, y]
#'       rint_k1_h_mat <- k_sw_rot[x, ] %*% rint_mat[1:3]
#'       bessel <- jc(1, 2 * (k_sw_norm[x] * rint_mat[4] *
#'                              cos(beta[x, y]))) / cos(beta[x, y])
#'       fb_a <- k_sw_norm[x] / 4 * R * rint_mat[4] * h *
#'         exp(2i * rint_k1_h_mat) * bessel * r0_diff_norm[y]
#'       sum(fb_a, na.rm = TRUE)
#'     }
#'     # Vectorize integrand function =============================================
#'     integrand_vectorized_c <- Vectorize(integrand_c)
#'     stochastic_TS_c <- function(n_k, n_segments, n_iterations) {
#'       phase_cyl <- vapply(seq_len(n_k),
#'         FUN = function(x) {
#'           cyl_phase <- array(
#'             vapply(seq_len(n_segments),
#'               FUN = function(y) {
#'                 phase_integrate(
#'                   x, y,
#'                   n_iterations,
#'                   integrand_vectorized_c,
#'                   phase_sd
#'                 )
#'               },
#'               FUN.VALUE = complex(n_iterations)
#'             ),
#'             dim = c(n_iterations, n_segments)
#'           )
#'           cyl_sum_phase <- rowSums(cyl_phase, na.rm = TRUE)
#'           cyl_sum_phase
#'         },
#'         FUN.VALUE = complex(n_iterations)
#'       )
#'       phase_cyl <- array(phase_cyl, dim = c(n_iterations, n_k))
#'       data.frame(
#'         f_bs = colMeans(phase_cyl),
#'         sigma_bs = colMeans(.sigma_bs(phase_cyl)),
#'         TS_mean = db(colMeans(.sigma_bs(phase_cyl))),
#'         TS_sd = db(apply(.sigma_bs(phase_cyl), 2, stats::sd))
#'       )
#'     }
#'     # Calculate linear scatter response ========================================
#'     backscatter_df <- stochastic_TS_c(
#'       n_k = length(k_sw_norm),
#'       n_segments = sub_params$n_segments - 1,
#'       n_iterations = sub_params$meta_params$n_iterations
#'     )
#'     backscatter_df
#'   }
#'   # Generate results dataframe by collating resampled results ==================
#'   results <- do.call(
#'     "rbind",
#'     lapply(seq_along(model$parameters),
#'       FUN = function(i) SDWBA_resampled_c(i)
#'     )
#'   )
#'   # Update scatterer object ====================================================
#'   methods::slot(object, "model")$SDWBA_curved$f_bs <- results$f_bs
#'   methods::slot(object, "model")$SDWBA_curved$sigma_bs <- results$sigma_bs
#'   methods::slot(object, "model")$SDWBA_curved$TS <- results$TS_mean
#'   methods::slot(object, "model")$SDWBA_curved$TS_sd <- results$TS_sd
#'   # Return object ==============================================================
#'   object
#' }
#' ################################################################################
#' # Primary scattering model for an elastic solid sphere (CAL)
#' ################################################################################
#' #' Calculates theoretical TS of a solid sphere of a certain material at a given
#' #' frequency.
#' #'
#' #' @description
#' #' This function is a wrapper around TS_calculate(...) that parametrizes the
#' #' remainder of the model, while also doing simple calculations that do not
#' #' need to be looped. This function provides a TS estimate at a given frequency.
#' #' @param object CAL-class object.
#' #' @usage
#' #' calibration(object)
#' #' @return The theoretical acoustic target strength (TS, dB re. 1 \eqn{m^2}) of
#' #' a solid sphere at a given frequency.
#' #' @references
#' #' MacLennan D. N. (1981). The theory of solid spheres as sonar calibration
#' #' targets. Scottish Fisheries Research No. 22, Department of Agriculture and
#' #' Fisheries for Scotland.
#' #' @export
#' calibration <- function(object) {
#'   # Extract model parameters/inputs ============================================
#'   model <- acousticTS::extract(
#'     object,
#'     "model_parameters"
#'   )$calibration
#'   ### Now we solve / calculate equations =======================================
#'   # Equations 6a -- weight wavenumber by radius ================================
#'   ka_sw <- model$parameters$acoustics$k_sw * model$body$radius
#'   ka_l <- model$parameters$acoustics$k_l * model$body$radius
#'   ka_t <- model$parameters$acoustics$k_t * model$body$radius
#'   # Set limit for iterations ===================================================
#'   m_limit <- round(max(ka_sw)) + 10
#'   # Create modal series number vector ==========================================
#'   m <- 0:m_limit
#'   # Convert these vectors into matrices ========================================
#'   ka_sw_m <- modal_matrix(ka_sw, m_limit)
#'   ka_l_m <- modal_matrix(ka_l, m_limit)
#'   ka_t_m <- modal_matrix(ka_t, m_limit)
#'   # Calculate Legendre polynomial ==============================================
#'   Pl <- Pn(m, cos(model$body$theta))
#'   # Calculate spherical Bessel functions of first kind =========================
#'   js_mat <- js(m, ka_sw_m)
#'   js_mat_l <- js(m, ka_l_m)
#'   js_mat_t <- js(m, ka_t_m)
#'   # Calculate spherical Bessel functions of second kind ========================
#'   ys_mat <- ys(m, ka_sw_m)
#'   # Calculate first derivative of spheric Bessel functions of first kind =======
#'   jsd_mat <- jsd(m, ka_sw_m)
#'   jsd_mat_l <- jsd(m, ka_l_m)
#'   jsd_mat_t <- jsd(m, ka_t_m)
#'   # Calculate first derivative of spheric Bessel functions of second kind ======
#'   ysd_mat <- ysd(m, ka_sw_m)
#'   # Calculate density contrast =================================================
#'   g <- model$body$density / model$medium$density
#'   # Tangent functions ==========================================================
#'   tan_sw <- -ka_sw_m * jsd_mat / js_mat
#'   tan_l <- -ka_l_m * jsd_mat_l / js_mat_l
#'   tan_t <- -ka_t_m * jsd_mat_t / js_mat_t
#'   tan_beta <- -ka_sw_m * ysd_mat / ys_mat
#'   tan_diff <- -js_mat / ys_mat
#'   # Difference terms ===========================================================
#'   along_m <- (m * m + m)
#'   tan_l_add <- tan_l + 1
#'   tan_t_div <- along_m - 1 - ka_t_m * ka_t_m / 2 + tan_t
#'   numerator <- (tan_l / tan_l_add) - (along_m / tan_t_div)
#'   denominator1 <- (along_m - ka_t_m * ka_t_m / 2 + 2 * tan_l) / tan_l_add
#'   denominator2 <- along_m * (tan_t + 1) / tan_t_div
#'   denominator <- denominator1 - denominator2
#'   ratio <- -0.5 * (ka_t_m * ka_t_m) * numerator / denominator
#'   # Additional trig functions ==================================================
#'   phi <- -ratio / g
#'   eta_tan <- tan_diff * (phi + tan_sw) / (phi + tan_beta)
#'   cos_eta <- 1 / sqrt(1 + eta_tan * eta_tan)
#'   sin_eta <- eta_tan * cos_eta
#'   # Fill in rest of Hickling (1962) equation ===================================
#'   f_j <- colSums(
#'     (2 * m + 1) * Pl[m + 1] * (sin_eta * (1i * cos_eta - sin_eta))
#'   )
#'   # Calculate linear backscatter coefficient ===================================
#'   f_bs <- abs(-2i * f_j / ka_sw) * model$body$radius / 2
#'   sigma_bs <- f_bs * f_bs
#'   TS <- 10 * log10(sigma_bs)
#'   # Add results to scatterer object ============================================
#'   methods::slot(
#'     object,
#'     "model"
#'   )$calibration <- data.frame(
#'     frequency = model$parameters$acoustics$frequency,
#'     ka = ka_sw,
#'     f_bs = f_bs,
#'     sigma_bs = sigma_bs,
#'     TS = TS
#'   )
#'   # Return object ==============================================================
#'   object
#' }
#' ################################################################################
#' ################################################################################
#' # GAS-FILLED SCATTERERS
#' ################################################################################
#' ################################################################################
#' # Anderson (1950) gas-filled fluid-sphere model
#' ################################################################################
#' #' Calculates the theoretical TS of a fluid sphere using an exact modal series
#' #' solution proposed by Andersen (1950).
#' #' @param object GAS- or SBF-class object.
#' #' @details
#' #' Calculates the theoretical TS of a fluid sphere at a given frequency using
#' #' an exact modal series solution.
#' #' @return
#' #' Target strength (TS, dB re: 1 m^2)
#' #' @references
#' #' Anderson, V.C. (1950). Sound scattering from a fluid sphere. Journal of the
#' #' Acoustical Society of America, 22, 426-431.
#' #' @export
#' MSS_anderson <- function(object) {
#'   # Extract model parameters/inputs ============================================
#'   model <- extract(object, "model_parameters")$MSS_anderson
#'   # Multiple acoustic wavenumber by radius =====================================
#'   ## Medium ====================================================================
#'   k1a <- model$parameters$acoustics$k_sw * model$body$radius
#'   ## Fluid within sphere =======================================================
#'   k2a <- model$parameters$acoustics$k_f * model$body$radius
#'   # Set limit for iterations ===================================================
#'   m_limit <- model$parameters$ka_limit
#'   m <- 0:m_limit
#'   # Convert ka vectors to matrices =============================================
#'   ka1_m <- modal_matrix(k1a, m_limit)
#'   ka2_m <- modal_matrix(k2a, m_limit)
#'   # Calculate modal series coefficient, b_m (or C_m) ===========================
#'   ## Material properties term ==================================================
#'   gh <- model$body$g * model$body$h
#'   # Numerator term =============================================================
#'   N1 <- (jsd(m, ka2_m) * ys(m, ka1_m)) /
#'     (js(m, ka2_m) * jsd(m, ka1_m))
#'   N2 <- (gh * ysd(m, ka1_m) / jsd(m, ka1_m))
#'   CN <- N1 - N2
#'   # Denominator term ===========================================================
#'   D1 <- (jsd(m, ka2_m) * js(m, ka1_m)) /
#'     (js(m, ka2_m) * jsd(m, ka1_m))
#'   D2 <- gh
#'   CD <- D1 - D2
#'   # Finalize modal series coefficient ==========================================
#'   C_m <- CN / CD
#'   b_m <- -1 / (1 + 1i * C_m)
#'   # Calculate linear scatter response ==========================================
#'   ## Sum across columns to complete modal series summation over frequency range=
#'   f_j <- colSums((2 * m + 1) * (-1)^m * b_m)
#'   f_sphere <- -1i / model$parameters$acoustics$k_sw * f_j
#'   # Calculate linear scattering coefficient, sigma_bs ==========================
#'   sigma_bs <- abs(f_sphere) * abs(f_sphere)
#'   # Calculate TS ===============================================================
#'   TS <- 10 * log10(sigma_bs)
#'   # Add results to scatterer object ============================================
#'   methods::slot(
#'     object,
#'     "model"
#'   )$MSS_anderson <- data.frame(
#'     frequency = model$parameters$acoustics$frequency,
#'     ka = k1a,
#'     f_bs = f_sphere,
#'     sigma_bs = sigma_bs,
#'     TS = TS
#'   )
#'   # Return object ==============================================================
#'   object
#' }
#'
#' #' Calculate Bessel function cache for the Goodman and Stern (1962) model
#' #' @param ka_matrix_m Modal ka matrix
#' #' @param m Modal vector
#' #' @return Cached Bessel function values
#' #' @keywords internal
#' #' @noRd
#' .calculate_bessel_cache <- function(ka_matrix_m, m) {
#'   # # Define all Bessel functions to apply =======================================
#'   # bessel_functions <- list(
#'   #   js = js, jsd = jsd, jsdd = jsdd,
#'   #   ys = ys, ysd = ysd, ysdd = ysdd,
#'   #   hs = hs, hsd = hsd
#'   # )
#'   # # Map out the function assignment ============================================
#'   # bessel_map <- list(
#'   #   k1a_shell = c("js", "jsd", "hs", "hsd"),
#'   #   kLa_shell = c("js", "jsd", "jsdd", "ys", "ysd", "ysdd"),
#'   #   kTa_shell = c("js", "jsd", "jsdd", "ys", "ysd", "ysdd"),
#'   #   kTa_fluid = c("js", "jsd", "jsdd", "ys", "ysd", "ysdd"),
#'   #   kLa_fluid = c("js", "jsd", "jsdd", "ys", "ysd", "ysdd"),
#'   #   k3a_fluid = c("js", "jsd")
#'   # )
#'   # # Pre-calculate the Bessel functions and their respective derivatives ========
#'   # bessel_cache <- lapply(row.names(ka_matrix), function(ka_m) {
#'   #   ka_series <- ka_matrix[ka_m, ]
#'   #   bessel_match <- bessel_map[[ka_m]]
#'   #
#'   #   # Calculate for only matching functions
#'   #   lapply(bessel_functions[bessel_match], function(func) {
#'   #     func(m, ka_series)
#'   #   })
#'   # })
#'   # # Add the names ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#'   # names(bessel_cache) <- row.names(ka_matrix)
#'   # bessel_cache
#'   # Define all Bessel functions to apply =======================================
#'   bessel_functions <- list(
#'     js = js_old, jsd = jsd_old, jsdd = jsdd_old,
#'     ys = ys_old, ysd = ysd_old, ysdd = ysdd_old,
#'     hs = hs_old, hsd = hsd_old
#'   )
#'   # Map out the function assignment ============================================
#'   bessel_map <- list(
#'     k1a_shell = c("js", "jsd", "hs", "hsd"),
#'     kLa_shell = c("js", "jsd", "jsdd", "ys", "ysd", "ysdd"),
#'     kTa_shell = c("js", "jsd", "jsdd", "ys", "ysd", "ysdd"),
#'     kTa_fluid = c("js", "jsd", "jsdd", "ys", "ysd", "ysdd"),
#'     kLa_fluid = c("js", "jsd", "jsdd", "ys", "ysd", "ysdd"),
#'     k1a_fluid = c("js", "jsd")
#'   )
#'   # Pre-calculate the Bessel functions and their respective derivatives ========
#'   bessel_cache <- lapply(names(ka_matrix_m), function(ka_m) {
#'     ka_series <- ka_matrix_m[[ka_m]]
#'     bessel_match <- bessel_map[[ka_m]]
#'
#'     # Calculate for only matching functions
#'     lapply(bessel_functions[bessel_match], function(func) {
#'       func(m, ka_series)
#'     })
#'   })
#'   # Add the names ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#'   names(bessel_cache) <- row.names(ka_matrix)
#'   bessel_cache
#' }
#' ################################################################################
#' # Goodman and Stern (1962) modal series solution for elastic-shelled spheres
#' ################################################################################
#' #' Calculates the theoretical TS of an elastic-shelled sphere using the modal
#' #' series solution from Goodman and Stern (1962).
#' #' @param object GAS- or SBF-class object.
#' #' @details
#' #' Calculates the theoretical TS of an elastic-shelled sphere using an exact
#' #' modal series solution
#' #' @return
#' #' Target strength (TS, dB re: 1 m^2)
#' #' @references
#' #' Goodman, R.R., and Stern, R. (1962). Reflection and transmission of sound by
#' #' elastic spherical shells. The Journal of the Acoustical Society of America,
#' #' 34, 338-344.
#' #' @export
#' MSS_goodman_stern <- function(object) {
#'   # Extract model parameters/inputs ============================================
#'   model <- extract(object, "model_parameters")$MSS_goodman_stern
#'   # Get requisite elastic properties ===========================================
#'   G <- model$shell$G
#'   lambda <- model$shell$lambda
#'   density_shell <- model$shell$density
#'   # Extract the required morphometrics =========================================
#'   radius_shell <- model$shell$radius
#'   radius_fluid <- model$fluid$radius
#'   # Get the internal fluid density =============================================
#'   density_fluid <- model$fluid$density
#'   # Calculate shell sound speeds in the longitudinal and transverse directions =
#'   sound_speed_longitudinal <- sqrt((lambda + 2 * G) / density_shell)
#'   sound_speed_transversal <- sqrt(G / density_shell)
#'   # Calculate the associated wavenumbers =======================================
#'   kL <- k(model$parameters$acoustics$frequency, sound_speed_longitudinal)
#'   kT <- k(model$parameters$acoustics$frequency, sound_speed_transversal)
#'   # Calculate the reindexed ka values ==========================================
#'   ka_matrix <- .calculate_ka_matrix(
#'     model$parameters$acoustics$frequency,
#'     model$medium$sound_speed,
#'     model$fluid$sound_speed,
#'     sound_speed_longitudinal,
#'     sound_speed_transversal,
#'     radius_shell,
#'     radius_fluid
#'   )
#'
#'   ANS <- goodman_stern_bm_cpp(
#'     ka_matrix["k1a_shell",],
#'     ka_matrix["kLa_shell",],
#'     ka_matrix["kTa_shell",],
#'     ka_matrix["kLa_fluid",],
#'     ka_matrix["kTa_fluid",],
#'     ka_matrix["k1a_fluid",],
#'     m,
#'     lambda,
#'     G,
#'     model$medium$density / density_shell,
#'     density_fluid / density_shell
#'   )
#'   # Get the modal series iterators =============================================
#'   m_limit <- model$parameters$m_limit
#'   m <- model$parameters$m
#'   # Expand the wavenumber matrix over the modal series limits ==================
#'   ka_matrix_m <- lapply(rownames(ka_matrix), function(ka) {
#'     modal_matrix(ka_matrix[ka, ], m_limit)
#'   })
#'   # Add the names ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#'   names(ka_matrix_m) <- rownames(ka_matrix)
#'   # Pre-calculate and cache the Bessel function outputs ========================
#'   bessel_cache <- .calculate_bessel_cache(ka_matrix_m, m)
#'   # Compute the alpha coefficient matrix =======================================
#'   alpha <- .goodman_stern_alpha(
#'     bessel_cache, ka_matrix_m, m,
#'     lambda, G, model$medium$density,
#'     density_shell, density_fluid
#'   )
#'   # Create the boundary conditions matrices ====================================
#'   boundary_matrices <- .goodman_stern_boundaries(
#'     alpha, ka_matrix, m
#'   )
#'
#'   complex_det <-function(M) {
#'     prod(eigen(M, only.values = TRUE)$values)
#'   }
#'
#'   # Compute the modal series coefficient (b_m) using Cramer's Rule =============
#'   b_m <- lapply(seq_along(boundary_matrices), function(freq_idx) {
#'     vapply(seq_along(m), function(m_idx) {
#'       boundary_m <- boundary_matrices[[freq_idx]][[m_idx]]
#'       # # Calculate determinants using QR decomposition for numerical stability ++
#'       # det_numerator <- prod(
#'       #   abs(Re(diag(qr(boundary_m$A_numerator)$qr)))
#'       # )
#'       # det_denominator <- prod(
#'       #   abs(Re(diag(qr(boundary_m$A_denominator)$qr)))
#'       # )
#'       det_numerator <- complex_det(boundary_m$A_numerator)
#'       det_denominator <- complex_det(boundary_m$A_denominator)
#'
#'       # Apply Cramer's rule: b_m = det(A_numerator) / det(A_denominator) +++++++
#'       if (det_denominator == 0) {
#'         warning(
#'           paste(
#'             "Zero denominator at frequency",
#'             freq_idx, "mode", m[m_idx]
#'           )
#'         )
#'         return(NA)
#'       }
#'
#'       det_numerator / det_denominator
#'     }, FUN.VALUE = complex(1))
#'   })
#'   # Calculate the linear scattering coefficient, f_bs ==========================
#'   f_bs <- -1i / model$parameters$acoustics$k_sw *
#'     vapply(seq_along(b_m), function(ka_idx) {
#'       sum((-1)^m * (2 * m + 1) * b_m[[ka_idx]])
#'     }, FUN.VALUE = complex(1))
#'
#'   # m <- 0:m_limit
#'   # prefactors <- ifelse(m == 0, -1, -1i^m * (2*m + 1))
#'   #
#'   # f_bs <- (1 / model$parameters$acoustics$k_sw) * vapply(seq_along(b_m), function(ka_idx) {
#'   #   sum(prefactors * b_m[[ka_idx]])
#'   # }, FUN.VALUE = complex(1))
#'   # Convert to sigma_bs ========================================================
#'   sigma_bs <- abs(f_bs)^2
#'   result_cpp <- ANS
#'   fbs_from_cpp <- vapply(seq_len(nrow(result_cpp)), function(i) {
#'     -1i / model$parameters$acoustics$k_sw[i] * sum((-1)^m * (2 * m + 1) * result_cpp[i, ])
#'   }, FUN.VALUE = complex(1))
#'   sigma_bs_from_cpp <- abs(fbs_from_cpp)^2
#'
#'   TS <- 10 * log10(sigma_bs)
#'   TS_from_cpp <- 10 * log10(sigma_bs_from_cpp)
#'
#'   plot(frequency * 1e-3, TS, col="red", type='l')
#'   lines(frequency * 1e-3, TS_from_cpp, col="blue")
#'   lines(frequency * 1e-3, shelled_sphere@model$high_pass_stanton$TS, col="gray40")
#'
#'
#'   f_bs_from_cpp <- -1i / model$parameters$acoustics$k_sw *
#'     vapply(seq_along(b_m), function(ka_idx) {
#'       sum((-1)^m * (2 * m + 1) * b_m[[ka_idx]])
#'     }, FUN.VALUE = numeric(1))
#'   # Define MSS slot for ESS-type scatterer =====================================
#'   methods::slot(object, "model")$MSS_goodman_stern <- data.frame(
#'     frequency = model$parameters$acoustics$frequency,
#'     ka_shell = model$parameters$acoustics$k_sw * radius_shell,
#'     ka_fluid = model$parameters$acoustics$k_sw * radius_fluid,
#'     f_bs = f_bs,
#'     sigma_bs = sigma_bs,
#'     TS = db(sigma_bs)
#'   )
#'   # Return object ==============================================================
#'   object
#' }
#' ################################################################################
#' # Kirchoff-Ray Mode approximation
#' ################################################################################
#' #' Calculates the theoretical TS using Kirchoff-ray Mode approximation.
#' #'
#' #' @param object Desired object/animal shape. Must be class "SBF".
#' #' @usage
#' #' KRM(object)
#' #' @details
#' #' Calculates the theoretical TS using the Kirchoff-ray Mode model.
#' #' @return
#' #' Target strength (TS, dB re: 1 m^2)
#' #' @references
#' #' Clay C.S. and Horne J.K. (1994). Acoustic models of fish: The Atlantic cod
#' #' (Gadus morhua). Journal of the Acoustical Society of AMerica, 96, 1661-1668.
#' #' @export
#' KRM <- function(object) {
#'   # Detect object class ========================================================
#'   scatterer_type <- class(object)
#'   # Extract model parameter inputs =============================================
#'   model <- extract(object, "model_parameters")$KRM
#'   # Extract body parameters ====================================================
#'   body <- extract(object, "body")
#'   # Calculate reflection coefficient for medium-body interface =================
#'   R12 <- reflection_coefficient(model$medium, model$body)
#'   # Calculate transmission coefficient and its reverse =========================
#'   T12T21 <- 1 - R12 * R12
#'   # Sum across body position vector ============================================
#'   rpos <- switch(scatterer_type,
#'     FLS = rbind(
#'       x = body$rpos[1, ],
#'       w = c(
#'         body$radius[2],
#'         body$radius[2:(length(body$radius) - 1)],
#'         body$radius[(length(body$radius)) - 1]
#'       ) * 2,
#'       zU = c(
#'         body$radius[2],
#'         body$radius[2:(length(body$radius) - 1)],
#'         body$radius[(length(body$radius)) - 1]
#'       ),
#'       zL = -c(
#'         body$radius[2],
#'         body$radius[2:(length(body$radius) - 1)],
#'         body$radius[(length(body$radius)) - 1]
#'       )
#'     ),
#'     SBF = body$rpos
#'   )
#'   body_rpos_sum <- along_sum(rpos, model$parameters$ns_b)
#'   # Approximate radius of body cylinders =======================================
#'   a_body <- switch(scatterer_type,
#'     FLS = body_rpos_sum[2, ] / 4,
#'     SBF = body_rpos_sum[2, ] / 4
#'   )
#'   # Combine wavenumber (k) and radii to calculate "ka" =========================
#'   ka_body <- matrix(
#'     data = rep(a_body,
#'       each = length(model$parameters$acoustics$k_sw)
#'     ),
#'     ncol = length(a_body),
#'     nrow = length(model$parameters$acoustics$k_b)
#'   ) * model$parameters$acoustics$k_b
#'   # Convert c-z coordinates to required u-v rotated coordinates ================
#'   uv_body <- acousticTS::body_rotation(
#'     body_rpos_sum,
#'     body$rpos,
#'     body$theta,
#'     length(model$parameters$acoustics$k_sw)
#'   )
#'   # Calculate body empirical phase shift function ==============================
#'   body_dorsal_sum <- switch(scatterer_type,
#'     FLS = matrix(
#'       data = rep(body_rpos_sum[3, ],
#'         each = length(model$parameters$acoustics$k_sw)
#'       ),
#'       ncol = length(body_rpos_sum[3, ]),
#'       nrow = length(model$parameters$acoustics$k_sw)
#'     ) / 2,
#'     SBF = matrix(
#'       data = rep(body_rpos_sum[3, ],
#'         each = length(model$parameters$acoustics$k_sw)
#'       ),
#'       ncol = length(body_rpos_sum[3, ]),
#'       nrow = length(model$parameters$acoustics$k_sw)
#'     ) / 2
#'   )
#'   Psi_b <- -pi * model$parameters$acoustics$k_b * body_dorsal_sum /
#'     (2 * (model$parameters$acoustics$k_b * body_dorsal_sum + 0.4))
#'   # Estimate natural log function (phase, etc.) ================================
#'   exp_body <- exp(-2i * model$parameters$acoustics$k_sw * uv_body$vbU) -
#'     T12T21 * exp(-2i * model$parameters$acoustics$k_sw * uv_body$vbU +
#'                    2i * model$parameters$acoustics$k_b *
#'                    (uv_body$vbU - uv_body$vbL) + 1i * Psi_b)
#'   # Resolve summation term =====================================================
#'   body_summation <- sqrt(ka_body) * uv_body$delta_u
#'   # Calculate linear scattering length (m) =====================================
#'   f_body <- rowSums(-((1i * (R12 / (2 * sqrt(pi)))) *
#'                         body_summation * exp_body))
#'   if (scatterer_type == "FLS") {
#'     # Define KRM slot for FLS-type scatterer ===================================
#'     methods::slot(object, "model")$KRM <- data.frame(
#'       frequency = model$parameters$acoustics$frequency,
#'       ka = model$parameters$acoustics$k_sw *
#'         stats::median(a_body, na.rm = TRUE),
#'       f_bs = f_body,
#'       sigma_bs = abs(f_body) * abs(f_body),
#'       TS = 20 * log10(abs(f_body))
#'     )
#'   } else if (scatterer_type == "SBF") {
#'     #### Repeat process for bladder ============================================
#'     # Extract bladder parameters ===============================================
#'     bladder <- acousticTS::extract(object, "bladder")
#'     # Calculate reflection coefficient for bladder =============================
#'     R23 <- acousticTS::reflection_coefficient(
#'       body,
#'       bladder
#'     )
#'     # Sum across body/swimbladder position vectors =============================
#'     bladder_rpos_sum <- acousticTS::along_sum(
#'       bladder$rpos,
#'       model$parameters$ns_sb
#'     )
#'     # Approximate radii of bladder discrete cylinders ==========================
#'     a_bladder <- bladder_rpos_sum[2, ] / 4
#'     # Combine wavenumber (k) and radii to calculate "ka" for bladder ===========
#'     ka_bladder <- matrix(
#'       data = rep(a_bladder, each = length(model$parameters$acoustics$k_sw)),
#'       ncol = length(a_bladder),
#'       nrow = length(model$parameters$acoustics$k_sw)
#'     ) * model$parameters$acoustics$k_sw
#'     # Calculate Kirchoff approximation empirical factor, A_sb ==================
#'     A_sb <- ka_bladder / (ka_bladder + 0.083)
#'     # Calculate empirical phase shift for a fluid cylinder, Psi_p ==============
#'     Psi_p <- ka_bladder / (40 + ka_bladder) - 1.05
#'     # Convert x-z coordinates to requisite u-v rotated coordinates =============
#'     uv_bladder <- acousticTS::bladder_rotation(
#'       bladder_rpos_sum,
#'       bladder$rpos,
#'       bladder$theta,
#'       length(model$parameters$acoustics$k_sw)
#'     )
#'     # Estimate natural log functions  ==========================================
#'     exp_bladder <- exp(-1i * (2 * model$parameters$acoustics$k_b *
#'                                 uv_bladder$v + Psi_p)) * uv_bladder$delta_u
#'     # Calculate the summation term =============================================
#'     bladder_summation <- A_sb * sqrt((ka_bladder + 1) * sin(bladder$theta))
#'     # Estimate backscattering length, f_fluid/f_soft ===========================
#'     f_bladder <- rowSums(-1i * (R23 * T12T21) / (2 * sqrt(pi)) *
#'                            bladder_summation * exp_bladder)
#'     # Estimate total backscattering length, f_bs ===============================
#'     f_bs <- f_body + f_bladder
#'     # Define KRM slot for FLS-type scatterer ===================================
#'     methods::slot(object, "model")$KRM <- data.frame(
#'       frequency = model$parameters$acoustics$frequency,
#'       ka = model$parameters$acoustics$k_sw *
#'         stats::median(a_body, na.rm = TRUE),
#'       f_body = f_body,
#'       f_bladder = f_bladder,
#'       f_bs = f_bs,
#'       sigma_bs = abs(f_bs) * abs(f_bs),
#'       TS = 20 * log10(abs(f_bs))
#'     )
#'   }
#'   # Return object ==============================================================
#'   object
#' }
#' ################################################################################
#' # Primary scattering model for an elastic shelled scatterers (ESS)
#' ################################################################################
#' ################################################################################
#' # Ray-based high pass approximation
#' ################################################################################
#'
#' #' Calculates the theoretical TS of a shelled organism using the non-modal
#' #' High Pass (HP) model
#' #'
#' #' @param object Desired animal object (Elastic Shelled).
#' #' @return
#' #' Target strength (TS, dB re: 1 m^2)
#' #' @references
#' #' Stanton, T.K. (1989). Simple approximate formulas for backscattering of
#' #' sound by spherical and elongated objects. The Journal of the Acoustical
#' #' Society of America, 86, 1499-1510.
#' #'
#' #' @export
#' high_pass_stanton <- function(object) {
#'   # Extract model parameters/inputs ==========================================
#'   model_params <- extract(object, "model_parameters")$high_pass_stanton
#'   medium <- model_params$medium
#'   acoustics <- model_params$acoustics
#'   shell <- model_params$shell
#'   # Multiply acoustic wavenumber by body radius ==============================
#'   k1a <- acoustics$k_sw * shell$radius
#'   # Calculate Reflection Coefficient =========================================
#'   R <- (shell$h * shell$g - 1) / (shell$h * shell$g + 1)
#'   # Calculate backscatter constant, alpha_pi
#'   alpha_pi <- (1 - shell$g * (shell$h * shell$h)) /
#'     (3 * shell$g * (shell$h * shell$h)) +
#'     (1 - shell$g) / (1 + 2 * shell$g)
#'   # Define approximation constants, G_c and F_c ==============================
#'   F_c <- 1
#'   G_c <- 1
#'   # Caclulate numerator term =================================================
#'   num <- ((shell$radius * shell$radius) * (k1a * k1a * k1a * k1a) *
#'             (alpha_pi * alpha_pi) * G_c)
#'   # Caclulate denominator term ===============================================
#'   dem <- (1 + (4 * (k1a * k1a * k1a * k1a) * (alpha_pi * alpha_pi)) /
#'             ((R * R) * F_c))
#'   # Calculate backscatter and return
#'   f_bs <- num / dem
#'   methods::slot(object, "model")$high_pass_stanton <- data.frame(
#'     frequency = acoustics$frequency,
#'     k1a = acoustics$k_sw * shell$radius,
#'     k_s = acoustics$k_b,
#'     f_bs = f_bs,
#'     sigma_bs = abs(f_bs),
#'     TS = 10 * log10(abs(f_bs))
#'   )
#'   object
#' }
#' ################################################################################
#' # Utility scattering models for multiple classes
#' ################################################################################
#' ################################################################################
#' # Fine cylinder modal series solution
#' ################################################################################
#' FCMS <- function(object) {
#'   # Extract model parameters/inputs ============================================
#'   model_params <- extract(object, "model_parameters")$FCMS
#'   parameters <- model_params$parameters
#'   acoustics <- parameters$acoustics
#'   medium <- model_params$medium
#'   body <- model_params$body
#'   # Pull out maximum 'm' specifically ==========================================
#'   m_max <- acoustics$m_max
#'   # Pre-allocate Neumann factors ===============================================
#'   nu <- lapply(
#'     seq_len(nrow(acoustics)),
#'     function(f) neumann(0:acoustics$m_max[f])
#'   )
#'   # Compute kL and ka ==========================================================
#'   k1L <- body$length_body * acoustics$k_sw
#'   k1a <- acoustics$k_sw * sin(body$theta) * body$radius_body
#'   k2a <- acoustics$k_sw * sin(body$theta) / body$h * body$radius_body
#'   # Combine material properties ================================================
#'   gh <- body$g * body$h
#'   # Resolve modal series coefficient calculation method ========================
#'   if (parameters$Bm_method == "Bm_fluid") {
#'     .bm_wrap <- function(f) {
#'       # Get m vector ===========================================================
#'       m <- 0:m_max[f]
#'       # Convert to modal matrices ==============================================
#'       k1a_m <- modal_matrix(k1a[f], max(m))
#'       k2a_m <- modal_matrix(k2a[f], max(m))
#'       # Compute numerator for the boundary coefficient Cm ======================
#'       Cm_num <- (jcd(m, k2a_m)*yc(m, k1a_m)) /
#'         (jc(m, k2a_m)*jcd(m, k1a_m)) -
#'         gh*(ycd(m, k1a_m)/jcd(m, k1a_m))
#'       # Compute denominator for the boundary coefficient Cm ====================
#'       Cm_denom <- (jcd(m, k2a_m)*jc(m, k1a_m)) /
#'         (jc(m, k2a_m)*jcd(m, k1a_m)) - gh
#'       # Resolve Cm =============================================================
#'       Cm <- Cm_num / Cm_denom
#'       # Return Bm
#'       1i^(m+1) * (-nu[[f]] * 1i^(m) / (1 + 1i * Cm))
#'     }
#'   } else if (parameters$Bm_method == "Bm_rigid") {
#'     .bm_wrap <- function(f) {
#'       # Get m vector ===========================================================
#'       m <- 0:m_max[f]
#'       # Convert to modal matrix ================================================
#'       k1a_m <- modal_matrix(k1a[f], max(m))
#'       # Return Bm ==============================================================
#'       (-1)^m * nu[[f]] * (jcd(m, k1a_m) / hcd(m, k1a_m))
#'     }
#'   } else {
#'     .bm_wrap <- function(f) {
#'       # Get m vector ===========================================================
#'       m <- 0:m_max[f]
#'       # Convert to modal matrix ================================================
#'       k1a_m <- modal_matrix(k1a[f], max(m))
#'       # Return Bm ==============================================================
#'       (-1)^m * nu[[f]] * (jc(m, k1a_m) / hc(m, k1a_m))
#'     }
#'   }
#'   fc_bm <- lapply(seq_len(nrow(acoustics)), FUN = function(f) .bm_wrap(f))
#'   # Compute the linear scattering coefficient, f_bs ============================
#'   f_bs <- vapply(
#'     seq_len(nrow(acoustics)),
#'     FUN = function(f) {
#'       if (parameters$Bm_method == "Bm_fluid") {
#'         -body$length_body / pi *
#'           (sin(k1L[f] * cos(body$theta_body)) /
#'              (k1L[f] * cos(body$theta_body))) *
#'           sum(fc_bm[[f]])
#'       } else {
#'         1i * body$length_body / pi *
#'           (sin(k1L[f] * cos(body$theta_body)) /
#'              (k1L[f] * cos(body$theta_body))) *
#'           sum(fc_bm[[f]])
#'       }
#'     },
#'     FUN.VALUE = complex(1)
#'   )
#'   # Calculate backscatter and return ===========================================
#'   # Compute sigma_bs +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#'   o_bs <- abs(f_bs)
#'   methods::slot(object, "model")$FCMS <- data.frame(
#'     frequency = acoustics$frequency,
#'     f_bs = f_bs,
#'     sigma_bs = o_bs,
#'     TS = 20 * log10(o_bs)
#'   )
#'   object
#' }
#'
#' ################################################################################
#' # Sphere modal series solution
#' ################################################################################
#' SPHMS <- function(object) {
#'   # Extract model parameters/inputs ============================================
#'   model_params <- extract(object, "model_parameters")$SPHMS
#'   parameters <- model_params$parameters
#'   acoustics <- parameters$acoustics
#'   medium <- model_params$medium
#'   body <- model_params$body
#'   # Pull out maximum 'm' specifically ==========================================
#'   m_limit <- acoustics$m_limit
#'   # Compute Delta (shell thickness) =============================================
#'   b_delta <- body$radius_shell - body$radius_fluid
#'   # Compute ka and kb ==========================================================
#'   k1a <- acoustics$k_sw * body$radius_shell
#'   k2a <- acoustics$k_shell * body$radius_shell
#'   k2b <- acoustics$k_shell * body$radius_fluid
#'   k3b <- acoustics$k_fluid * body$radius_fluid
#'   # Define the boundary expansion coefficient method ===========================
#'   .bm_wrap <- switch(
#'     parameters$Bm_method,
#'     Bm_rigid = function(f) .sphms_bm_rigid(k1a[f], m_limit[f]),
#'     Bm_pressure_release = function(f) .sphms_bm_prelease(k1a[f], m_limit[f]) ,
#'     Bm_fluid = function(f) {
#'       .sphms_bm_fluid(k1a[f], k2a[f], body$g31, body$h31, m_limit[f])
#'     },
#'     Bm_shelled_pressure_release = function(f) {
#'       .sphms_bm_shelled_prelease(
#'         k1a[f], k2a[f], k2b[f], body$g31, body$h31, m_limit[f]
#'       )
#'     },
#'     Bm_shelled_fluid = function(f) {
#'       .sphms_bm_shelled_fluid(
#'         k1a[f], k2a[f], k2b[f], k3b[f], body$g21, body$g32, body$h21, body$h32,
#'         m_limit[f]
#'       )
#'     }
#'   )
#'   sph_bm <- vapply(seq_len(nrow(acoustics)),
#'                    FUN = function(f) .bm_wrap(f),
#'                    FUN.VALUE = complex(1))
#'   # Compute the linear scattering coefficient, f_bs ============================
#'   f_bs <- -1i / acoustics$k_sw * sph_bm
#'   # Calculate backscatter and return ===========================================
#'   # Compute sigma_bs +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#'   o_bs <- abs(f_bs)
#'   methods::slot(object, "model")$SPHMS <- data.frame(
#'     frequency = acoustics$frequency,
#'     f_bs = f_bs,
#'     sigma_bs = o_bs,
#'     TS = 20 * log10(o_bs)
#'   )
#'   object
#' }
#'
#' .sphms_bm_rigid <- function(k1a, m_limit) {
#'   # Expand k1a into matrix =====================================================
#'   k1a_m <- modal_matrix(k1a, m_limit)
#'   # Expand modal series iterator ===============================================
#'   m <- 0:m_limit
#'   # Calculate the expansion coefficient, A_m ===================================
#'   Am <- -(jsd(m_limit, k1a_m) / hsd(m_limit, k1a_m))
#'   # Return the entire boundary modal term ======================================
#'   sum((2 * m + 1) * (-1)^m * Am)
#' }
#'
#' .sphms_bm_prelease <- function(k1a, m_limit) {
#'   # Expand k1a into matrix =====================================================
#'   k1a_m <- modal_matrix(k1a, m_limit)
#'   # Expand modal series iterator ===============================================
#'   m <- 0:m_limit
#'   # Calculate the expansion coefficient, A_m ===================================
#'   Am <- -(js(m_limit, k1a_m) / hs(m_limit, k1a_m))
#'   # Return the entire boundary modal term ======================================
#'   sum((2 * m + 1) * (-1)^m * Am)
#' }
#'
#' .sphms_bm_fluid <- function(k1a, k2a, g31, h31, m_limit) {
#'   # Expand ka into matrices ====================================================
#'   k1a_m <- modal_matrix(k1a, m_limit)
#'   k3a_m <- modal_matrix(k3a, m_limit)
#'   # Expand modal series iterator ===============================================
#'   m <- 0:m_limit
#'   # Get material properties product ============================================
#'   gh <- g31 * h31
#'   # Calculate expansion ceofficient, C_m =======================================
#'   # Numerator ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#'   cm_num <- (jsd(m, k3a_m) * ys(m, k1a_m)) / (js(m, k3a_m) * jsd(m, k1a_m)) -
#'     (gh * ysd(m, k1a_m) / jsd(m, k1a_m))
#'   # Denominator ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#'   cm_denom <- (jsd(m, k3a_m) * js(m, k1a_m)) / (js(m, k3a_m) * jsd(m, k1a_m)) -
#'     gh
#'   # Quotient +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#'   cm <- cm_num / cm_denom
#'   # Calculate the expansion coefficient, A_m ===================================
#'   Am <- (-1 / (1 + 1i * cm))
#'   # Return the entire boundary modal term ======================================
#'   sum((2 * m + 1) * (-1)^m * Am)
#' }
#'
#' .sphms_bm_shelled_prelease <- function(k1a, k2a, k2b, g31, h31, m_limit) {
#'   # Expand ka into matrices ====================================================
#'   k1a_m <- modal_matrix(k1a, m_limit)
#'   k2a_m <- modal_matrix(k2a, m_limit)
#'   k2b_m <- modal_matrix(k2b, m_limit)
#'   # Expand modal series iterator ===============================================
#'   m <- 0:m_limit
#'   # Get material properties product ============================================
#'   gh <- g31 * h31
#'   # Calculate the expansion coefficient, A_m ===================================
#'   # Compute elements for boundary matrix +++++++++++++++++++++++++++++++++++++++
#'   a11 <- -hs(m, k1a_m)
#'   a12 <- js(m, k1a_m)
#'   a21 <- -g21 * h21 * hsd(m, k1a_m)
#'   b1 <- js(m, k1a_m)
#'   b2 <- jsd(m, k1a_m) * g21 * h21
#'   d1 <- js(m, k2a_m) * ys(m, k2b_m) - js(m, k2b_m) * ys(m, k2a_m)
#'   d2 <- jsd(m, k2a_m) * ys(m, k2b_m) - js(m, k2b_m) * ysd(m, k2a_m)
#'   # Calculate coefficient ++++++++++++++++++++++++++++++++++++++++++++++++++++++
#'   Am <- (b1 * d2 - d1 * b2) / (a11 * d2 - d1 * a21)
#'   # Return the entire boundary modal term ======================================
#'   sum((2 * m + 1) * (-1)^m * Am)
#' }
#'
#' .sphms_bm_shelled_fluid <- function(
#'     k1a, k2a, k2b, k3b, g21, g32, h21, h32, m_limit
#' ) {
#'   # Expand modal series iterator ===============================================
#'   m <- 0:m_limit
#'   # Expand ka into matrices ====================================================
#'   k1a_m <- modal_matrix(k1a, m_limit)
#'   k2a_m <- modal_matrix(k2a, m_limit)
#'   k2b_m <- modal_matrix(k2b, m_limit)
#'   k3b_m <- modal_matrix(k3b, m_limit)
#'   # Get material properties product ============================================
#'   gh <- g31 * h31
#'   # Calculate the expansion coefficient, A_m ===================================
#'   # Compute elements for boundary matrix +++++++++++++++++++++++++++++++++++++++
#'   a11 <- -hs(m, k1a_m)
#'   a12 <- js(m, k2a_m)
#'   a13 <- ys(m, k2a_m)
#'   a21 <- -g21 * h21 * hsd(m, k1a_m)
#'   a22 <- jsd(m, k2a_m)
#'   a23 <- ysd(m, k2a_m)
#'   a32 <- js(m, k2b_m) * jsd(m, k3b_m) - g32 * h32 * jsd(m, k2b_m) * js(m, k3b_m)
#'   a33 <- ys(m, k2b_m) * jsd(m, k3b_m) - g32 * h32 * ysd(m, k2b_m) * js(m, k3b_m)
#'   b1 <- js(m, k1a_m)
#'   b2 <- jsd(m, k1a_m) * g21 * h21
#'   # Calculate coefficient ++++++++++++++++++++++++++++++++++++++++++++++++++++++
#'   Am <- (b1*a22*a33 + a13*b2*a32 - a12*b2*a33 - b1*a23*a32) /
#'     (a11*a22*a33 + a13*a21*a32 - a12*a21*a33 - a11*a23*a32)
#'   # Return the entire boundary modal term ======================================
#'   sum((2 * m + 1) * (-1)^m * Am)
#' }
#'
#' ################################################################################
#' # Prolate spheroid modal series solution
#' ################################################################################
#' #' Calculates the theoretical target strength of a prolate spheroid using a
#' #' modal series solution
#' #'
#' #' @param object Prolate spheroid object.
#' #' @return
#' #' Target strength (TS, dB re: 1 m^2)
#' #' @references
#' #' Furusawa, M. (1988). Prolate spheroidal models for predicting general trends
#' #' of fish target strength. Journal of the Acoustical Society of Japan, 9:
#' #' 13-24.
#' #'
#' #' @export
#' PSMS <- function(object) {
#'   # Extract model parameters/inputs ============================================
#'   model_params <- extract(object, "model_parameters")$PSMS
#'   parameters <- model_params$parameters
#'   acoustics <- parameters$acoustics
#'   medium <- model_params$medium
#'   body <- model_params$body
#'   # Pre-allocate Neumann factors ===============================================
#'   nu <- lapply(
#'     seq_len(nrow(acoustics)),
#'     function(f) neumann(0:acoustics$m_max[f])
#'   )
#'   # Pre-allocate azimuthal phase relations =====================================
#'   azimuth <- lapply(
#'     seq_len(nrow(acoustics)),
#'     function(f) {
#'       cos(0:acoustics$m_max[f] * (body$phi_body - body$phi_scatter))
#'     }
#'   )
#'   # Compute the expansion matrix Amn ===========================================
#'   if (parameters$Amn_method == "Amn_liquid_simplify") {
#'     .amn_wrap <- function(f) {
#'       # Compute using the simplified expression from Eq. 5 (Furusawa, 1988) ++++
#'       liquid_spheroidal_simplified_expansion(
#'         acoustics$m_max[f], acoustics$n_max[f], acoustics$chi_sw[f],
#'         acoustics$chi_body[f], body$xi, medium$density, body$density
#'       )$amn
#'     }
#'   } else if (parameters$Amn_method == "Amn_liquid") {
#'     .amn_wrap <- function(f) {
#'       # Compute the kernel matrices for a liquid-filled spheroidal scatterer +++
#'       kernels <- liquid_spheroidal_kernels(
#'         acoustics$m_max[f], acoustics$n_max[f], acoustics$chi_sw[f],
#'         acoustics$chi_body[f], body$theta_body, body$xi,
#'         body$density, medium$density
#'       )
#'       # Solve for the expansion coefficient Amn ++++++++++++++++++++++++++++++++
#'       solve_liquid_spheroidal_Amn(
#'         kernels$K1_kernel, kernels$K3_kernel
#'       )
#'     }
#'   } else if (parameters$Amn_method == "Amn_pressure_release") {
#'     .amn_wrap <- function(f) {
#'       # Compute the external radial incident wave function +++++++++++++++++++++
#'       R_incident <- radial_external_incoming_matrix(
#'         acoustics$m_max[f],
#'         acoustics$n_max[f], acoustics$chi_sw[f], body$xi
#'       )
#'       # Compute the external radial scattering wave function +++++++++++++++++++
#'       R_scattering <- radial_external_scattering_matrix(
#'         acoustics$m_max[f],
#'         acoustics$n_max[f], acoustics$chi_sw[f], body$xi
#'       )
#'       # Compute the expansion matrix +++++++++++++++++++++++++++++++++++++++++++
#'       - R_incident$value / R_scattering$value
#'     }
#'   } else {
#'     .amn_wrap <- function(f) {
#'       # Compute the first deriv of the external radial incident wave function ++
#'       dR_incident <- radial_external_incoming_matrix(
#'         acoustics$m_max[f],
#'         acoustics$n_max[f], acoustics$chi_sw[f], body$xi
#'       )
#'       # Compute the first deriv of the external radial scattering wave function
#'       dR_scattering <- radial_external_scattering_matrix(
#'         acoustics$m_max[f],
#'         acoustics$n_max[f], acoustics$chi_sw[f], body$xi
#'       )
#'       # Compute the expansion matrix +++++++++++++++++++++++++++++++++++++++++++
#'       - dR_incident$derivative / dR_scattering$derivative
#'     }
#'   }
#'   ps_amn <- lapply(
#'     seq_len(nrow(acoustics)),
#'     FUN = function(f) .amn_wrap(f)
#'   )
#'   # Compute the linear scattering coefficient, f_bs ============================
#'   f_bs <- vapply(
#'     seq_len(nrow(acoustics)),
#'     FUN = function(f){
#'       .psms_fbs(
#'         nu[[f]], azimuth[[f]], acoustics$m_max[f], acoustics$n_max[f],
#'         acoustics$chi_sw[f],
#'         body$theta_body, body$theta_scatter, ps_amn[[f]],
#'         ifelse(
#'           parameters$Amn_method %in% c("Amn_liquid", "Amn_liquid_simplified"),
#'           TRUE,
#'           FALSE
#'           )
#'       )
#'     },
#'     FUN.VALUE = complex(1)
#'   )
#'   # Calculate backscatter and return ===========================================
#'   # Compute sigma_bs +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#'   o_bs <- abs(-2i / acoustics$k_sw * f_bs)
#'   methods::slot(object, "model")$PSMS <- data.frame(
#'     frequency = acoustics$frequency,
#'     f_bs = f_bs,
#'     sigma_bs = o_bs,
#'     TS = 20 * log10(o_bs)
#'   )
#'   object
#' }
#'
#' .psms_fbs <- function(nu, azimuth, m_max, n_max, chi_sw, theta_body,
#'                       theta_scatter, ps_amn, list2mat = FALSE) {
#'   # Pre-compute the product of the regular and scatter Smn matrices ============
#'   smn_matrix <- outer(
#'     0:m_max, 0:n_max,
#'     Vectorize(function(m, n) {
#'       if (n < m) return(NA)
#'       # Regular
#'       Smn(m, n, chi_sw, cos(theta_body), normalize = TRUE)$value *
#'         # Scattering
#'         Smn(m, n, chi_sw, cos(theta_scatter), normalize = TRUE)$value
#'     })
#'   )
#'
#'   # Convert list to matrix if defined ==========================================
#'   if (list2mat) {
#'     # Get the matrix dimensions and pre-allocated 'Amn_matrix' +++++++++++++++++
#'     smn_matrix_dims <- dim(smn_matrix)
#'     nrows <- smn_matrix_dims[1]
#'     ncols <- smn_matrix_dims[2]
#'     Amn_matrix <- matrix(0 + 0i,
#'                          nrow = nrows,
#'                          ncol = ncols)
#'     # Reshape Amn list into matrix +++++++++++++++++++++++++++++++++++++++++++++
#'     seq_dim <- seq_len(min(length(ps_amn), nrows))
#'     Amn_rows <- vapply(
#'       seq_dim, FUN = function(r) {
#'         vec <- complex(length = ncols)
#'         am_vec <- as.vector(ps_amn[[r]])
#'         if (length(am_vec) > 0L) {
#'           end_col <- min(smn_matrix_dims[2], r + length(am_vec) - 1L)
#'           cols <- r:end_col
#'           vec[cols] <- am_vec[seq_len(length(cols))]
#'         }
#'         vec
#'       }, FUN.VALUE = complex(smn_matrix_dims[2]), USE.NAMES = FALSE
#'     )
#'     Amn_matrix[seq_dim, ] <- t(Amn_rows)
#'   } else {
#'     Amn_matrix <- ps_amn
#'   }
#'   # Calculate the linear scattering coefficient ================================
#'   sum(nu * smn_matrix * Amn_matrix * azimuth, na.rm = TRUE)
#' }
#'
#' ################################################################################
#' # Model registry/API
#' ################################################################################
#' #' Model registry: maps model names to their functions and related initializers
#' #' @export
# model_registry <- list(
#   # DCM = DCM,
#   DWBA = DWBA,
#   DWBA_curved = DWBA_curved,
#   SDWBA = SDWBA,
#   SDWBA_curved = SDWBA_curved,
#   calibration = calibration,
#   ESSMS = ESSMS,
#   SPHMS = SPHMS,
#   KRM = KRM,
#   high_pass_stanton = high_pass_stanton,
#   PSMS = PSMS,
#   FCMS = FCMS
# )
