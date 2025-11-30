#' Calculates the theoretical TS of a fluid-like scatterer at a given frequency
#' using the stochastic distorted Born wave approximation (DWBA) model.
#'
#' @param object FLS-class scatterer.
#' @references
#' Demer, D.A., and Conti, S.G. 2003. Reconciling theoretical versus empirical
#' target strengths of krill: effects of phase variability on the distorted-wave
#' Born approximation. ICES J. Mar. Sci., 60: 429-434.
#' @export
SDWBA <- function(object) {
  # Extract model parameters/inputs ============================================
  model <- acousticTS::extract(object, "model_parameters")$SDWBA
  body <- acousticTS::extract(object, "body")
  theta <- body$theta
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
  SDWBA_resampled <- function(i) {
    # Parse phase-specific parameters ==========================================
    sub_params <- model$parameters[[i]]
    # Calculate rotation matrix and update wavenumber matrix ===================
    k_sw_rot <- sub_params$acoustics$k_sw %*% rotation_matrix
    # Calculate Euclidean norms ================================================
    k_sw_norm <- vecnorm(k_sw_rot)
    # Update position matrices  ================================================
    r0 <- rbind(
      sub_params$body_params$rpos[1:3, ],
      sub_params$body_params$radius
    )
    # Calculate position matrix lags  ==========================================
    r0_diff <- t(diff(t(r0)))
    # Multiply wavenumber and body matrices ====================================
    r0_diff_k <- t(vapply(seq_along(k_sw_norm),
                          FUN = function(x) {
                            colSums(r0_diff[1:3, ] *
                                      k_sw_rot[x, ])
                          },
                          FUN.VALUE = numeric(ncol(r0_diff))
    ))
    # Calculate Euclidean norms ================================================
    r0_diff_norm <- sqrt(colSums(r0_diff[1:3, ] * r0_diff[1:3, ]))
    # Estimate angles between body cylinders ===================================
    alpha <- acos(r0_diff_k / (k_sw_norm %*% t(r0_diff_norm)))
    beta <- abs(alpha - pi / 2)
    # Call in metrics ==========================================================
    phase_sd <- sub_params$meta_params$phase_sd
    # r0_diff_h <- r0_diff / h
    # r0_h <- r0 / h
    # Define integrand =========================================================
    integrand <- function(s, x, y) {
      # integrand <- function( s , x ) {
      # rint_mat <- s * r0_diff_h[ , y ] + r0_h[ , y ]
      rint_mat <- s * r0_diff[, y] + r0[, y]
      # rint_k1_h_mat <- k_sw_rot[ x , ] %*% rint_mat[ 1 : 3 ]
      rint_k1_h_mat <- k_sw_rot[x, ] %*% rint_mat[1:3] / h
      bessel <- jc(1, 2 * (k_sw_norm[x] * rint_mat[4] / h *
                             cos(beta[x, y]))) / cos(beta[x, y])
      # bessel <- jc( 1 , 2 * ( k_sw_norm[ x ] * rint_mat[ 4 ] *
      #                           cos( beta[ x , y ] ) ) ) / cos( beta[ x , y] )
      # fb_a <- k_sw_norm[ x ] / 4 * R * rint_mat[ 4 ] * h *
      #   exp( 2i * rint_k1_h_mat ) * bessel[y] * r0_diff_norm[ y ]
      fb_a <- k_sw_norm[x] / 4 * R * rint_mat[4] *
        exp(2i * rint_k1_h_mat) * bessel * r0_diff_norm[y]
      sum(fb_a, na.rm = TRUE)
    }
    # Vectorize integrand function =============================================
    integrand_vectorized <- Vectorize(integrand)
    stochastic_TS <- function(n_k, n_segments, n_iterations) {
      phase_cyl <- vapply(seq_len(n_k),
                          FUN = function(x) {
                            cyl_phase <- array(
                              vapply(seq_len(n_segments),
                                     FUN = function(y) {
                                       phase_integrate(
                                         x, y,
                                         n_iterations,
                                         integrand_vectorized,
                                         phase_sd
                                       )
                                     },
                                     FUN.VALUE = complex(n_iterations)
                              ),
                              dim = c(n_iterations, n_segments)
                            )
                            cyl_sum_phase <- rowSums(cyl_phase, na.rm = TRUE)
                            cyl_sum_phase
                          },
                          FUN.VALUE = complex(n_iterations)
      )
      phase_cyl <- array(phase_cyl, dim = c(n_iterations, n_k))
      data.frame(
        f_bs = colMeans(phase_cyl),
        sigma_bs = colMeans(sigma_bs(phase_cyl)),
        TS_mean = db(colMeans(sigma_bs(phase_cyl))),
        TS_sd = db(apply(sigma_bs(phase_cyl), 2, stats::sd))
      )
    }
    # Calculate linear scatter response ========================================
    backscatter_df <- stochastic_TS(
      n_k = length(k_sw_norm),
      n_segments = sub_params$n_segments - 1,
      n_iterations = sub_params$meta_params$n_iterations
    )
    backscatter_df
  }
  # Generate results dataframe by collating resampled results ==================
  results <- do.call(
    "rbind",
    lapply(seq_along(model$parameters),
           FUN = function(i) SDWBA_resampled(i)
    )
  )
  # Update scatterer object ====================================================
  methods::slot(object, "model")$SDWBA$f_bs <- results$f_bs
  methods::slot(object, "model")$SDWBA$sigma_bs <- results$sigma_bs
  methods::slot(object, "model")$SDWBA$TS <- results$TS_mean
  methods::slot(object, "model")$SDWBA$TS_sd <- results$TS_sd
  # Return object ==============================================================
  object
}
#' Calculates the theoretical TS of a fluid-like scatterer at a given frequency
#' using the stochastic distorted Born wave approximation (DWBA) model using
#' Eq. (6) from Stanton et al. (1998).
#' @param object FLS-class scatterer.
#' @references
#' Stanton, T.K., Chu, D., and Wiebe, P.H. 1998. Sound scattering by several
#' zooplankton groups. II. Scattering models. J. Acoust. Soc. Am., 103, 236-253.
#' @rdname SDWBA_curved
#' @export
SDWBA_curved <- function(object) {
  # Extract model parameters/inputs ============================================
  model <- acousticTS::extract(object, "model_parameters")$SDWBA_curved
  body <- acousticTS::extract(object, "body")
  # Parse position matrices ====================================================
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
  SDWBA_resampled_c <- function(i) {
    # Parse phase-specific parameters ==========================================
    sub_params <- model$parameters[[i]]
    # Calculate rotation matrix and update wavenumber matrix ===================
    k_sw_rot <- sub_params$acoustics$k_sw %*% rotation_matrix
    # Calculate Euclidean norms ================================================
    k_sw_norm <- vecnorm(k_sw_rot)
    # Update position matrices  ================================================
    r0 <- rbind(
      sub_params$body_params$rpos[1:3, ],
      sub_params$body_params$radius
    )
    # Calculate position matrix lags  ==========================================
    r0_diff <- t(diff(t(r0)))
    # Multiply wavenumber and body matrices ====================================
    r0_diff_k <- t(vapply(seq_along(k_sw_norm),
                          FUN = function(x) {
                            colSums(r0_diff[1:3, ] *
                                      k_sw_rot[x, ])
                          },
                          FUN.VALUE = numeric(ncol(r0_diff))
    ))
    # Calculate Euclidean norms ================================================
    r0_diff_norm <- sqrt(colSums(r0_diff[1:3, ] * r0_diff[1:3, ]))
    # Estimate angles between body cylinders ===================================
    alpha <- acos(r0_diff_k / (k_sw_norm %*% t(r0_diff_norm)))
    beta <- abs(alpha - pi / 2)
    # Call in metrics ==========================================================
    phase_sd <- sub_params$meta_params$phase_sd
    r0_diff_h <- r0_diff / h
    r0_h <- r0 / h
    # Define integrand =========================================================
    integrand_c <- function(s, x, y) {
      rint_mat <- s * r0_diff_h[, y] + r0_h[, y]
      rint_k1_h_mat <- k_sw_rot[x, ] %*% rint_mat[1:3]
      bessel <- jc(1, 2 * (k_sw_norm[x] * rint_mat[4] *
                             cos(beta[x, y]))) / cos(beta[x, y])
      fb_a <- k_sw_norm[x] / 4 * R * rint_mat[4] * h *
        exp(2i * rint_k1_h_mat) * bessel * r0_diff_norm[y]
      sum(fb_a, na.rm = TRUE)
    }
    # Vectorize integrand function =============================================
    integrand_vectorized_c <- Vectorize(integrand_c)
    stochastic_TS_c <- function(n_k, n_segments, n_iterations) {
      phase_cyl <- vapply(seq_len(n_k),
                          FUN = function(x) {
                            cyl_phase <- array(
                              vapply(seq_len(n_segments),
                                     FUN = function(y) {
                                       phase_integrate(
                                         x, y,
                                         n_iterations,
                                         integrand_vectorized_c,
                                         phase_sd
                                       )
                                     },
                                     FUN.VALUE = complex(n_iterations)
                              ),
                              dim = c(n_iterations, n_segments)
                            )
                            cyl_sum_phase <- rowSums(cyl_phase, na.rm = TRUE)
                            cyl_sum_phase
                          },
                          FUN.VALUE = complex(n_iterations)
      )
      phase_cyl <- array(phase_cyl, dim = c(n_iterations, n_k))
      data.frame(
        f_bs = colMeans(phase_cyl),
        sigma_bs = colMeans(sigma_bs(phase_cyl)),
        TS_mean = db(colMeans(sigma_bs(phase_cyl))),
        TS_sd = db(apply(sigma_bs(phase_cyl), 2, stats::sd))
      )
    }
    # Calculate linear scatter response ========================================
    backscatter_df <- stochastic_TS_c(
      n_k = length(k_sw_norm),
      n_segments = sub_params$n_segments - 1,
      n_iterations = sub_params$meta_params$n_iterations
    )
    backscatter_df
  }
  # Generate results dataframe by collating resampled results ==================
  results <- do.call(
    "rbind",
    lapply(seq_along(model$parameters),
           FUN = function(i) SDWBA_resampled_c(i)
    )
  )
  # Update scatterer object ====================================================
  methods::slot(object, "model")$SDWBA_curved$f_bs <- results$f_bs
  methods::slot(object, "model")$SDWBA_curved$sigma_bs <- results$sigma_bs
  methods::slot(object, "model")$SDWBA_curved$TS <- results$TS_mean
  methods::slot(object, "model")$SDWBA_curved$TS_sd <- results$TS_sd
  # Return object ==============================================================
  object
}
