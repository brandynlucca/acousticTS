################################################################################
# Primary scattering models for fluid-like scatterers (FLS)
################################################################################
#' Calculates the theoretical TS of a target using the deformed cylinder
#' model (DCM).
#'
#' @param object FLS-class scatterer.
#' @usage
#' DCM(object)
#' @return
#' Target strength (TS, dB re. 1 \eqn{m^2}) of a FLS-class object.
#' #' @references
#' Stanton, T.K., Chu, D., and Wiebe, P.H. 1998. Sound scattering by several
#' zooplankton groups. II. Scattering models. Journal of the Acoustical Society
#' of America, 103(1), 236-253.
#' @export
DCM <- function(object){
  # Extract model parameters/inputs ============================================
  model_params <- extract(object, "model_parameters")$DCM
  medium <- model_params$medium
  acoustics <- model_params$acoustics
  body <- extract(object, "model_parameters")$DCM$body
  # Multiply acoustic wavenumber by body radius ================================
  k1a <- acoustics$k_sw * body$radius
  k2a <- acoustics$k_b * body$radius
  # Calculate Reflection Coefficient ===========================================
  R <- R12(body, medium)
  # Calculate Transmission Coefficients ========================================
  T12 <- 1 - R; T21 <- 1 + R
  # Calculate phase shift term for the ray model ===============================
  mu <- (-pi / 2 * k1a) / (k1a + 0.4)
  # Calculate Interference Term between echoes from front/back interfaces ======
  I0 <- 1 - T12 * T21 * exp(4i * k2a) * exp(1i * mu * k1a)
  # Calculate half-width of directivity pattern ================================
  w12 <- sqrt(body$radius_curvature * body$radius)
  # Apply shift to orientation required for this model =========================
  theta_shift <- body$theta - pi / 2
  # Calculate Taylor expansion of directivity function =========================
  D_theta <- exp(-body$alpha_B *
                   (2 * theta_shift * body$radius_curvature / body$length)^2)
  # Calculate linear backscatter function, f_bs ================================
  f_bs <- 0.5 * w12 * R * exp(-2i * k1a) * I0 * D_theta
  # Calculate TS ===============================================================
  slot(object, "model")$DCM <- data.frame(f_bs = f_bs,
                                          sigma_bs = abs(f_bs)^2,
                                          TS = 20 * log10(abs(f_bs)))
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
#' @import stats
#' @export
DWBA <- function(object){
  # Extract model parameters/inputs ============================================
  model_params <- extract(object, "model_parameters")$DWBA
  medium <- model_params$medium
  acoustics <- model_params$parameters$acoustics
  body <- extract(object, "body")
  # Material properties calculation ============================================
  R <- 1 / (body$g * body$h^2) + 1 / body$g - 2
  # Calculate rotation matrix and update wavenumber matrix =====================
  rotation_matrix <- matrix(c(cos(body$theta), sin(body$theta), 0), 1)
  k_sw_rotation <- acoustics$k_sw %*% rotation_matrix
  # Calculate Euclidean norms ==================================================
  k_sw_norm <- vecnorm(k_sw_rotation)
  # Update position matrices  ==================================================
  rpos <- rbind(body$rpos, a=body$radius)
  # Calculate position matrix lags  ============================================
  rpos_diff <- t(diff(t(rpos)))
  # Multiply wavenumber and body matrices ======================================
  rpos_diff_k <- t(sapply(1:length(k_sw_norm),
                          FUN = function(x) {
                            colSums(rpos_diff[1:3, ] * k_sw_rotation[x, ])
                          }))
  # Calculate Euclidean norms ==================================================
  rpos_diff_norm <- sqrt(colSums(rpos_diff[1:3, ]^2))
  # Estimate angles between body cylinders =====================================
  alpha <- acos(rpos_diff_k / (k_sw_norm %*% t(rpos_diff_norm)))
  beta <- abs(alpha - pi / 2)
  # Draw phase deviation from Normal distributions =============================
  # Define integrand ===========================================================
  integrand <- function(s, x) {
    rint_mat <- s * rpos_diff + rpos[, 1:ncol(rpos_diff)]
    rint_k1_h_mat <- k_sw_rotation[x, ] %*% rint_mat[1:3, ] / body$h
    bessel <- ja(1, 2 * (k_sw_norm[x] * rint_mat[4, ] / body$h * cos(beta[x, ]))) / cos(beta[x, ])
    fb_a <- k_sw_norm[x] / 4 * R * rint_mat[4, ] * exp(2i * rint_k1_h_mat) * bessel * rpos_diff_norm
    return(sum(fb_a, na.rm=T))
  }
  # Vectorize integrand function ===============================================
  integrand_vec <- Vectorize(integrand)
  # Integrate over position matrix =============================================
  sigma <- sapply(1:length(k_sw_norm),
                  FUN = function(x) {
                    # Real =====================================================
                    Ri <- integrate(function(s) Re(integrand_vec(s, x)), 0, 1)$value
                    # Real =====================================================
                    Ii <- 1i * integrate(function(s) Im(integrand_vec(s, x)), 0, 1)$value
                    return(Ri + Ii)
                  })
  # Return object ==============================================================
  slot(object, "model")$DWBA <- data.frame(f_bs = sigma,
                                           sigma_bs = sigma^2,
                                           TS = 10 * log10(abs(sigma)^2))
  return(object)
}
################################################################################
# Primary scattering model for an elastic solid sphere (CAL)
################################################################################
#' Calculates theoretical TS of a solid sphere of a certain material at a given
#' frequency.
#'
#' @description
#' This function is a wrapper around TS_calculate(...) that parametrizes the
#' remainder of the model, while also doing simple calculations that do not
#' need to be looped. This function provides a TS estimate at a given frequency.
#' @param object CAL-class object.
#' @usage
#' calibration(object)
#' @return The theoretical acoustic target strength (TS, dB re. 1 \eqn{m^2})  of a
#' solid sphere at a given frequency.
#' @references
#' MacLennan D. N. (1981). The theory of solid spheres as sonar calibration
#' targets. Scottish Fisheries Research No. 22, Department of Agriculture and
#' Fisheries for Scotland.
#' @export
calibration <- function(object) {
  # Extract model parameters/inputs ============================================
  model_params <- extract(object, "model_parameters")$calibration
  medium <- model_params$medium
  acoustics <- model_params$parameters$acoustics
  model_body <- model_params$scatterer$body
  # Equations 6a ===============================================================
  ka_sw <- acoustics$k_sw * model_body$radius #medium
  ka_l <- acoustics$k_l * model_body$radius #longitudinal
  ka_t <- acoustics$k_t * model_body$radius #transversal
  # Set limit for iterations ===================================================
  limit <- round(max(ka_sw)) + 10
  m <- 0:limit
  # Convert ka vectors to matrices =============================================
  ka_sw_mat <- matrix(data = rep(ka_sw, each = limit + 1),
                      ncol = length(ka_sw),
                      nrow = limit + 1)
  ka_l_mat <- matrix(data = rep(ka_l, each = limit + 1),
                     ncol = length(ka_l),
                     nrow = limit + 1)
  ka_t_mat <- matrix(data = rep(ka_t, each = limit + 1),
                     ncol = length(ka_t),
                     nrow = limit + 1)
  # Calculate Legendre polynomial ==============================================
  Pl <- rep(c(1, -1), length(m)*2)[1:length(m)]
  # Calculate spherical Bessel functions of first kind =========================
  js_mat <- js(m, ka_sw_mat)
  js_mat_l <- js(m, ka_l_mat)
  js_mat_t <- js(m, ka_t_mat)
  # Calculate spherical Bessel functions of second kind ========================
  ys_mat <- ys(m, ka_sw_mat)
  # Calculate first derivative of spheric Bessel functions of first kind =======
  jsd_mat <- jsd(m, ka_sw_mat)
  jsd_mat_l <- jsd(m, ka_l_mat)
  jsd_mat_t <- jsd(m, ka_t_mat)
  # Calculate first derivative of spheric Bessel functions of second kind ======
  ysd_mat <- ysd(m, ka_sw_mat)
  # Calculate density contrast =================================================
  g <- model_body$density / medium$density
  # Tangent functions ==========================================================
  tan_sw <- -ka_sw_mat * jsd_mat / js_mat
  tan_l <- -ka_l_mat * jsd_mat_l / js_mat_l
  tan_t <- -ka_t_mat * jsd_mat_t / js_mat_t
  tan_beta <- -ka_sw_mat * ysd_mat / ys_mat
  tan_diff <- -js_mat / ys_mat
  # Difference terms ===========================================================
  along_m <- (m^2 + m)
  tan_l_add <- tan_l + 1
  tan_t_div <- along_m - 1 - ka_t_mat^2 / 2 + tan_t
  numerator <- (tan_l / tan_l_add) - (along_m / tan_t_div)
  denominator1 <- (along_m - ka_t_mat^2 / 2 + 2 * tan_l) / tan_l_add
  denominator2 <- along_m * (tan_t + 1) / tan_t_div
  denominator <- denominator1 - denominator2
  ratio <- -0.5 * ka_t_mat^2 * numerator / denominator
  # Additional trig functions ==================================================
  phi <- -ratio / g
  eta_tan <- tan_diff * (phi + tan_sw) / (phi + tan_beta)
  cos_eta <- 1 / sqrt(1 + eta_tan^2)
  sin_eta <- eta_tan * cos_eta
  # Fill in rest of Hickling (1962) equation ===================================
  f_j <- colSums((2 * m + 1) * Pl * (sin_eta * (1i * cos_eta - sin_eta)))
  # Calculate linear backscatter coefficient ===================================
  f_bs <- -2i * f_j / ka_sw
  sigma_bs <- (abs(f_bs) * model_body$radius / 2)^2
  TS <- 10 * log10(sigma_bs)
  slot(object, "model")$calibration <- data.frame(f_bs = f_bs,
                                                  sigma_bs = sigma_bs,
                                                  TS = TS)
  return(object)
}
################################################################################
################################################################################
# GAS-FILLED SCATTERERS
################################################################################
################################################################################
# Anderson (1950) gas-filled fluid-sphere model
################################################################################
#' Calculates the theoretical TS of a fluid sphere using an exact modal series
#' solution proposed by Andersen (1950).
#'
#' @param object GAS- or SBF-class object.
#' @usage
#' anderson_model(object)
#' @details
#' Calculates the theoretical TS of a fluid sphere at a given frequency using
#' an exact modal series solution.
#' @return
#' Target strength (TS, dB re: 1 m^2)
#' @references
#' Anderson, V.C. (1950). Sound scattering from a fluid sphere. Journal of the
#' Acoustical Society of America, 22, 426-431.
#' @export
anderson_model <- function(object) {
  # Extract model parameters/inputs ========================================
  body <- extract(object, "body")
  acoustics <- extract(object, "model_parameters")$anderson$acoustics
  # Multiple acoustic wavenumber by radius =================================
  k1a <- acoustics$k_sw * body$radius
  k2a <- acoustics$k_sw * body$radius
  # Set limit for iterations ===============================================
  ka_limit <- extract(object, "model_parameters")$anderson$modal$ka_limit
  m <- 0:ka_limit
  # Convert ka vectors to matrices =========================================
  ka_sw_mat <- matrix(data = rep(k1a, each = ka_limit + 1),
                      ncol = length(k1a),
                      nrow = ka_limit + 1)
  ka_b_mat <- matrix(data = rep(k2a, each = ka_limit + 1),
                     ncol = length(k2a),
                     nrow = ka_limit + 1)
  # Calculate modal series coefficient, b_m (or C_m) =======================
  gh <- body$g * body$h
  # Numerator term =========================================================
  num1 <- (jsd(m, ka_b_mat) * ys(m, ka_sw_mat)) / (js(m, ka_b_mat) * jsd(m, ka_sw_mat))
  num2 <- (gh * ysd(m, ka_sw_mat) / jsd(m, ka_sw_mat))
  c_num <- num1 - num2
  # Denominator term =======================================================
  dem1 <- (jsd(m, ka_b_mat) * js(m, ka_sw_mat)) / (js(m, ka_b_mat) * jsd(m, ka_sw_mat))
  dem2 <- gh
  c_dem <- dem1 - dem2
  # Now define coefficient =================================================
  C_m <- c_num / c_dem
  b_m <- -1 / (1 + 1i * C_m)
  # Sum across all columns to complete the modal series summation ==========
  f_j <- colSums((2 * m + 1) * (-1)^m * b_m)
  # Calculate the linear backscatter term for the gas-filled sphere
  f_sphere <- -1i / acoustics$k_sw * f_j
  # Return object ======================================================
  slot(object, "model")$fluid_sphere$anderson <- data.frame(f_sphere = f_sphere,
                                                            f_bs = abs(f_sphere),
                                                            sigma_bs = abs(f_sphere)^2,
                                                            TS = 20 * log10(abs(f_sphere)))
  return(object)
}
################################################################################
# Kirchoff-Ray Mode approximation
################################################################################
#' Calculates the theoretical TS using Kirchoff-ray Mode approximation.
#'
#' @param object Desired object/animal shape. Must be class "SBF".
#' @usage
#' KRM(object)
#' @details
#' Calculates the theoretical TS using the Kirchoff-ray Mode model.
#' @return
#' Target strength (TS, dB re: 1 m^2)
#' @references
#' Clay C.S. and Horne J.K. (1994). Acoustic models of fish: The Atlantic cod
#' (Gadus morhua). Journal of the Acoustical Society of AMerica, 96, 1661-1668.
#' @export
KRM <- function(object){
  # Extract model parameters/inputs ========================================
  body <- extract(object, "body")
  bladder <- extract(object, "bladder")
  model_params <- extract(object, "model_parameters")$KRM
  medium <- model_params$medium
  n_iteration <- model_params$parameters$parameters
  acoustics <- model_params$parameters$acoustics
  model_body <- model_params$scatterer$body
  model_bladder <- model_params$scatterer$bladder
  # Calculate reflection/transmission coefficients =========================
  RBC <- reflection_coefficient(body, bladder)
  RWB <- reflection_coefficient(medium, body)
  TT <- 1 - RWB^2
  # Calculate backscatter from swimbladder =================================
  # Sum across body/swimbladder position vectors ===========================
  bladder_rpos_sum <- along_sum(bladder$rpos, n_iteration$ncyl_sb)
  body_rpos_sum <- along_sum(body$rpos, n_iteration$ncyl_b)
  # Approximate radius of discrete cylinders ===============================
  a_bladder <- bladder_rpos_sum[2, ] / 4
  a_body <- body_rpos_sum[2, ] / 4
  # Combine wavenumber (k) and radii to calculate "ka" =====================
  # Generate appropriately sized matrix to be usable later =================
  ka_bladder <- matrix(data = rep(a_bladder, each = length(acoustics$k_sw)),
                       ncol = length(a_bladder),
                       nrow = length(acoustics$k_sw)) * acoustics$k_sw
  ka_body <- matrix(data = rep(a_body, each=length(acoustics$k_sw)),
                    ncol = length(a_body),
                    nrow = length(acoustics$k_b)) * acoustics$k_b
  # Calculate Kirchoff approximation empirical factor, A_sb ================
  A_sb <- ka_bladder / (ka_bladder + 0.083)
  # Calculate empirical phase shift for a fluid cylinder, Psi_p ============
  Psi_p <- ka_bladder / (40 + ka_bladder) - 1.05
  body_dorsal_sum <- matrix(data = rep(body_rpos_sum[3, ], each = length(acoustics$k_sw)),
                            ncol = length(body_rpos_sum[3, ]),
                            nrow = length(acoustics$k_sw)) / 2
  Psi_b <- -pi * acoustics$k_b * (body_dorsal_sum) /
    (2 * (acoustics$k_b * (body_dorsal_sum) + 0.4))
  # Convert x-z coordinates to requisite u-v rotated coordinates ===========
  uv_bladder <- bladder_rotation(bladder_rpos_sum,
                                 bladder$rpos,
                                 model_bladder$theta,
                                 length(acoustics$k_sw))
  uv_body <- body_rotation(body_rpos_sum,
                           body$rpos,
                           model_body$theta,
                           length(acoustics$k_sw))
  # Estimate natural log functions  ========================================
  exp_bladder <- exp(-1i * (2 * acoustics$k_b * uv_bladder$v + Psi_p)) *
    uv_bladder$delta_u
  exp_body <- exp(-2i * acoustics$k_sw * uv_body$vbU) * TT *
    exp(-2i * acoustics$k_sw * uv_body$vbU + 2i * acoustics$k_b *
          (uv_body$vbU - uv_body$vbL) + 1i * Psi_b)

  # Calculate the summation term ===========================================
  summation_soft <- A_sb * sqrt((ka_bladder + 1) * sin(model_bladder$theta))
  summation_fluid <- sqrt(ka_body) * uv_body$delta_u
  # Estimate backscattering length, f_fluid/f_soft =========================
  f_soft <- rowSums(-1i * (RBC * TT) / (2 * sqrt(pi)) * summation_soft * exp_bladder)
  f_fluid <- rowSums(-1i * (RWB / (2 * sqrt(pi))) * summation_fluid * exp_body)
  # Estimate total backscattering length, f_bs =============================
  f_bs <- f_soft + f_fluid
  # Return object ======================================================
  slot(object, "model")$KRM <- data.frame(f_fluid = f_fluid,
                                              f_soft = f_soft,
                                              f_bs = f_bs,
                                              TS = 20*log10(abs(f_bs)))
  return(object)
}
################################################################################
# Primary scattering model for an elastic shelled scatterers (ESS)
################################################################################
################################################################################
# Ray-based high pass approximation
################################################################################

#' Calculates the theoretical TS of a shelled organism using the non-modal
#' High Pass (HP) model
#'
#' @param object Desired animal object (Elastic Shelled).
#' @return
#' Target strength (TS, dB re: 1 m^2)
#' @references
#' Lavery, A.C., Wiebe, P.H., Stanton, T.K., Lawson, G.L., Benfield, M.C.,
#' Copley, N. 2007. Determining dominant scatterers of sound in mixed
#' zooplankton popuilations. The Journal of the Acoustical Society of America,
#' 122(6): 3304-3326.
#'
#' @export
stanton_high_pass <- function(object) {
  # Extract model parameters/inputs ============================================
  model_params <- extract(object, "model_parameters")$stanton_high_pass
  medium <- model_params$medium
  acoustics <- model_params$acoustics
  shell <- model_params$shell
  # Multiply acoustic wavenumber by body radius ================================
  k1a <- acoustics$k_sw * shell$radius
  # Calculate Reflection Coefficient ===========================================
  R <- R12(shell, medium)
  # Calculate backscatter constant, alpha_pi
  alpha_pi <- (1 - shell$g * shell$h^2) / (3 * shell$g * shell$h^2) +
    (1 - shell$g) / (1 + 2 * shell$g)
  # Define approximation constants, G_c and F_c ================================
  F_c <- 1; G_c <- 1
  # Caclulate numerator term ===================================================
  num <- (shell$radius^2 * k1a^4 * alpha_pi^2 * G_c)
  # Caclulate denominator term =================================================
  dem <- (1 + (4 * k1a^4 * alpha_pi^2) / (R^2 * F_c))
  # Calculate backscatter and return
  f_bs <- num / dem
  slot(object, "model")$stanton_high_pass <- data.frame(f_bs = f_bs,
                                                        sigma_bs = abs(f_bs),
                                                        TS = 10 * log10(abs(f_bs)))
  return(object)
}
