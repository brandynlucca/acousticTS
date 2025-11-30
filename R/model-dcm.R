################################################################################
# Primary scattering models for fluid-like scatterers (FLS)
################################################################################
#' Calculates the theoretical TS of a target using the deformed cylinder
#' model (DCM).
#'
#' @param object FLS-class scatterer.
#' @usage
#' DCM( object )
#' @return
#' Target strength (TS, dB re. 1 \eqn{m^2}) of a FLS-class object.
#' #' @references
#' Stanton, T.K., Chu, D., and Wiebe, P.H. 1998. Sound scattering by several
#' zooplankton groups. II. Scattering models. Journal of the Acoustical Society
#' of America, 103(1), 236-253.
#' @export
DCM <- function(object) {
  # Extract model parameters/inputs ============================================
  model <- acousticTS::extract(object, "model_parameters")$DCM
  body <- acousticTS::extract(object, "body")
  body$density <- model$medium$density * body$g
  body$sound_speed <- model$medium$sound_speed * body$h
  # Multiply acoustic wavenumber by body radius ================================
  k1a <- model$parameters$acoustics$k_sw * model$body$radius
  k2a <- model$parameters$acoustics$k_b * model$ body$radius
  # Calculate Reflection Coefficient ===========================================
  R12 <- acousticTS::reflection_coefficient(model$medium, body)
  # Calculate Transmission Coefficients ========================================
  T12 <- 1 - R12
  T21 <- 1 + R12
  # Calculate phase shift term for the ray model ===============================
  mu <- (-pi / 2 * k1a) / (k1a + 0.4)
  # Calculate Interference Term between echoes from front/back interfaces ======
  I0 <- 1 - T12 * T21 * base::exp(4i * k2a) * base::exp(1i * mu * k1a)
  # Calculate half-width of directivity pattern ================================
  w12 <- base::sqrt(model$body$radius_curvature * model$body$radius)
  # Apply shift to orientation required for this model =========================
  theta_shift <- body$theta - pi / 2
  # Calculate Taylor expansion of directivity function =========================
  D_theta <- base::exp(-model$body$alpha_B *
                         (2 * theta_shift * model$body$radius_curvature /
                            model$body$length)^2)
  # Calculate linear backscatter function, f_bs ================================
  f_bs <- 0.5 * w12 * R12 * base::exp(-2i * k1a) * I0 * D_theta
  # Update scatterer object ====================================================
  methods::slot(object, "model")$DCM <- base::data.frame(
    frequency = model$parameters$acoustics$frequency,
    ka = k1a,
    f_bs = f_bs,
    sigma_bs = base::abs(f_bs)^2,
    TS = 10 * base::log10(base::abs(f_bs)^2)
  )
  # Return object ==============================================================
  object
}
