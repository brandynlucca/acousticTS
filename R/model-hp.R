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
#' Stanton, T.K. (1989). Simple approximate formulas for backscattering of
#' sound by spherical and elongated objects. The Journal of the Acoustical
#' Society of America, 86, 1499-1510.
#'
#' @export
high_pass_stanton <- function(object) {
  # Extract model parameters/inputs ==========================================
  model_params <- extract(object, "model_parameters")$high_pass_stanton
  medium <- model_params$medium
  acoustics <- model_params$acoustics
  shell <- model_params$shell
  # Multiply acoustic wavenumber by body radius ==============================
  k1a <- acoustics$k_sw * shell$radius
  # Calculate Reflection Coefficient =========================================
  R <- (shell$h * shell$g - 1) / (shell$h * shell$g + 1)
  # Calculate backscatter constant, alpha_pi
  alpha_pi <- (1 - shell$g * (shell$h * shell$h)) /
    (3 * shell$g * (shell$h * shell$h)) +
    (1 - shell$g) / (1 + 2 * shell$g)
  # Define approximation constants, G_c and F_c ==============================
  F_c <- 1
  G_c <- 1
  # Caclulate numerator term =================================================
  num <- ((shell$radius * shell$radius) * (k1a * k1a * k1a * k1a) *
            (alpha_pi * alpha_pi) * G_c)
  # Caclulate denominator term ===============================================
  dem <- (1 + (4 * (k1a * k1a * k1a * k1a) * (alpha_pi * alpha_pi)) /
            ((R * R) * F_c))
  # Calculate backscatter and return
  f_bs <- num / dem
  methods::slot(object, "model")$high_pass_stanton <- data.frame(
    frequency = acoustics$frequency,
    k1a = acoustics$k_sw * shell$radius,
    k_s = acoustics$k_b,
    f_bs = f_bs,
    sigma_bs = abs(f_bs),
    TS = 10 * log10(abs(f_bs))
  )
  object
}
