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
#' @return The theoretical acoustic target strength (TS, dB re. 1 \eqn{m^2}) of
#' a solid sphere at a given frequency.
#' @references
#' MacLennan D. N. (1981). The theory of solid spheres as sonar calibration
#' targets. Scottish Fisheries Research No. 22, Department of Agriculture and
#' Fisheries for Scotland.
#' @export
calibration <- function(object) {
  # Extract model parameters/inputs ============================================
  model <- acousticTS::extract(
    object,
    "model_parameters"
  )$calibration
  ### Now we solve / calculate equations =======================================
  # Equations 6a -- weight wavenumber by radius ================================
  ka_sw <- model$parameters$acoustics$k_sw * model$body$radius
  ka_l <- model$parameters$acoustics$k_l * model$body$radius
  ka_t <- model$parameters$acoustics$k_t * model$body$radius
  # Set limit for iterations ===================================================
  m_limit <- round(max(ka_sw)) + 10
  # Create modal series number vector ==========================================
  m <- 0:m_limit
  # Convert these vectors into matrices ========================================
  ka_sw_m <- modal_matrix(ka_sw, m_limit)
  ka_l_m <- modal_matrix(ka_l, m_limit)
  ka_t_m <- modal_matrix(ka_t, m_limit)
  # Calculate Legendre polynomial ==============================================
  Pl <- Pn(m, cos(model$body$theta))
  # Calculate spherical Bessel functions of first kind =========================
  js_mat <- js(m, ka_sw_m)
  js_mat_l <- js(m, ka_l_m)
  js_mat_t <- js(m, ka_t_m)
  # Calculate spherical Bessel functions of second kind ========================
  ys_mat <- ys(m, ka_sw_m)
  # Calculate first derivative of spheric Bessel functions of first kind =======
  jsd_mat <- jsd(m, ka_sw_m)
  jsd_mat_l <- jsd(m, ka_l_m)
  jsd_mat_t <- jsd(m, ka_t_m)
  # Calculate first derivative of spheric Bessel functions of second kind ======
  ysd_mat <- ysd(m, ka_sw_m)
  # Calculate density contrast =================================================
  g <- model$body$density / model$medium$density
  # Tangent functions ==========================================================
  tan_sw <- -ka_sw_m * jsd_mat / js_mat
  tan_l <- -ka_l_m * jsd_mat_l / js_mat_l
  tan_t <- -ka_t_m * jsd_mat_t / js_mat_t
  tan_beta <- -ka_sw_m * ysd_mat / ys_mat
  tan_diff <- -js_mat / ys_mat
  # Difference terms ===========================================================
  along_m <- (m * m + m)
  tan_l_add <- tan_l + 1
  tan_t_div <- along_m - 1 - ka_t_m * ka_t_m / 2 + tan_t
  numerator <- (tan_l / tan_l_add) - (along_m / tan_t_div)
  denominator1 <- (along_m - ka_t_m * ka_t_m / 2 + 2 * tan_l) / tan_l_add
  denominator2 <- along_m * (tan_t + 1) / tan_t_div
  denominator <- denominator1 - denominator2
  ratio <- -0.5 * (ka_t_m * ka_t_m) * numerator / denominator
  # Additional trig functions ==================================================
  phi <- -ratio / g
  eta_tan <- tan_diff * (phi + tan_sw) / (phi + tan_beta)
  cos_eta <- 1 / sqrt(1 + eta_tan * eta_tan)
  sin_eta <- eta_tan * cos_eta
  # Fill in rest of Hickling (1962) equation ===================================
  f_j <- colSums(
    (2 * m + 1) * Pl[m + 1] * (sin_eta * (1i * cos_eta - sin_eta))
  )
  # Calculate linear backscatter coefficient ===================================
  f_bs <- abs(-2i * f_j / ka_sw) * model$body$radius / 2
  sigma_bs <- f_bs * f_bs
  TS <- 10 * log10(sigma_bs)
  # Add results to scatterer object ============================================
  methods::slot(
    object,
    "model"
  )$calibration <- data.frame(
    frequency = model$parameters$acoustics$frequency,
    ka = ka_sw,
    f_bs = f_bs,
    sigma_bs = sigma_bs,
    TS = TS
  )
  # Return object ==============================================================
  object
}
