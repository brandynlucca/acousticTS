################################################################################
# Primary scattering model for an elastic solid sphere (CAL)
################################################################################
#' Solid elastic (calibration) sphere modal series (SOEMS) solution
#'
#' Calculates the far-field scattering amplitude and related quantities for a
#' solid elastic (calibration) sphere using a modal series solution.
#'
#' @section Usage:
#' This model is accessed via:
#' \preformatted{
#' target_strength(
#'   ...,
#'   model="calibration",
#'   sound_speed_sw,
#'   density_sw,
#'   adaptive = TRUE
#' )
#' }
#'
#' @section Arguments:
#' \describe{
#'    \item{\code{sound_speed_sw}}{Seawater sound speed (\eqn{m~s^{-1}}).}
#'    \item{\code{density_sw}}{Seawater density (\eqn{kg~m^{-3}}).}
#'    \item{\code{adaptive}}{Logical. If \code{TRUE}, extend the partial-wave
#'    sum beyond the initial \eqn{\mathrm{round}(ka)+10} modal cap until the
#'    tail term falls below the internal convergence threshold. If
#'    \code{FALSE}, use the original fixed modal cutoff only.}
#' }
#'
#' @section Theory:
#'
#' The calibration sphere model computes acoustic scattering from a solid
#' elastic sphere by expanding the incident and scattered fields in terms of
#' spherical Bessel and Hankel functions and Legendre polynomials. Both
#' compressional and shear waves within the sphere are included, and the
#' appropriate boundary conditions are enforced at the sphere-water interface.
#'
#' The dimensionless frequency parameter is defined as \eqn{q = ka}, where
#' \eqn{k} is the wavenumber in water and \eqn{a} is the sphere radius. The
#' longitudinal and transverse wave numbers inside the sphere are
#'
#'  \deqn{
#'    q_1 = \frac{qc}{c_1} \\
#'    q_2 = \frac{qc}{c_2},
#'  }
#'
#' where \eqn{c_1} and \eqn{c_2} are the longitudinal and transverse sound
#' speeds in the sphere, respectively.
#'
#' The far-field backscattering form function is given by:
#'
#'  \deqn{
#'    f_\infty(q) = -\frac{2}{q} \sum\limits_{\ell=0}^{\infty}
#'    (-1)^\ell (2\ell+1)
#'    \sin \eta_\ell \exp(i \eta_\ell)
#'  }
#'
#' The phase angle for each mode is then given by:
#'
#'  \deqn{
#'    \tan \eta_\ell = -\frac{B_2 j_\ell'(q) - B_1 j_\ell(q)}{B_2 y_\ell'(q) -
#'    B_1 y_\ell(q)}
#'  }
#'
#' where the phase angle \eqn{\eta_\ell} is determined by the boundary
#' conditions and the material properties of the sphere and surrounding fluid.
#' The phase angle \eqn{\eta_\ell} is calculated using a series of
#' coefficients. The auxilliary quantities used for determining the boundary
#' conditions and subsequently \eqn{\eta_\ell} are reported in MacLennan (1981).
#'
#' @seealso
#' \code{\link{target_strength}}, \code{\link{CAL}}, \code{\link{Sphere}},
#' \code{\link{sphere}}
#'
#' @references
#'
#' Hickling, R. (1962). Analysis of echoes from a solid elastic sphere in
#' water. The Journal of the Acoustical Society of America, 34: 1582-1592.
#'
#' MacLennan D. N. (1981). The theory of solid spheres as sonar calibration
#' targets. Scottish Fisheries Research No. 22, Department of Agriculture and
#' Fisheries for Scotland.
#'
#'
#' @name SOEMS
#' @aliases soems SOEMS calibration CALIBRATION
#' @docType data
#' @keywords models acoustics
NULL

#' Initialize CAL-class object for modeling.
#' @param object CAL-class object.
#' @param frequency Frequency vector (Hz).
#' @param sound_speed_sw Seawater sound speed (m/s).
#' @param density_sw Seawater density (kg/m³).
#' @noRd
calibration_initialize <- function(
    object,
    frequency,
    sound_speed_sw = .SEAWATER_SOUND_SPEED_DEFAULT,
    density_sw = .SEAWATER_DENSITY_DEFAULT,
    adaptive = TRUE
  ) {
  # Parse shape ================================================================
  shape <- acousticTS::extract(
    object,
    "shape_parameters"
  )
  # Parse body =================================================================
  body <- acousticTS::extract(
    object,
    "body"
  )
  # Define medium parameters ===================================================
  medium_params <- data.frame(
    sound_speed = sound_speed_sw,
    density = density_sw
  )
  # Define model parameters recipe =============================================
  model_params <- list(
    acoustics = data.frame(
      frequency = frequency,
      # Wavenumber (medium) ====================================================
      k_sw = acousticTS::wavenumber(
        frequency,
        sound_speed_sw
      ),
      # Wavenumber (longitudinal axis) =========================================
      k_l = acousticTS::wavenumber(
        frequency,
        body$sound_speed_longitudinal
      ),
      # Wavenumber (tranversal axis) ===========================================
      k_t = acousticTS::wavenumber(
        frequency,
        body$sound_speed_transversal
      )
    ),
    parameters = data.frame(
      ncyl_b = shape$n_segments,
      adaptive = adaptive
    )
  )
  # Define body parameters recipe ==============================================
  body_params <- data.frame(
    diameter = body$diameter,
    radius = body$radius,
    sound_speed_longitudinal = body$sound_speed_longitudinal,
    sound_speed_transversal = body$sound_speed_transversal,
    density = body$density,
    theta = body$theta
  )
  # Add model parameters slot to scattering object =============================
  methods::slot(
    object,
    "model_parameters"
  )$calibration <- list(
    parameters = model_params,
    medium = medium_params,
    body = body_params
  )
  # Add model results slot to scattering object ================================
  methods::slot(
    object,
    "model"
  )$calibration <- data.frame(
    frequency = frequency,
    sigma_bs = rep(
      NA,
      length(frequency)
    )
  )
  # Output =====================================================================
  return(object)
}

#' Calculates theoretical TS of a solid sphere of a certain material at a given
#' frequency.
#'
#' @param object CAL-class object.
#' @noRd
calibration <- function(object) {
  # Extract model parameters/inputs ============================================
  model <- acousticTS::extract(
    object,
    "model_parameters"
  )$calibration
  # Convergence settings =======================================================
  modal_tol <- 1e-10
  max_modal_order <- 500L
  adaptive <- isTRUE(model$parameters$parameters$adaptive[[1]])
  ### Now we solve / calculate equations =======================================
  # Equations 6a -- weight wavenumber by radius ================================
  ka_sw <- model$parameters$acoustics$k_sw * model$body$radius
  ka_l <- model$parameters$acoustics$k_l * model$body$radius
  ka_t <- model$parameters$acoustics$k_t * model$body$radius
  # Helper for a single partial-wave term ======================================
  modal_term <- function(m,
                         ka_sw,
                         ka_l,
                         ka_t,
                         theta,
                         density_body,
                         density_sw) {
    # Calculate Legendre polynomial ============================================
    Pl <- as.numeric(Pn(m, cos(theta)))
    # Calculate spherical Bessel functions of first kind =======================
    js_mat <- js(m, ka_sw)
    js_mat_l <- js(m, ka_l)
    js_mat_t <- js(m, ka_t)
    # Calculate spherical Bessel functions of second kind ======================
    ys_mat <- ys(m, ka_sw)
    # Calculate first derivative of spherical Bessel functions of first kind ===
    jsd_mat <- jsd(m, ka_sw)
    jsd_mat_l <- jsd(m, ka_l)
    jsd_mat_t <- jsd(m, ka_t)
    # Calculate first derivative of spherical Bessel functions of second kind ==
    ysd_mat <- ysd(m, ka_sw)
    # Calculate density contrast ===============================================
    g <- density_body / density_sw
    # Tangent functions ========================================================
    tan_sw <- -ka_sw * jsd_mat / js_mat
    tan_l <- -ka_l * jsd_mat_l / js_mat_l
    tan_t <- -ka_t * jsd_mat_t / js_mat_t
    tan_beta <- -ka_sw * ysd_mat / ys_mat
    tan_diff <- -js_mat / ys_mat
    # Difference terms =========================================================
    along_m <- (m * m + m)
    tan_l_add <- tan_l + 1
    tan_t_div <- along_m - 1 - ka_t * ka_t / 2 + tan_t
    numerator <- (tan_l / tan_l_add) - (along_m / tan_t_div)
    denominator1 <- (along_m - ka_t * ka_t / 2 + 2 * tan_l) / tan_l_add
    denominator2 <- along_m * (tan_t + 1) / tan_t_div
    denominator <- denominator1 - denominator2
    ratio <- -0.5 * (ka_t * ka_t) * numerator / denominator
    # Additional trig functions ===============================================
    phi <- -ratio / g
    eta_tan <- tan_diff * (phi + tan_sw) / (phi + tan_beta)
    cos_eta <- 1 / sqrt(1 + eta_tan * eta_tan)
    sin_eta <- eta_tan * cos_eta
    # Hickling (1962) / MacLennan (1981) modal term ===========================
    (2 * m + 1) * Pl * (sin_eta * (1i * cos_eta - sin_eta))
  }
  # Compute form function ======================================================
  f_j <- mapply(FUN = function(ka_sw,
                               ka_l,
                               ka_t,
                               theta,
                               density_body,
                               density_sw) {
      # Start with the common MacLennan / echoSMs truncation guess =============
      ml <- as.integer(round(ka_sw) + 10L)
      m <- 0:ml
      terms <- vapply(
        X = m,
        FUN = modal_term,
        FUN.VALUE = complex(1L),
        ka_sw = ka_sw,
        ka_l = ka_l,
        ka_t = ka_t,
        theta = theta,
        density_body = density_body,
        density_sw = density_sw
      )
      # Extend the modal sum until the tail term is negligible =================
      if (adaptive) {
        while (Mod(tail(terms, 1L)) > modal_tol && ml < max_modal_order) {
          ml <- ml + 1L
          terms <- c(
            terms,
            modal_term(
              m = ml,
              ka_sw = ka_sw,
              ka_l = ka_l,
              ka_t = ka_t,
              theta = theta,
              density_body = density_body,
              density_sw = density_sw
            )
          )
        }
      }
      sum(terms)
    }, ka_sw, ka_l, ka_t, model$body$theta, model$body$density,
    model$medium$density)
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
