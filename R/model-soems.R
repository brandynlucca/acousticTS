################################################################################
# Primary scattering model for an elastic solid sphere (CAL) ===================
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
#' @return No value; this help topic documents the calibration model.
#' @examples
#' subset(available_models(), model == "calibration")
#'
#' @name SOEMS
#' @aliases soems SOEMS calibration CALIBRATION
#' @keywords models acoustics internal
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
  # Define model parameters recipe =============================================
  model_params <- list(
    acoustics = .init_acoustics_df(
      frequency,
      k_sw = sound_speed_sw,
      k_l = body$sound_speed_longitudinal,
      k_t = body$sound_speed_transversal
    ),
    parameters = .model_data_frame(
      ncyl_b = shape$n_segments,
      adaptive = adaptive
    )
  )
  # Define body parameters recipe ==============================================
  body_params <- .model_data_frame(
    diameter = body$diameter,
    radius = body$radius,
    sound_speed_longitudinal = body$sound_speed_longitudinal,
    sound_speed_transversal = body$sound_speed_transversal,
    density = body$density,
    theta = body$theta
  )
  .init_model_slots(
    object = object,
    model_name = "calibration",
    frequency = frequency,
    model_parameters = list(
      parameters = model_params,
      medium = .init_medium_params(sound_speed_sw, density_sw),
      body = body_params
    )
  )
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
  adaptive <- isTRUE(model$parameters$parameters$adaptive[[1]])
  ### Now we solve / calculate equations =======================================
  # Equations 6a -- weight wavenumber by radius ================================
  ka_sw <- model$parameters$acoustics$k_sw * model$body$radius
  ka_l <- model$parameters$acoustics$k_l * model$body$radius
  ka_t <- model$parameters$acoustics$k_t * model$body$radius
  f_j <- calibration_modal_sum_cpp(
    ka_sw = ka_sw,
    ka_l = ka_l,
    ka_t = ka_t,
    theta = model$body$theta,
    density_ratio = model$body$density / model$medium$density,
    adaptive = adaptive
  )
  # Calculate linear backscatter coefficient ===================================
  f_bs <- abs(-2i * f_j / ka_sw) * model$body$radius / 2
  sigma_bs <- f_bs * f_bs
  TS <- 10 * log10(sigma_bs)
  # Add results to scatterer object ============================================
  methods::slot(
    object,
    "model"
  )$calibration <- .model_data_frame(
    frequency = model$parameters$acoustics$frequency,
    ka = ka_sw,
    f_bs = f_bs,
    sigma_bs = sigma_bs,
    TS = TS
  )
  # Return object ==============================================================
  object
}
