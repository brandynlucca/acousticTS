################################################################################
# Two-ray cylinder model (TRCM) for elongated scatterers
################################################################################
#' Two-ray cylinder model (TRCM) for elongated scatterers
#'
#' @description
#' Computes the far-field scattering amplitude and related quantities for
#' elongated fluid-like scatterers using the two-ray approximation model, as
#' described by Stanton et al. (1993, 1998). The two-ray model is a
#' high-frequency ray-based approximation that accounts for interference
#' between reflections from the front and back interfaces of the scatterer,
#' along with a directivity pattern that depends on the scatterer's orientation
#' and curvature. The model is computationally efficient and captures the
#' essential physics of acoustic scattering from elongated bodies such as
#' zooplankton, particularly at high frequencies where the acoustic wavelength
#' is much smaller than the organism size. The model supports both straight
#' and bent cylinders, with the bent cylinder formulation incorporating the
#' radius of curvature to account for body shape effects on the directivity
#' pattern. For more details, see the [expanded documentation on the two-ray cylinder
#' model](https://brandynlucca.github.io/acousticTS/articles/trcm/trcm-theory.html).
#'
#' @section Usage:
#' This model is accessed via:
#' \preformatted{
#' target_strength(
#'   ...,
#'   model = "TRCM",
#'   radius_curvature,
#'   radius_curvature_ratio,
#'   radius_cylinder_fun,
#'   sound_speed_sw,
#'   density_sw
#' )
#' }
#'
#' @section Arguments:
#' \describe{
#'   \item{\code{radius_curvature}}{Radius of curvature for bent cylinders
#'   (\eqn{m}). If \code{NULL}, the model assumes a straight cylinder unless
#'   \code{radius_curvature_ratio} is specified.}
#'   \item{\code{radius_curvature_ratio}}{Ratio of radius of curvature to body
#'   length (\eqn{\rho_c / L}). Used to compute \code{radius_curvature} if not
#'   explicitly provided. Default is \code{NULL}.}
#'   \item{\code{radius_cylinder_fun}}{Method for selecting the representative
#'   radius when the cylinder has variable radius. One of \code{"center"}
#'   (default), \code{"mean"}, \code{"median"}, or \code{"max"}.}
#'   \item{\code{sound_speed_sw}}{Seawater sound speed (\eqn{m~s^{-1}}).}
#'   \item{\code{density_sw}}{Seawater density (\eqn{kg~m^{-3}}).}
#' }
#'
#' @seealso
#' \code{\link{target_strength}}, \code{\link{FLS}}, \code{\link{Cylinder}},
#' \code{\link{cylinder}}
#'
#' @references
#' Stanton, T.K., Chu, D., Wiebe, P.H., and Clay, C.S. (1993a). Average
#' echoes from randomly oriented random-length finite cylinders: Zooplankton
#' models. The Journal of the Acoustical Society of America, 94: 3463-3472.
#'
#' Stanton, T.K., Chu, D., Wiebe, P.H., Martin, L.V., and Eastwood, R.L.
#' (1998). Sound scattering by several zooplankton groups. I. Experimental
#' determination of dominant scattering mechanisms. The Journal of the
#' Acoustical Society of America, 103: 225-235.
#'
#' Stanton, T.K., Chu, D., and Wiebe, P.H. (1998). Sound scattering by several
#' zooplankton groups. II. Scattering models. The Journal of the Acoustical
#' Society of America, 103: 236-253.
#'
#' @name TRCM
#' @aliases trcm TRCM
#' @docType data
#' @keywords models acoustics internal
NULL

#' Initialize Scatterer-class object for the two-ray model
#'
#' @param object FLS-class object.
#' @param frequency Transmit frequency (Hz).
#' @param sound_speed_sw Seawater sound speed (m/s).
#' @param density_sw Seawater density (kg/m³).
#' @noRd
trcm_initialize <- function(object,
                            frequency,
                            sound_speed_sw = 1500,
                            density_sw = 1026,
                            stationary_phase = FALSE) {
  # Define medium parameters early for contrast calculations ===================
  medium_params <- .init_medium_params(sound_speed_sw, density_sw)
  # Parse shape ================================================================
  shape_params <- acousticTS::extract(object, "shape_parameters")
  # Validate scatterer and shape types =========================================
  if (!methods::is(object, "FLS")) {
    stop(
      "The two-ray model requires the scatterer to be fluid-like ('FLS'). ",
      "Input scatterer is type '", class(object), "'."
    )
  }
  if (shape_params$shape != "Cylinder") {
    stop(
      "The two-ray model requires scatterer to be shape-type 'Cylinder'. ",
      "Input scatterer is shape-type '", shape_params$shape, "'."
    )
  }
  # Parse body =================================================================
  body <- acousticTS::extract(object, "body")
  body_params <- list(
    length = shape_params$length,
    radius = shape_params$radius,
    theta = body$theta,
    is_bent = FALSE,
    theta_shift = body$theta - pi / 2
  )
  # Determine if cylinder is bent or straight ==================================
  is_bent <- !is.null(shape_params$radius_curvature_ratio) &&
    !is.na(shape_params$radius_curvature_ratio)
  # Extract radius of curvature if present =====================================
  if (is_bent) {
    # Interpret stored ratio as curvature radius scaled by body length +++++++++
    body_params$radius_curvature <- shape_params$radius_curvature_ratio *
      shape_params$length
    body_params$is_bent <- TRUE
  } else {
    body_params$radius_curvature <- NA
  }
  # Fill out material properties ===============================================
  body_params$h <- body$h %||%
    if (!is.null(body$sound_speed)) body$sound_speed / medium_params$sound_speed else NA
  body_params$g <- body$g %||%
    if (!is.null(body$density)) body$density / medium_params$density else NA
  # Define model parameters ====================================================
  model_params <- list(
    acoustics = transform(
      .init_acoustics_df(
        frequency,
        k_sw = sound_speed_sw,
        k_b = sound_speed_sw * body_params$h
      ),
      lambda = sound_speed_sw / frequency
    ),
    stationary_phase = stationary_phase
  )
  .init_model_slots(
    object = object,
    model_name = "TRCM",
    frequency = frequency,
    model_parameters = list(
      parameters = model_params,
      medium = medium_params,
      body = body_params
    )
  )
}

#' Helper function for handling the straight-cylinder form for TRCM.
#' @keywords internal
#' @noRd
.trcm_straight <- function(k1, k2a, l, a, r, I, theta_shift) {
  # Calculate wavenumber-length and -radius products ===========================
  k1l <- k1 * l
  k1a <- k1 * a
  # Update the interference term ===============================================
  I_term <- 1 - I * exp(4i * k2a * cos(theta_shift))
  # Evaluate sinc function directivity =========================================
  Delta <- k1 * l * sin(theta_shift)
  sinc_directivity <- ifelse(abs(Delta) < 1e-10, 1, sin(Delta) / Delta)
  # Compute the linear scattering length, f_bs =================================
  (-1i / (2 * sqrt(pi))) * exp(1i * pi / 4) *
    exp(-2i * k1 * a * cos(theta_shift)) * l * sqrt(k1 * a * cos(theta_shift)) *
    r * sinc_directivity * I_term
}

#' Helper function for handling the Fresnel integration for equivalent length
#' @keywords internal
#' @noRd
.trcm_equivalent_length_fresnel <- function(k1, l, a, rho_c) {
  # Get the maximum curvature ==================================================
  gamma_max <- l / (2 * rho_c)
  z_max <- rho_c * (1 - cos(gamma_max))
  A <- z_max / a
  B <- l / a
  # Integrate over curvature to get the equivalent length ======================
  integrand <- function(x, k1, A, B, a) {
    exp(1i * 8 * k1 * A * x^2 / (B^2 * a))
  }
  # Integration is done over all k1 ++++++++++++++++++++++++++++++++++++++++++++
  lebc_scalar <- function(k1, A, B, a, lower, upper) {
    real_part <- stats::integrate(
      function(x) Re(integrand(x, k1, A, B, a)),
      lower,
      upper
    )$value
    imag_part <- stats::integrate(
      function(x) Im(integrand(x, k1, A, B, a)),
      lower,
      upper
    )$value
    real_part + 1i * imag_part
  }
  vapply(
    k1,
    lebc_scalar,
    complex(1),
    A = A,
    B = B,
    a = a,
    lower = -l / 2,
    upper = l / 2
  )
}

#' Helper function for handling the bent cylinder stationary phase
#' @keywords internal
#' @noRd
.trcm_equivalent_length_stationary_phase <- function(lambda, l, rho_c) {
  sqrt(rho_c * lambda / 2) * exp(1i * pi / 4)
}

#' Helper function for handling the bent-cylinder form for TRCM.
#' @keywords internal
#' @noRd
.trcm_curved <- function(
    k1, k2a, lambda, l, a, rho_c, r, I, theta_shift, stationary_phase
  ) {
  # Stationary phase correction ================================================
  if (stationary_phase) {
    # Simplified equivalent bent cylinder length +++++++++++++++++++++++++++++++
    return(
      .trcm_equivalent_length_stationary_phase(lambda, l, rho_c) *
        .trcm_straight(k1, k2a, l, a, r, I, 0) / l
    )
  }
  lebc <- .trcm_equivalent_length_fresnel(k1, l, a, rho_c)
  # Compute the linear scattering strength, f_bs ===============================
  lebc * .trcm_straight(k1, k2a, l, a, r, I, theta_shift) / l
}

#' Calculates the theoretical TS of a fluid-like straight or bent cylinder at a
#' given frequency using the two-ray path finite cylinder model.
#'
#' @param object FLS-class scatterer.
#' @noRd
TRCM <- function(object) {
  # Extract model parameters ===================================================
  model <- acousticTS::extract(object, "model_parameters")$TRCM
  acoustics <- model$parameters$acoustics
  body <- model$body
  # Calculate ka and k2a =======================================================
  k1a <- acoustics$k_sw * body$radius
  k2a <- acoustics$k_b * body$radius
  # Calculate the reflection coefficient, R ====================================
  r <- (body$g * body$h - 1) / (body$g * body$h + 1)
  # Calculate transmission coefficients ========================================
  T12 <- 1 - r
  T21 <- 1 + r
  # Calculate phase advancement term ===========================================
  mu <- (-pi / 2 * k1a) / (k1a + 0.4)
  # Calculate interference term ================================================
  I_term <- T12 * T21 * exp(1i * mu)
  # Compute the linear scattering length, f_bs =================================
  if (body$is_bent) {
    fbs <- .trcm_curved(
      acoustics$k_sw, k2a, acoustics$lambda, body$length, body$radius,
      body$radius_curvature, r, I_term, body$theta_shift,
      model$parameters$stationary_phase
    )
  } else {
    fbs <- .trcm_straight(
      acoustics$k_sw, k2a, body$length, body$radius, r, I_term, body$theta_shift
    )
  }
  # Update scatterer object ====================================================
  methods::slot(object, "model")$TRCM <- data.frame(
    frequency = acoustics$frequency,
    ka = k1a,
    f_bs = fbs,
    sigma_bs = abs(fbs)^2,
    TS = 10 * log10(abs(fbs)^2)
  )
  # Return object ==============================================================
  object
}
