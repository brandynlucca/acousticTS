################################################################################
# Elastic cylinder modal series solution
################################################################################
#' Elastic cylinder modal series (ECMS) solution
#'
#' @description
#' Computes backscatter from a finite elastic cylinder by combining the
#' phase-shift solution for an infinite elastic cylinder with the finite-length
#' coherence factor used for normal and near-normal incidence.
#'
#' @section Usage:
#' This model is accessed via:
#' \preformatted{
#' target_strength(
#'   ...,
#'   model = "ECMS",
#'   sound_speed_sw,
#'   density_sw,
#'   density_body,
#'   sound_speed_longitudinal_body,
#'   sound_speed_transversal_body,
#'   m_limit
#' )
#' }
#'
#' @section Arguments:
#' \describe{
#'   \item{\code{sound_speed_sw}}{Seawater sound speed (\eqn{m~s^{-1}}).}
#'   \item{\code{density_sw}}{Seawater density (\eqn{kg~m^{-3}}).}
#'   \item{\code{density_body}}{Elastic-cylinder density (\eqn{kg~m^{-3}}).
#'   If omitted, the model uses the density stored on the scatterer body.}
#'   \item{\code{sound_speed_longitudinal_body}}{Longitudinal wave speed in the
#'   elastic cylinder (\eqn{m~s^{-1}}).}
#'   \item{\code{sound_speed_transversal_body}}{Transversal (shear) wave speed
#'   in the elastic cylinder (\eqn{m~s^{-1}}).}
#'   \item{\code{m_limit}}{Optional model truncation limit used to cap the
#'   number of retained cylindrical modes.}
#' }
#'
#' @references
#' Faran, J.J. (1951). Sound scattering by solid cylinders and spheres.
#' \emph{The Journal of the Acoustical Society of America}, 23: 405-418.
#'
#' Stanton, T.K. (1988). Sound scattering by cylinders of finite length. II.
#' Elastic cylinders. \emph{The Journal of the Acoustical Society of America},
#' 83: 64-67.
#'
#' @name ECMS
#' @aliases ecms ECMS
#' @docType data
#' @keywords models acoustics
NULL

#' Initialize scatterer-class object for the elastic cylinder modal series model.
#' @noRd
ecms_initialize <- function(object,
                            frequency,
                            sound_speed_sw = .SEAWATER_SOUND_SPEED_DEFAULT,
                            density_sw = .SEAWATER_DENSITY_DEFAULT,
                            density_body = NULL,
                            sound_speed_longitudinal_body = NULL,
                            sound_speed_transversal_body = NULL,
                            m_limit = NULL) {
  shape <- extract(object, "shape_parameters")

  if (!methods::is(object, "FLS")) {
    stop(
      "ECMS currently requires a fluid-like ('FLS') cylinder scatterer to ",
      "supply the geometry."
    )
  }

  if (shape$shape != "Cylinder") {
    stop(
      "The elastic cylinder modal series solution requires scatterer to be ",
      "shape-type 'Cylinder'. Input scatterer is shape-type '",
      shape$shape, "'."
    )
  }

  body <- extract(object, "body")
  rho_body <- density_body %||% body$density
  cL_body <- sound_speed_longitudinal_body %||% body$sound_speed_longitudinal
  cT_body <- sound_speed_transversal_body %||% body$sound_speed_transversal

  if (is.null(rho_body)) {
    stop(
      "ECMS requires 'density_body', either stored on the scatterer or passed ",
      "directly to target_strength()."
    )
  }
  if (is.null(cL_body) || is.null(cT_body)) {
    stop(
      "ECMS requires both 'sound_speed_longitudinal_body' and ",
      "'sound_speed_transversal_body'."
    )
  }

  if (abs(body$theta - pi / 2) > pi / 18) {
    warning(
      "ECMS is intended for broadside or near-broadside incidence."
    )
  }

  k_sw <- wavenumber(frequency, sound_speed_sw)
  k_l <- wavenumber(frequency, cL_body)
  k_t <- wavenumber(frequency, cT_body)
  ka_max <- pmax(k_sw, k_l, k_t) * shape$radius * abs(sin(body$theta))
  if (is.null(m_limit)) {
    m_limit <- ceiling(ka_max) + 10
  }

  methods::slot(object, "model_parameters")$ECMS <- list(
    parameters = list(
      acoustics = data.frame(
        frequency = frequency,
        k_sw = k_sw,
        k_l = k_l,
        k_t = k_t,
        m_limit = m_limit
      )
    ),
    medium = data.frame(
      sound_speed = sound_speed_sw,
      density = density_sw
    ),
    body = list(
      length = shape$length,
      radius = shape$radius,
      theta = body$theta,
      density = rho_body,
      sound_speed_longitudinal = cL_body,
      sound_speed_transversal = cT_body
    )
  )

  methods::slot(object, "model")$ECMS <- data.frame(
    frequency = frequency,
    sigma_bs = rep(NA_real_, length(frequency))
  )

  object
}

#' @noRd
.ecms_safe_divide <- function(numerator, denominator, eps = 1e-12) {
  denominator_adj <- ifelse(
    abs(denominator) < eps,
    ifelse(Re(denominator) < 0, -eps, eps),
    denominator
  )
  numerator / denominator_adj
}

#' @noRd
.ecms_backscatter_scalar <- function(k_sw,
                                     k_l,
                                     k_t,
                                     length_body,
                                     radius_body,
                                     theta_body,
                                     density_sw,
                                     density_body,
                                     m_limit) {
  sin_theta <- abs(sin(theta_body))
  if (sin_theta < 1e-10) {
    stop("ECMS is not implemented for end-on incidence.")
  }

  x1 <- k_l * radius_body * sin_theta
  x2 <- k_t * radius_body * sin_theta
  x3 <- k_sw * radius_body * sin_theta

  m <- 0:m_limit
  nu <- neumann(m)
  m2 <- m^2

  J1 <- jc(m, x1)
  J2 <- jc(m, x2)
  J3 <- jc(m, x3)
  Y3 <- yc(m, x3)

  tan_alpha1 <- -.ecms_safe_divide(x1 * jcd(m, x1), J1)
  tan_alpha2 <- -.ecms_safe_divide(x2 * jcd(m, x2), J2)
  tan_alpha3 <- -.ecms_safe_divide(x3 * jcd(m, x3), J3)
  tan_beta3  <- -.ecms_safe_divide(x3 * ycd(m, x3), Y3)
  tan_delta3 <- -.ecms_safe_divide(J3, Y3)

  X2 <- x2^2 / 2
  denom_a2 <- m2 - X2 + tan_alpha2

  term1 <- .ecms_safe_divide(tan_alpha1, tan_alpha1 + 1)
  term2 <- .ecms_safe_divide(m2, denom_a2)
  term3 <- .ecms_safe_divide(m2 - X2 + tan_alpha1, tan_alpha1 + 1)
  term4 <- .ecms_safe_divide(m2 * (tan_alpha2 + 1), denom_a2)

  tan_zeta <- .ecms_safe_divide(
    (-X2) * (term1 - term2),
    term3 - term4
  )
  tan_phi <- -(density_sw / density_body) * tan_zeta
  tan_eta <- tan_delta3 * .ecms_safe_divide(
    tan_phi + tan_alpha3,
    tan_phi + tan_beta3
  )

  cos_eta <- 1 / sqrt(1 + tan_eta^2)
  sin_eta <- tan_eta * cos_eta
  modal_sum <- sum(
    (-1)^m * nu * sin_eta * (cos_eta - 1i * sin_eta),
    na.rm = TRUE
  )

  A <- k_sw * length_body * cos(theta_body)
  sinc_factor <- if (abs(A) < 1e-10) 1 else sin(A) / A

  -(length_body / pi) * sinc_factor * modal_sum
}

#' Elastic cylinder modal series solution.
#' @noRd
ECMS <- function(object) {
  model <- extract(object, "model_parameters")$ECMS
  acoustics <- model$parameters$acoustics
  body <- model$body
  medium <- model$medium

  f_bs <- mapply(
    FUN = .ecms_backscatter_scalar,
    k_sw = acoustics$k_sw,
    k_l = acoustics$k_l,
    k_t = acoustics$k_t,
    m_limit = acoustics$m_limit,
    MoreArgs = list(
      length_body = body$length,
      radius_body = body$radius,
      theta_body = body$theta,
      density_sw = medium$density,
      density_body = body$density
    )
  )

  sigma_bs <- abs(f_bs)^2

  methods::slot(object, "model")$ECMS <- data.frame(
    frequency = acoustics$frequency,
    ka = acoustics$k_sw * body$radius * abs(sin(body$theta)),
    f_bs = f_bs,
    sigma_bs = sigma_bs,
    TS = 10 * log10(sigma_bs)
  )

  object
}
