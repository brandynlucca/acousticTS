################################################################################
# Elastic cylinder modal series solution
################################################################################
#' Elastic cylinder modal series (ECMS) solution
#'
#' @description
#' Computes backscatter from straight and uniformly bent finite elastic
#' cylinders by combining the phase-shift solution for an infinite elastic
#' cylinder with the finite-length coherence factor used for normal and
#' near-normal incidence, plus the Fresnel coherence correction used for
#' uniformly bent cylinders. In the current package class system, the elastic
#' cylinder is carried most naturally by an \code{ESS}-class cylinder within
#' the broader elastic-based \code{ELA} family, with the elastic-cylinder
#' properties stored on the shell slot. Legacy \code{FLS}-class cylinders are
#' also accepted as geometry carriers for backward compatibility.
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
#' Stanton, T.K. (1989). Sound scattering by cylinders of finite length. III.
#' Deformed cylinders. \emph{The Journal of the Acoustical Society of America},
#' 86: 691-705.
#'
#' Gorska, N., Ona, E., and Korneliussen, R. (2005). Acoustic backscattering by
#' Atlantic mackerel as being representative of fish that lack a swimbladder.
#' Backscattering by individual fish. \emph{ICES Journal of Marine Science},
#' 62: 984-995.
#'
#' @name ECMS
#' @aliases ecms ECMS
#' @docType data
#' @keywords models acoustics internal
NULL

#' Initialize scatterer-class object for the elastic cylinder modal series
#' model.
#' @noRd
.ecms_validate_scope <- function(object, shape) {
  # Limit ECMS to cylindrical ESS and legacy cylindrical FLS objects ===========
  if (!(methods::is(object, "ESS") || methods::is(object, "FLS"))) {
    stop(
      "ECMS requires a cylindrical 'ESS' or legacy cylindrical 'FLS' scatterer."
    )
  }

  # Reject non-cylindrical geometries ==========================================
  if (shape$shape != "Cylinder") {
    stop(
      "The elastic cylinder modal series solution requires scatterer to be ",
      "shape-type 'Cylinder'. Input scatterer is shape-type '",
      shape$shape, "'."
    )
  }
}

#' Resolve the active elastic component for ECMS
#' @noRd
.ecms_resolve_component <- function(object) {
  # ESS objects store the elastic layer in the shell slot ======================
  if (methods::is(object, "ESS")) {
    return(extract(object, "shell"))
  }

  # Legacy cylindrical FLS objects store the elastic body directly ============
  extract(object, "body")
}

#' Resolve the ECMS body material properties
#' @noRd
.ecms_resolve_materials <- function(object,
                                    sound_speed_sw,
                                    density_sw,
                                    density_body = NULL,
                                    sound_speed_longitudinal_body = NULL,
                                    sound_speed_transversal_body = NULL) {
  # Recover the active elastic component and hydrate any stored contrasts ======
  elastic_comp <- .ecms_resolve_component(object)
  rho_body <- density_body %||% elastic_comp$density
  cL_body <- sound_speed_longitudinal_body %||%
    elastic_comp$sound_speed_longitudinal
  cT_body <- sound_speed_transversal_body %||%
    elastic_comp$sound_speed_transversal

  # Validate the required elastic material inputs ==============================
  if (is.null(rho_body)) {
    stop(
      "ECMS requires the elastic-cylinder density, either stored on the ",
      "scatterer or passed directly as 'density_body'."
    )
  }
  if (is.null(cL_body) || is.null(cT_body)) {
    stop(
      "ECMS requires both the longitudinal and transversal wave speeds, ",
      "either stored on the scatterer or passed directly to target_strength()."
    )
  }

  list(
    component = elastic_comp,
    density = rho_body,
    sound_speed_longitudinal = cL_body,
    sound_speed_transversal = cT_body
  )
}

#' Warn when ECMS inputs leave the intended incidence regime
#' @noRd
.ecms_warn_incidence_regime <- function(shape_core, elastic_comp) {
  # Warn when incidence departs too far from broadside =========================
  if (abs(elastic_comp$theta - pi / 2) > pi / 18) {
    warning("ECMS is intended for broadside or near-broadside incidence.")
  }

  # Bent-cylinder corrections share the same incidence limitation =============
  if (!is.null(shape_core$radius_curvature_ratio) &&
      !is.na(shape_core$radius_curvature_ratio) &&
      abs(elastic_comp$theta - pi / 2) > pi / 18) {
    warning(
      "ECMS bent-cylinder correction is intended for broadside or near-",
      "broadside incidence."
    )
  }
}

#' Resolve the stored ECMS cylinder geometry
#' @noRd
.ecms_resolve_geometry <- function(shape_core, elastic_comp) {
  # Resolve the effective cylinder radius from the shape metadata first ========
  radius_body <- if (!is.null(shape_core$radius)) {
    shape_core$radius
  } else {
    elastic_comp$radius
  }
  if (length(radius_body) > 1) {
    radius_body <- max(radius_body, na.rm = TRUE)
  }

  # Resolve the body length from shape metadata or the stored profile ==========
  length_body <- shape_core$length
  if (is.null(length_body)) {
    length_body <- diff(range(elastic_comp$rpos["x", ], na.rm = TRUE))
  }

  # Preserve the curvature ratio when the caller supplied a bent cylinder ======
  curvature_ratio <- if (!is.null(shape_core$radius_curvature_ratio)) {
    shape_core$radius_curvature_ratio
  } else {
    NA_real_
  }

  list(
    length = length_body,
    radius = radius_body,
    curvature_ratio = curvature_ratio
  )
}

#' Resolve the ECMS modal truncation limit
#' @noRd
.ecms_resolve_m_limit <- function(m_limit,
                                  frequency,
                                  sound_speed_sw,
                                  cL_body,
                                  cT_body,
                                  radius_body,
                                  theta_body) {
  # Preserve an explicit truncation choice when the caller supplied one ========
  if (!is.null(m_limit)) {
    return(m_limit)
  }

  # Otherwise use the standard size-parameter truncation heuristic ============
  k_sw <- wavenumber(frequency, sound_speed_sw)
  k_l <- wavenumber(frequency, cL_body)
  k_t <- wavenumber(frequency, cT_body)
  ka_max <- pmax(k_sw, k_l, k_t) * radius_body * abs(sin(theta_body))

  ceiling(ka_max) + 10
}

#' Build the stored ECMS body metadata
#' @noRd
.ecms_body_parameters <- function(geometry,
                                  elastic_comp,
                                  rho_body,
                                  cL_body,
                                  cT_body) {
  # Convert curvature ratios into physical radii for the stored metadata ======
  list(
    length = geometry$length,
    radius = geometry$radius,
    theta = elastic_comp$theta,
    density = rho_body,
    sound_speed_longitudinal = cL_body,
    sound_speed_transversal = cT_body,
    is_bent = !is.null(geometry$curvature_ratio) &&
      !is.na(geometry$curvature_ratio),
    radius_curvature = if (!is.null(geometry$curvature_ratio) &&
                           !is.na(geometry$curvature_ratio)) {
      geometry$curvature_ratio * geometry$length
    } else {
      NA_real_
    }
  )
}

#' Initialize scatterer-class object for the elastic cylinder modal series
#' model.
#' @noRd
ecms_initialize <- function(object,
                            frequency,
                            sound_speed_sw = .SEAWATER_SOUND_SPEED_DEFAULT,
                            density_sw = .SEAWATER_DENSITY_DEFAULT,
                            density_body = NULL,
                            sound_speed_longitudinal_body = NULL,
                            sound_speed_transversal_body = NULL,
                            m_limit = NULL) {
  # Extract the shape metadata and active cylinder geometry ====================
  shape <- extract(object, "shape_parameters")
  shape_core <- if (methods::is(object, "ESS")) shape$shell else shape
  .ecms_validate_scope(object, shape)
  # Resolve the elastic-cylinder material properties ===========================
  material <- .ecms_resolve_materials(
    object = object,
    sound_speed_sw = sound_speed_sw,
    density_sw = density_sw,
    density_body = density_body,
    sound_speed_longitudinal_body = sound_speed_longitudinal_body,
    sound_speed_transversal_body = sound_speed_transversal_body
  )
  elastic_comp <- material$component
  rho_body <- material$density
  cL_body <- material$sound_speed_longitudinal
  cT_body <- material$sound_speed_transversal
  .ecms_warn_incidence_regime(shape_core, elastic_comp)
  # Build the acoustic size parameters and truncation ==========================
  k_sw <- wavenumber(frequency, sound_speed_sw)
  k_l <- wavenumber(frequency, cL_body)
  k_t <- wavenumber(frequency, cT_body)
  geometry <- .ecms_resolve_geometry(shape_core, elastic_comp)
  m_limit <- .ecms_resolve_m_limit(
    m_limit = m_limit,
    frequency = frequency,
    sound_speed_sw = sound_speed_sw,
    cL_body = cL_body,
    cT_body = cT_body,
    radius_body = geometry$radius,
    theta_body = elastic_comp$theta
  )
  # Store the ECMS initialization recipe =======================================
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
    body = .ecms_body_parameters(
      geometry = geometry,
      elastic_comp = elastic_comp,
      rho_body = rho_body,
      cL_body = cL_body,
      cT_body = cT_body
    )
  )

  methods::slot(object, "model")$ECMS <- data.frame(
    frequency = frequency,
    sigma_bs = rep(NA_real_, length(frequency))
  )
  # Return the initialized scatterer ===========================================
  object
}

#' @noRd
.ecms_safe_divide <- function(numerator, denominator, eps = 1e-12) {
  # Replace near-zero denominators with signed epsilons ========================
  denominator_adj <- ifelse(
    abs(denominator) < eps,
    ifelse(Re(denominator) < 0, -eps, eps),
    denominator
  )
  # Return the stabilized quotient =============================================
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
  # Reject unsupported end-on incidence ========================================
  sin_theta <- abs(sin(theta_body))
  if (sin_theta < 1e-10) {
    stop("ECMS is not implemented for end-on incidence.")
  }
  # Define the exterior and interior size parameters ===========================
  x1 <- k_l * radius_body * sin_theta
  x2 <- k_t * radius_body * sin_theta
  x3 <- k_sw * radius_body * sin_theta

  m <- 0:m_limit
  nu <- neumann(m)
  m2 <- m^2
  # Evaluate the cylindrical Bessel functions ==================================
  J1 <- jc(m, x1)
  J2 <- jc(m, x2)
  J3 <- jc(m, x3)
  Y3 <- yc(m, x3)
  # Build the elastic-cylinder phase-shift algebra =============================
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
  # Convert the phase shifts into the modal sum ================================
  cos_eta <- 1 / sqrt(1 + tan_eta^2)
  sin_eta <- tan_eta * cos_eta
  modal_sum <- sum(
    (-1)^m * nu * sin_eta * (cos_eta - 1i * sin_eta),
    na.rm = TRUE
  )
  # Apply the finite-length coherence factor ===================================
  A <- k_sw * length_body * cos(theta_body)
  sinc_factor <- if (abs(A) < 1e-10) 1 else sin(A) / A

  -(length_body / pi) * sinc_factor * modal_sum
}

#' @noRd
.ecms_equivalent_length_fresnel <- function(k1, l, a, rho_c) {
  .trcm_equivalent_length_fresnel(k1 = k1, l = l, a = a, rho_c = rho_c)
}

#' Elastic cylinder modal series solution.
#' @noRd
ECMS <- function(object) {
  # Extract the stored ECMS inputs =============================================
  model <- extract(object, "model_parameters")$ECMS
  acoustics <- model$parameters$acoustics
  body <- model$body
  medium <- model$medium
  # Evaluate the straight-cylinder backscatter amplitudes ======================
  straight_f_bs <- mapply(
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
  # Apply the bent-cylinder coherence correction when needed ===================
  f_bs <- if (isTRUE(body$is_bent)) {
    lebc <- .ecms_equivalent_length_fresnel(
      k1 = acoustics$k_sw,
      l = body$length,
      a = body$radius,
      rho_c = body$radius_curvature
    )
    lebc * straight_f_bs / body$length
  } else {
    straight_f_bs
  }
  # Store the final ECMS results ===============================================
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
