################################################################################
# Finite cylinder modal series solution
################################################################################
#' Finite cylinder modal series (FCMS) solution
#'
#' @description
#' Calculates the far-field scattering amplitude and related quantities for a
#' finite cylinder using the modal series solution, supporting various boundary
#' conditions (rigid, pressure-release, liquid-filled, and gas-filled).
#'
#' @section Usage:
#' This model is accessed via:
#' \preformatted{
#' target_strength(
#'   ...,
#'   model="fcms",
#'   boundary,
#'   sound_speed_sw,
#'   density_sw,
#'   m_limit
#' )
#' }
#'
#' @section Arguments:
#' \describe{
#'   \item{\code{boundary}}{Boundary condition at a cylindrical surface.
#'     One of \code{"fixed_rigid"}, \code{"pressure_release"},
#'     \code{"liquid_filled"}, or \code{"gas_filled"}. See the
#'     boundary conditions documentation for more
#'     details on these different boundary conditions.}
#'   \item{\code{sound_speed_sw}}{Seawater sound speed (\eqn{m~s^{-1}}).}
#'   \item{\code{density_sw}}{Seawater density (\eqn{kg~m^{-3}}).}
#'   \item{\code{m_limit}}{Optional model truncation limit used to cap the
#'   number of modes in the numerical calculation.}
#' }
#'
#'
#' @section Theory:
#' The modal series solution for a finite cylinder expresses the backscattering
#' amplitude as:
#'
#' \deqn{
#'   f_{bs} = -\frac{L}{\pi} \frac{\sin(k L \cos \theta)}{k L \cos \theta}
#'   \sum_{m=0}^{\infty} i^{m+1} B_m
#' }
#'
#' where \eqn{L} is the cylinder length, \eqn{k} is the wavenumber in the
#' surrounding medium, and \eqn{\theta} is the angle between the cylinder axis
#' and the incident wave direction. The coefficients \eqn{B_m} depend on the
#' boundary condition at the cylinder surface.
#'
#' **Boundary Conditions and Modal Coefficients**
#'
#' - **Rigid (fixed) cylinder:** The normal velocity at the surface is zero.
#' The modal coefficient is:
#'   \deqn{
#'     B_m = (-1)^m \epsilon_m \frac{J_m'(K a)}{H_m^{(1)'}(K a)}
#'   }
#'   where \eqn{J_m} and \eqn{H_m^{(1)}} are the cylindrical Bessel and Hankel
#'   functions of order \eqn{m}, the prime denotes differentiation with respect
#'   to the argument, \eqn{K = k \sin \theta}, and \eqn{a} is the cylinder
#'   radius. The Neumann factor is \eqn{\epsilon_0 = 1}, \eqn{\epsilon_m = 2}
#'   for \eqn{m \geq 1}.
#'
#' - **Pressure-release cylinder:** The acoustic pressure at the surface is
#' zero. The modal coefficient is:
#'   \deqn{
#'     B_m = (-1)^m \epsilon_m \frac{J_m(K a)}{H_m^{(1)}(K a)}
#'   }
#'
#' - **Fluid-filled (or gas-filled) cylinder:** Both pressure and normal
#' velocity are nonzero at the surface. The modal coefficient is:
#'   \deqn{
#'     B_m = -\epsilon_m / (1 + i C_m)
#'   }
#'   where
#'   \deqn{
#'     C_m = \frac{
#'       \left[ J_m'(K' a) Y_m(K a) \right] / \left[ J_m(K' a) J_m'(K a) \right]
#'       - g h \left[ Y_m'(K a) / J_m'(K a) \right]
#'     }{
#'       \left[ J_m'(K' a) J_m(K a) \right] / \left[ J_m(K' a) J_m'(K a) \right]
#'       - g h
#'     }
#'   }
#'   Here, \eqn{Y_m} is the cylindrical Bessel function of the second kind,
#'   \eqn{K' = K / h}, \eqn{g} is the density contrast (target to medium), and
#'   \eqn{h} is the sound speed contrast (target to medium).
#'
#' **Modal Truncation**
#'
#' The modal sum is truncated at a maximum order determined by
#' \eqn{m_{\max} = \max(\lceil k a \rceil) + 10}, which is sufficient for
#' convergence in most practical cases.
#'
#' @seealso
#' \code{\link{target_strength}}, \code{\link{FLS}}, \code{\link{GAS}},
#' \code{\link{Cylinder}}, \code{\link{cylinder}}
#'
#' @references
#'
#' Stanton, T.K. (1988). Sound scattering by cylinders of finite length. I.
#' Fuid cylinders. The Journal of the Acoustical Society of America, 83: 55-63.
#'
#' Stanton, T.K. (1989). Sound scattering by cylinders of finite length. III.
#' Deformed cylinders. The Journal of the Acoustical Society of America, 85:
#' 232-237.
#'
#' @name FCMS
#' @aliases fcms FCMS
#' @docType data
#' @keywords models acoustics internal
NULL

#' Initialize Scatterer-class object for the modal series solution for a finite
#' cylinder
#' @param object Scatterer-class object.
#' @param frequency Frequency vector (Hz).
#' @param boundary Boundary condition.
#' @param sound_speed_sw Seawater sound speed (m/s).
#' @param density_sw Seawater density (kg/m³).
#' @noRd
fcms_initialize <- function(object,
                            frequency,
                            boundary = "liquid_filled",
                            sound_speed_sw = 1500,
                            density_sw = 1026,
                            m_limit = NULL) {
  # Detect object class ========================================================
  scatterer_type <- class(object)
  # Detect object shape ========================================================
  scatterer_shape <- extract(object, "shape_parameters")
  # Validate shape +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  if (scatterer_shape$shape != "Cylinder") {
    stop(
      "The modal series solution for a finite cylinder requires scatterer to ",
      "be shape-type 'Cylinder'. Input scatterer is shape-type ",
      paste0("'", scatterer_shape, "'.")
    )
  }
  # Parse body =================================================================
  body <- .hydrate_contrasts(
    extract(object, "body"),
    sound_speed_sw, density_sw
  )
  # Define model parameters recipe =============================================
  model_params <- list(
    acoustics = .init_acoustics_df(
      frequency,
      k_sw = sound_speed_sw,
      k_f = body$h * sound_speed_sw
    )
  )
  # Determine expansion coefficient Bm method ==================================
  # Validate method ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  if (!(boundary %in% c(
    "liquid_filled", "fixed_rigid", "pressure_release", "gas_filled"
  ))) {
    stop(
      "Only the following values for 'boundary' are available in this ",
      "implementation of the finite cylinder modal series solution: ",
      "'liquid_filled' (default), 'gas_filled', 'fixed_rigid',
      'pressure_release'."
    )
  }
  # Assign method ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  model_params$Bm_method <- switch(boundary,
    liquid_filled = "Bm_fluid",
    gas_filled = "Bm_fluid",
    fixed_rigid = "Bm_rigid",
    pressure_release = "Bm_pressure_release"
  )
  # Compute body parameters ====================================================
  body_params <- list(
    # Prolate spheroidal coordinate 'xi' +++++++++++++++++++++++++++++++++++++++
    length_body = scatterer_shape$length,
    # Radius 'radius' ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    radius_body = scatterer_shape$radius,
    # Theta angle 'theta' (incident direction) +++++++++++++++++++++++++++++++++
    theta_body = body$theta,
    # Material properties
    g_body = body$g,
    h_body = body$h
  )
  # Define limits for 'm' modal series iterator ================================
  # ka + 10 seems to be an appropriate modal truncation ++++++++++++++++++++++++
  if (!is.null(m_limit)) {
    model_params$acoustics$m_limit <- m_limit
  } else {
    model_params$acoustics$m_limit <- ceiling(
      model_params$acoustics$k_sw * scatterer_shape$radius
    ) + 10
  }
  .init_model_slots(
    object = object,
    model_name = "FCMS",
    frequency = frequency,
    model_parameters = list(
      parameters = model_params,
      medium = .init_medium_params(sound_speed_sw, density_sw),
      body = body_params
    )
  )
}

#' Finite cylinder modal series (FCMS) solution
#'
#' Calculates the far-field scattering amplitude and related quantities for a
#' finite cylinder using the modal series solution, supporting various boundary
#' conditions.
#'
#' @param object Scatterer-object with a Cylinder-class shape.
#' @noRd
FCMS <- function(object) {
  # Extract model parameters/inputs ============================================
  model_params <- acousticTS::extract(object, "model_parameters")$FCMS
  parameters <- model_params$parameters
  acoustics <- parameters$acoustics
  medium <- model_params$medium
  body <- model_params$body
  # Pull out maximum 'm' specifically ==========================================
  m_max <- max(acoustics$m_limit)
  # Pre-allocate Neumann factors ===============================================
  nu <- neumann(0:m_max)
  # Compute kL and ka ==========================================================
  k1L <- body$length_body * acoustics$k_sw
  k1a <- acoustics$k_sw * sin(body$theta) * body$radius_body
  k2a <- acoustics$k_sw * sin(body$theta) / body$h * body$radius_body
  # Combine material properties ================================================
  gh <- body$g * body$h
  # Resolve modal series coefficient calculation method ========================
  Bm <- switch(parameters$Bm_method,
    Bm_fluid = .fcms_bm_fluid(k1a, k2a, gh, nu, acoustics$m_limit),
    Bm_rigid = .fcms_bm_fixed_rigid(k1a, nu, acoustics$m_limit),
    Bm_pressure_release = .fcms_bm_pressure_release(k1a, nu, acoustics$m_limit)
  )
  # Convert to matrix if needed ================================================
  if (!is.matrix(Bm)) {
    Bm <- t(as.matrix(Bm))
  }
  # Compute the linear scattering coefficient, f_bs ============================
  if (parameters$Bm_method == "Bm_fluid") {
    f_bs <- -body$length_body / pi *
      sin(k1L * cos(body$theta_body)) / (k1L * cos(body$theta_body)) *
      colSums(Bm, na.rm = TRUE)
  } else {
    f_bs <- 1i * body$length_body / pi *
      sin(k1L * cos(body$theta_body)) / (k1L * cos(body$theta_body)) *
      colSums(Bm, na.rm = TRUE)
  }
  # Calculate backscatter and return ===========================================
  # Compute sigma_bs +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  o_bs <- abs(f_bs)
  methods::slot(object, "model")$FCMS <- data.frame(
    frequency = acoustics$frequency,
    f_bs = f_bs,
    sigma_bs = o_bs,
    TS = 20 * log10(o_bs)
  )
  object
}

#' Helper function that calculates the boundary conditions for a fixed rigid
#' finite cylinder.
#' @noRd
.fcms_bm_fixed_rigid <- function(k1a, nu, m_limit) {
  # Precompute the rigid-cylinder modal weights ================================
  weights <- (-1)^(0:max(m_limit)) * nu
  # Evaluate the rigid boundary-condition coefficients =========================
  .modal_series_apply(
    m_limit = m_limit,
    FUN = function(k1a, ml) {
      jcd(0:ml, k1a) / hcd(0:ml, k1a)
    },
    k1a = k1a
  ) * weights
}


#' Helper function that calculates the boundary conditions for a
#' pressure-release finite cylinder.
#' @noRd
.fcms_bm_pressure_release <- function(k1a, nu, m_limit) {
  # Precompute the soft-cylinder modal weights =================================
  weights <- (-1)^(0:max(m_limit)) * nu
  # Evaluate the pressure-release coefficients =================================
  .modal_series_apply(
    m_limit = m_limit,
    FUN = function(k1a, ml) {
      jc(0:ml, k1a) / hc(0:ml, k1a)
    },
    k1a = k1a
  ) * weights
}

#' Helper function that calculates the boundary conditions for a fluid-filled
#' finite cylinder.
#' @noRd
.fcms_bm_fluid <- function(k1a, k2a, gh, nu, m_limit) {
  # Resolve the largest retained modal order ===================================
  m_max <- max(m_limit)
  # Build the fluid-cylinder Cm coefficients ===================================
  Cm <- .modal_series_apply(
    m_limit = m_limit,
    FUN = function(k1a, k2a, gh, ml) {
      m <- 0:ml
      cm_num <- (jcd(m, k2a) * yc(m, k1a)) / (jc(m, k2a) * jcd(m, k1a)) -
        gh * (ycd(m, k1a) / jcd(m, k1a))
      cm_denom <- (jcd(m, k2a) * jc(m, k1a)) /
        (jc(m, k2a) * jcd(m, k1a)) - gh
      cm_num / cm_denom
    },
    k1a = k1a,
    k2a = k2a,
    gh = gh
  )
  # Convert Cm into the retained Bm coefficients ===============================
  1i^((0:m_max) + 1) * (-nu * 1i^(0:m_max) / (1 + 1i * Cm))
}
