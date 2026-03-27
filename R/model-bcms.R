################################################################################
# Bent cylinder modal series solution
################################################################################
#' Bent cylinder modal series (BCMS) solution
#'
#' @description
#' Computes backscatter from straight and uniformly bent finite cylinders by
#' combining the exact finite-cylinder modal-series kernel with the equivalent
#' coherent-length correction used for uniformly bent cylinders near normal
#' incidence.
#'
#' @section Usage:
#' This model is accessed via:
#' \preformatted{
#' target_strength(
#'   ...,
#'   model = "BCMS",
#'   boundary,
#'   sound_speed_sw,
#'   density_sw,
#'   m_limit
#' )
#' }
#'
#' @section Arguments:
#' \describe{
#'   \item{\code{boundary}}{Boundary condition at the cylinder surface. One of
#'   \code{"fixed_rigid"}, \code{"pressure_release"},
#'   \code{"liquid_filled"}, or \code{"gas_filled"}.}
#'   \item{\code{sound_speed_sw}}{Seawater sound speed (\eqn{m~s^{-1}}).}
#'   \item{\code{density_sw}}{Seawater density (\eqn{kg~m^{-3}}).}
#'   \item{\code{m_limit}}{Optional model truncation limit used to cap the
#'   number of retained cylindrical modes.}
#' }
#'
#' @references
#' Stanton, T.K. (1988). Sound scattering by cylinders of finite length. I.
#' Fluid cylinders. \emph{The Journal of the Acoustical Society of America},
#' 83: 55-63.
#'
#' Stanton, T.K. (1989). Sound scattering by cylinders of finite length. III.
#' Deformed cylinders. \emph{The Journal of the Acoustical Society of America},
#' 85: 232-237.
#'
#' Stanton, T.K., Chu, D., Wiebe, P.H., and Clay, C.S. (1993). Average echoes
#' from randomly oriented random-length finite cylinders: zooplankton models.
#' \emph{The Journal of the Acoustical Society of America}, 94: 3463-3472.
#'
#' @name BCMS
#' @aliases bcms BCMS
#' @docType data
#' @keywords models acoustics internal
NULL

#' Initialize scatterer-class object for the bent cylinder modal series model.
#' @noRd
.bcms_resolve_boundary <- function(object, boundary) {
  # Derive the default boundary from the scatterer class when omitted ==========
  if (is.null(boundary)) {
    return(switch(class(object),
      FLS = "liquid_filled",
      GAS = "gas_filled",
      stop("BCMS currently supports 'FLS' and 'GAS' cylinder scatterers only.")
    ))
  }

  # Reject unsupported boundary aliases ========================================
  if (!(boundary %in% c(
    "liquid_filled", "fixed_rigid", "pressure_release", "gas_filled"
  ))) {
    stop(
      "Only 'liquid_filled', 'gas_filled', 'fixed_rigid', and ",
      "'pressure_release' are available in this implementation of BCMS."
    )
  }

  boundary
}

#' Warn when BCMS inputs leave the intended approximation regime
#' @noRd
.bcms_warn_approximation_regime <- function(shape, body, stationary_phase) {
  # Warn when the bent-cylinder correction is used away from broadside =========
  if (!is.null(shape$radius_curvature_ratio) &&
    !is.na(shape$radius_curvature_ratio) &&
    abs(body$theta - pi / 2) > pi / 18) {
    warning(
      "BCMS bent-cylinder correction is intended for broadside or near-",
      "broadside incidence."
    )
  }

  # Remind callers that the stationary-phase flag is inactive ==================
  if (isTRUE(stationary_phase)) {
    warning(
      "'stationary_phase' is ignored for BCMS. BCMS uses the Fresnel-integral ",
      "bent-cylinder modal-series form only."
    )
  }
}

#' Resolve the BCMS modal truncation limit
#' @noRd
.bcms_resolve_m_limit <- function(m_limit, frequency, sound_speed_sw, radius) {
  # Preserve an explicit truncation choice when the caller supplied one ========
  if (!is.null(m_limit)) {
    return(m_limit)
  }

  # Otherwise use the package default size-based truncation ====================
  ceiling(wavenumber(frequency, sound_speed_sw) * radius) + 10
}

#' Build the stored BCMS body metadata
#' @noRd
.bcms_body_parameters <- function(shape, body) {
  # Resolve the stored curvature metadata from the shape description ============
  curvature_ratio <- shape$radius_curvature_ratio

  list(
    length = shape$length,
    radius = shape$radius,
    theta = body$theta,
    g = body$g,
    h = body$h,
    is_bent = !is.null(curvature_ratio) && !is.na(curvature_ratio),
    radius_curvature = if (!is.null(curvature_ratio) && !is.na(curvature_ratio)) {
      curvature_ratio * shape$length
    } else {
      NA_real_
    }
  )
}

#' Initialize scatterer-class object for the bent cylinder modal series model.
#' @noRd
bcms_initialize <- function(object,
                            frequency,
                            boundary = NULL,
                            sound_speed_sw = .SEAWATER_SOUND_SPEED_DEFAULT,
                            density_sw = .SEAWATER_DENSITY_DEFAULT,
                            stationary_phase = FALSE,
                            m_limit = NULL) {
  # Extract the stored shape metadata ==========================================
  shape <- extract(object, "shape_parameters")
  # Validate the required cylinder geometry ====================================
  if (shape$shape != "Cylinder") {
    stop(
      "The bent cylinder modal series solution requires scatterer to be ",
      "shape-type 'Cylinder'. Input scatterer is shape-type '",
      shape$shape, "'."
    )
  }
  # Resolve the supported boundary alias =======================================
  boundary <- .bcms_resolve_boundary(object, boundary)
  # Hydrate the body contrasts and warn on limited regimes =====================
  body <- .hydrate_contrasts(
    extract(object, "body"),
    sound_speed_sw, density_sw
  )
  .bcms_warn_approximation_regime(shape, body, stationary_phase)
  # Resolve the modal truncation limit =========================================
  m_limit <- .bcms_resolve_m_limit(
    m_limit = m_limit,
    frequency = frequency,
    sound_speed_sw = sound_speed_sw,
    radius = shape$radius
  )
  # Store the BCMS initialization recipe =======================================
  .init_model_slots(
    object = object,
    model_name = "BCMS",
    frequency = frequency,
    model_parameters = list(
      parameters = list(
        acoustics = transform(
          .init_acoustics_df(
            frequency,
            k_sw = sound_speed_sw,
            k_f = sound_speed_sw * body$h
          ),
          lambda = sound_speed_sw / frequency,
          m_limit = m_limit
        ),
        Bm_method = switch(boundary,
          liquid_filled = "Bm_fluid",
          gas_filled = "Bm_fluid",
          fixed_rigid = "Bm_rigid",
          pressure_release = "Bm_pressure_release"
        ),
        boundary = boundary
      ),
      medium = .init_medium_params(sound_speed_sw, density_sw),
      body = .bcms_body_parameters(shape, body)
    )
  )
}

#' @noRd
.bcms_straight_fbs <- function(acoustics, body, bm_method) {
  # Precompute the modal indexing terms ========================================
  m_max <- max(acoustics$m_limit)
  nu <- neumann(0:m_max)
  # Build the angular and size parameters ======================================
  k1L <- body$length * acoustics$k_sw
  k1a <- acoustics$k_sw * sin(body$theta) * body$radius
  k2a <- acoustics$k_sw * sin(body$theta) / body$h * body$radius
  gh <- body$g * body$h
  # Evaluate the boundary-condition coefficients ===============================
  Bm <- switch(bm_method,
    Bm_fluid = .fcms_bm_fluid(k1a, k2a, gh, nu, acoustics$m_limit),
    Bm_rigid = .fcms_bm_fixed_rigid(k1a, nu, acoustics$m_limit),
    Bm_pressure_release = .fcms_bm_pressure_release(k1a, nu, acoustics$m_limit)
  )
  # Coerce the modal coefficients to a matrix ==================================
  if (!is.matrix(Bm)) {
    Bm <- t(as.matrix(Bm))
  }
  # Combine the coherent prefactor and modal sum ===============================
  prefactor <- body$length / pi *
    sin(k1L * cos(body$theta)) / (k1L * cos(body$theta))

  if (bm_method == "Bm_fluid") {
    -prefactor * colSums(Bm, na.rm = TRUE)
  } else {
    1i * prefactor * colSums(Bm, na.rm = TRUE)
  }
}

#' @noRd
.bcms_equivalent_length_fresnel <- function(k1, l, a, rho_c) {
  .trcm_equivalent_length_fresnel(k1 = k1, l = l, a = a, rho_c = rho_c)
}

#' Bent cylinder modal series solution.
#' @noRd
BCMS <- function(object) {
  # Extract the stored BCMS inputs =============================================
  model <- extract(object, "model_parameters")$BCMS
  acoustics <- model$parameters$acoustics
  body <- model$body
  # Evaluate the straight-cylinder reference amplitude =========================
  straight_fbs <- .bcms_straight_fbs(
    acoustics = acoustics,
    body = body,
    bm_method = model$parameters$Bm_method
  )
  # Apply the bent-cylinder coherence correction when needed ===================
  f_bs <- if (body$is_bent) {
    lebc <- .bcms_equivalent_length_fresnel(
      k1 = acoustics$k_sw,
      l = body$length,
      a = body$radius,
      rho_c = body$radius_curvature
    )
    lebc * straight_fbs / body$length
  } else {
    straight_fbs
  }
  # Store the final BCMS results ===============================================
  sigma_bs <- abs(f_bs)^2

  methods::slot(object, "model")$BCMS <- data.frame(
    frequency = acoustics$frequency,
    ka = acoustics$k_sw * body$radius * sin(body$theta),
    f_bs = f_bs,
    sigma_bs = sigma_bs,
    TS = 10 * log10(sigma_bs)
  )

  object
}
