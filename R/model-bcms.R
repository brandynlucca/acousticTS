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
#' @keywords models acoustics
NULL

#' Initialize scatterer-class object for the bent cylinder modal series model.
#' @noRd
bcms_initialize <- function(object,
                            frequency,
                            boundary = NULL,
                            sound_speed_sw = .SEAWATER_SOUND_SPEED_DEFAULT,
                            density_sw = .SEAWATER_DENSITY_DEFAULT,
                            stationary_phase = FALSE,
                            m_limit = NULL) {
  shape <- extract(object, "shape_parameters")

  if (shape$shape != "Cylinder") {
    stop(
      "The bent cylinder modal series solution requires scatterer to be ",
      "shape-type 'Cylinder'. Input scatterer is shape-type '",
      shape$shape, "'."
    )
  }

  if (is.null(boundary)) {
    boundary <- switch(
      class(object),
      FLS = "liquid_filled",
      GAS = "gas_filled",
      stop(
        "BCMS currently supports 'FLS' and 'GAS' cylinder scatterers only."
      )
    )
  }

  if (!(boundary %in% c(
    "liquid_filled", "fixed_rigid", "pressure_release", "gas_filled"
  ))) {
    stop(
      "Only 'liquid_filled', 'gas_filled', 'fixed_rigid', and ",
      "'pressure_release' are available in this implementation of BCMS."
    )
  }

  body <- extract(object, "body")
  contrasts <- .derive_contrasts(body, sound_speed_sw, density_sw)
  body$h <- contrasts$h
  body$g <- contrasts$g

  if (!is.null(shape$radius_curvature_ratio) &&
      !is.na(shape$radius_curvature_ratio) &&
      abs(body$theta - pi / 2) > pi / 18) {
    warning(
      "BCMS bent-cylinder correction is intended for broadside or near-",
      "broadside incidence."
    )
  }

  if (isTRUE(stationary_phase)) {
    warning(
      "'stationary_phase' is ignored for BCMS. BCMS uses the Fresnel-integral ",
      "bent-cylinder modal-series form only."
    )
  }

  medium_params <- data.frame(
    sound_speed = sound_speed_sw,
    density = density_sw
  )

  if (is.null(m_limit)) {
    m_limit <- ceiling(
      wavenumber(frequency, sound_speed_sw) * shape$radius
    ) + 10
  }

  methods::slot(object, "model_parameters")$BCMS <- list(
    parameters = list(
      acoustics = data.frame(
        frequency = frequency,
        lambda = sound_speed_sw / frequency,
        k_sw = wavenumber(frequency, sound_speed_sw),
        k_f = wavenumber(frequency, sound_speed_sw * body$h),
        m_limit = m_limit
      ),
      Bm_method = switch(
        boundary,
        liquid_filled = "Bm_fluid",
        gas_filled = "Bm_fluid",
        fixed_rigid = "Bm_rigid",
        pressure_release = "Bm_pressure_release"
      ),
      boundary = boundary
    ),
    medium = medium_params,
    body = list(
      length = shape$length,
      radius = shape$radius,
      theta = body$theta,
      g = body$g,
      h = body$h,
      is_bent = !is.null(shape$radius_curvature_ratio) &&
        !is.na(shape$radius_curvature_ratio),
      radius_curvature = if (!is.null(shape$radius_curvature_ratio) &&
                             !is.na(shape$radius_curvature_ratio)) {
        shape$radius_curvature_ratio * shape$length
      } else {
        NA_real_
      }
    )
  )

  methods::slot(object, "model")$BCMS <- data.frame(
    frequency = frequency,
    sigma_bs = rep(NA_real_, length(frequency))
  )

  object
}

#' @noRd
.bcms_straight_fbs <- function(acoustics, body, bm_method) {
  m_max <- max(acoustics$m_limit)
  nu <- neumann(0:m_max)

  k1L <- body$length * acoustics$k_sw
  k1a <- acoustics$k_sw * sin(body$theta) * body$radius
  k2a <- acoustics$k_sw * sin(body$theta) / body$h * body$radius
  gh <- body$g * body$h

  Bm <- switch(
    bm_method,
    Bm_fluid = .fcms_bm_fluid(k1a, k2a, gh, nu, acoustics$m_limit),
    Bm_rigid = .fcms_bm_fixed_rigid(k1a, nu, acoustics$m_limit),
    Bm_pressure_release = .fcms_bm_pressure_release(k1a, nu, acoustics$m_limit)
  )

  if (!is.matrix(Bm)) {
    Bm <- t(as.matrix(Bm))
  }

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
  model <- extract(object, "model_parameters")$BCMS
  acoustics <- model$parameters$acoustics
  body <- model$body

  straight_fbs <- .bcms_straight_fbs(
    acoustics = acoustics,
    body = body,
    bm_method = model$parameters$Bm_method
  )

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
