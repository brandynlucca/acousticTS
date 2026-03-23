################################################################################
# Viscous-elastic spherical scattering model
################################################################################
#' Viscous-elastic spherical model (VESMS)
#'
#' @description
#' Computes backscatter from a gas-filled spherical inclusion surrounded by an
#' elastic shell and an outer viscous layer, following the wideband
#' backscattering model used by Khodabandeloo et al. (2021). The model is
#' parameterized on `ESS` objects, where the `fluid` slot represents the inner
#' gas sphere and the `shell` slot represents the elastic shell. The external
#' viscous layer is supplied as model-specific arguments.
#'
#' @section Usage:
#' This model is accessed via:
#' \preformatted{
#' target_strength(
#'   ...,
#'   model = "VESMS",
#'   sound_speed_sw,
#'   density_sw,
#'   sound_speed_viscous,
#'   density_viscous,
#'   shear_viscosity_viscous,
#'   bulk_viscosity_viscous,
#'   radius_viscous,
#'   viscous_thickness,
#'   m_limit
#' )
#' }
#'
#' @section Arguments:
#' \describe{
#'   \item{\code{sound_speed_sw}}{Seawater sound speed (\eqn{m~s^{-1}}).}
#'   \item{\code{density_sw}}{Seawater density (\eqn{kg~m^{-3}}).}
#'   \item{\code{sound_speed_viscous}}{Compressional sound speed in the outer
#'   viscous layer (\eqn{m~s^{-1}}).}
#'   \item{\code{density_viscous}}{Density of the outer viscous layer
#'   (\eqn{kg~m^{-3}}).}
#'   \item{\code{shear_viscosity_viscous}}{Shear viscosity of the outer viscous
#'   layer (\eqn{kg~m^{-1}~s^{-1}}).}
#'   \item{\code{bulk_viscosity_viscous}}{Optional bulk viscosity of the outer
#'   viscous layer. When omitted, it defaults to the supplied shear viscosity to
#'   match the reference implementation.}
#'   \item{\code{radius_viscous}}{Optional outer radius of the viscous layer
#'   (\eqn{R_2}, m).}
#'   \item{\code{viscous_thickness}}{Optional thickness of the viscous layer
#'   relative to the shell outer radius. Supply either this or
#'   \code{radius_viscous}, not both.}
#'   \item{\code{m_limit}}{Optional truncation limit for the modal summation.
#'   A value of 2 retains the \eqn{m = 0, 1, 2} terms. When omitted, the model
#'   uses \eqn{\max(2, \mathrm{round}(k_1 R_2) + 10)} at each frequency.}
#' }
#'
#' @details
#' When \code{radius_viscous} and \code{viscous_thickness} are both omitted, the
#' model estimates the viscous-layer radius from the neutral-buoyancy relation
#' used by Khodabandeloo et al. (2021, Eq. 18).
#'
#' @references
#' Khodabandeloo, B., Agersted, M.D., Klevjer, T., Macaulay, G.J., and Melle,
#' W. (2021). Estimating target strength and physical characteristics of
#' gas-bearing mesopelagic fish from wideband in situ echoes using a
#' viscous-elastic scattering model. \emph{The Journal of the Acoustical
#' Society of America}, 149, 673-691.
#'
#' Feuillade, C., and Nero, R.W. (1998). A viscous-elastic swimbladder model for
#' describing enhanced-frequency resonance scattering from fish. \emph{The
#' Journal of the Acoustical Society of America}, 103, 3245-3255.
#'
#' Anson, D.S., and Chivers, R.C. (1993). An irregular frequencies solution to
#' acoustic scattering by elastic spheres. \emph{The Journal of the Acoustical
#' Society of America}, 93, 118-123.
#'
#' @name VESMS
#' @aliases vesms VESMS
#' @docType data
#' @keywords models acoustics
NULL

#' @noRd
.vesms_is_spherical <- function(object) {
  shape <- extract(object, "shape_parameters")
  shape_tag <- shape$shape %||% ""
  any(grepl("sphere", tolower(as.character(shape_tag))))
}

#' @noRd
.vesms_validate_outer_radius <- function(radius_viscous, radius_shell) {
  if (!is.numeric(radius_viscous) ||
      length(radius_viscous) != 1L ||
      !is.finite(radius_viscous) ||
      radius_viscous <= radius_shell) {
    stop(
      "VESMS requires the viscous-layer radius to be a finite scalar larger ",
      "than the shell radius."
    )
  }

  radius_viscous
}

#' @noRd
.vesms_neutral_radius <- function(radius_gas,
                                  density_sw,
                                  density_viscous,
                                  density_gas) {
  denom <- density_viscous - density_sw
  if (!is.finite(denom) || abs(denom) <= sqrt(.Machine$double.eps)) {
    stop(
      "Neutral-buoyancy estimation of the viscous-layer radius is undefined ",
      "when 'density_viscous' equals the surrounding-medium density. Supply ",
      "'radius_viscous' or 'viscous_thickness' explicitly."
    )
  }

  radius_gas * (1 + (density_sw - density_gas) / denom)^(1 / 3)
}

#' Initialize ESS-class object for the viscous-elastic spherical model.
#' @param object ESS-class object.
#' @param frequency Frequency vector (Hz).
#' @param sound_speed_sw Seawater sound speed (m/s).
#' @param density_sw Seawater density (kg/m^3).
#' @param sound_speed_viscous Compressional sound speed in the viscous layer.
#' @param density_viscous Density of the viscous layer.
#' @param shear_viscosity_viscous Shear viscosity of the viscous layer.
#' @param bulk_viscosity_viscous Optional bulk viscosity of the viscous layer.
#' @param radius_viscous Optional outer radius of the viscous layer.
#' @param viscous_thickness Optional viscous-layer thickness.
#' @param m_limit Optional modal truncation limit.
#' @noRd
vesms_initialize <- function(object,
                             frequency,
                             sound_speed_sw = .SEAWATER_SOUND_SPEED_DEFAULT,
                             density_sw = .SEAWATER_DENSITY_DEFAULT,
                             sound_speed_viscous,
                             density_viscous,
                             shear_viscosity_viscous,
                             bulk_viscosity_viscous = NULL,
                             radius_viscous = NULL,
                             viscous_thickness = NULL,
                             m_limit = NULL) {
  if (!methods::is(object, "ESS")) {
    stop(
      "VESMS currently requires an elastic-shelled scatterer ('ESS'). Input ",
      "scatterer is type '", class(object), "'."
    )
  }

  if (!.vesms_is_spherical(object)) {
    stop("VESMS currently requires a spherical ESS geometry.")
  }

  if (!is.null(radius_viscous) && !is.null(viscous_thickness)) {
    stop(
      "Specify at most one of 'radius_viscous' or 'viscous_thickness' for ",
      "VESMS."
    )
  }

  shell <- .extract_material_props(
    extract(object, "shell"),
    sound_speed_sw,
    density_sw
  )
  fluid <- .extract_material_props(
    extract(object, "fluid"),
    sound_speed_sw,
    density_sw
  )

  shell$radius <- extract(object, "shell")$radius
  shell$shell_thickness <- extract(object, "shell")$shell_thickness
  fluid$radius <- extract(object, "fluid")$radius

  if (is.null(shell$radius) || is.null(fluid$radius) ||
      !is.finite(shell$radius) || !is.finite(fluid$radius) ||
      fluid$radius >= shell$radius) {
    stop(
      "VESMS requires a valid shell radius larger than the inner gas radius."
    )
  }

  if (is.null(shell$density) || !is.finite(shell$density) ||
      is.null(shell$G) || !is.finite(shell$G)) {
    stop(
      "VESMS requires the ESS shell to have finite density and shear modulus."
    )
  }

  if (is.null(shell$lambda) || !is.finite(shell$lambda)) {
    shell$lambda <- tryCatch(
      lame(K = shell$K, E = shell$E, G = shell$G, nu = shell$nu),
      error = function(e) NA_real_
    )
  }

  if (!is.finite(shell$lambda)) {
    stop(
      "VESMS requires Lam", "\u00E9", "'s first parameter for the shell. ",
      "Provide sufficient shell elastic properties when building the ESS ",
      "object."
    )
  }

  if (is.null(fluid$density) || !is.finite(fluid$density) ||
      is.null(fluid$sound_speed) || !is.finite(fluid$sound_speed)) {
    stop(
      "VESMS requires the ESS inner fluid slot to represent a gas core with ",
      "finite density and sound speed."
    )
  }

  if (missing(sound_speed_viscous) ||
      missing(density_viscous) ||
      missing(shear_viscosity_viscous)) {
    stop(
      "VESMS requires 'sound_speed_viscous', 'density_viscous', and ",
      "'shear_viscosity_viscous'."
    )
  }

  if (!is.finite(sound_speed_viscous) || sound_speed_viscous <= 0 ||
      !is.finite(density_viscous) || density_viscous <= 0 ||
      !is.finite(shear_viscosity_viscous) || shear_viscosity_viscous <= 0) {
    stop(
      "VESMS requires positive finite viscous-layer sound speed, density, and ",
      "shear viscosity."
    )
  }

  bulk_viscosity_viscous <- if (is.null(bulk_viscosity_viscous)) {
    shear_viscosity_viscous
  } else {
    bulk_viscosity_viscous
  }

  if (!is.finite(bulk_viscosity_viscous) || bulk_viscosity_viscous < 0) {
    stop("VESMS requires a finite non-negative bulk viscosity.")
  }

  radius_viscous <- if (!is.null(radius_viscous)) {
    .vesms_validate_outer_radius(radius_viscous, shell$radius)
  } else if (!is.null(viscous_thickness)) {
    .vesms_validate_outer_radius(shell$radius + viscous_thickness, shell$radius)
  } else {
    .vesms_validate_outer_radius(
      .vesms_neutral_radius(
        radius_gas = fluid$radius,
        density_sw = density_sw,
        density_viscous = density_viscous,
        density_gas = fluid$density
      ),
      shell$radius
    )
  }

  m_limit <- if (is.null(m_limit)) {
    pmax(2L, round(wavenumber(frequency, sound_speed_sw) * radius_viscous) + 10L)
  } else {
    as.integer(m_limit)
  }

  if (length(m_limit) == 1L) {
    m_limit <- rep(m_limit, length(frequency))
  }

  if (length(m_limit) != length(frequency) ||
      any(!is.finite(m_limit)) ||
      any(m_limit < 0)) {
    stop(
      "VESMS requires 'm_limit' to be a non-negative scalar or a vector with ",
      "one value per frequency."
    )
  }

  acoustics <- data.frame(
    frequency = frequency,
    k_sw = wavenumber(frequency, sound_speed_sw),
    m_limit = m_limit
  )

  methods::slot(object, "model_parameters")$VESMS <- list(
    parameters = list(acoustics = acoustics),
    medium = data.frame(
      sound_speed = sound_speed_sw,
      density = density_sw
    ),
    viscous = list(
      sound_speed = sound_speed_viscous,
      density = density_viscous,
      shear_viscosity = shear_viscosity_viscous,
      bulk_viscosity = bulk_viscosity_viscous,
      radius = radius_viscous
    ),
    shell = shell,
    fluid = fluid
  )

  methods::slot(object, "model")$VESMS <- data.frame(
    frequency = frequency,
    sigma_bs = rep(NA_real_, length(frequency))
  )

  object
}

#' @noRd
VESMS <- function(object) {
  model <- extract(object, "model_parameters")$VESMS
  acoustics <- model$parameters$acoustics
  medium <- model$medium
  viscous <- model$viscous
  shell <- model$shell
  fluid <- model$fluid
  f_bs <- vesms_backscatter_cpp(
    frequency = acoustics$frequency,
    m_limit = as.integer(acoustics$m_limit),
    density_sw = medium$density,
    density_viscous = viscous$density,
    density_shell = shell$density,
    density_gas = fluid$density,
    sound_speed_sw = medium$sound_speed,
    sound_speed_viscous = viscous$sound_speed,
    sound_speed_gas = fluid$sound_speed,
    radius_viscous = viscous$radius,
    radius_shell = shell$radius,
    radius_gas = fluid$radius,
    shear_viscosity_viscous = viscous$shear_viscosity,
    bulk_viscosity_viscous = viscous$bulk_viscosity,
    shear_modulus_shell = shell$G,
    lambda_shell = shell$lambda
  )

  sigma_bs <- .sigma_bs(f_bs)

  methods::slot(object, "model")$VESMS <- data.frame(
    frequency = acoustics$frequency,
    ka_viscous = acoustics$k_sw * viscous$radius,
    ka_shell = acoustics$k_sw * shell$radius,
    ka_gas = acoustics$k_sw * fluid$radius,
    f_bs = f_bs,
    sigma_bs = sigma_bs,
    TS = db(sigma_bs)
  )

  object
}
