################################################################################
# Prolate spheroid modal series solution =======================================
################################################################################
#' Prolate spheroidal modal series (PSMS) solution
#'
#' @details
#' Calculates the far-field scattering amplitude and related quantities for a
#' prolate spheroid using the modal series solution, supporting various boundary
#' conditions (rigid, pressure-release, liquid-filled, and gas-filled).
#'
#' @section Usage:
#' This model is accessed via:
#' \preformatted{
#' target_strength(
#'   ...,
#'   model="psms",
#'   phi_body,
#'   boundary,
#'   adaptive,
#'   simplify_Amn,
#'   precision,
#'   n_integration,
#'   sound_speed_sw,
#'   density_sw
#' )
#' }
#'
#' @section Arguments:
#' \describe{
#'   \item{\code{phi_body}}{Incident roll angle (radians).}
#'   \item{\code{boundary}}{Boundary condition at the spheroid surface.
#'     One of \code{"fixed_rigid"}, \code{"pressure_release"},
#'     \code{"liquid_filled"}, or \code{"gas_filled"}. See the
#'     boundary conditions documentation for more
#'     details on these different boundary conditions.}
#'   \item{\code{adaptive}}{Apply adaptive modal support and quadrature
#'   selection. Defaults to \code{FALSE}.}
#'   \item{\code{simplify_Amn}}{When \code{FALSE} (default), solve the coupled
#'   pressure and normal-velocity boundary equations for liquid- or gas-filled
#'   spheroids. When \code{TRUE}, retain the diagonal approximation described in
#'   Theory. This approximation can differ near gas resonances. Ignored for
#'   other boundary conditions.}
#'   \item{\code{precision}}{A literal argument that allows for levels of
#'   precision when calculating the expansion matrix, \eqn{A_{mn}}. There are
#'   two choices: \code{"double"} for double-precision (64-bit) and
#'   \code{"quad"} for quadruple-precision (128-bit). Use \code{"quad"} when
#'   strong gas contrast or high acoustic size makes double precision unstable.}
#'   \item{\code{n_integration}}{Number of Gauss-Legendre integration points
#'   used to calculate angular overlap integrals for full liquid-filled
#'   solutions in either precision. With \code{adaptive = FALSE}, an explicit
#'   value is used as supplied and \code{NULL} selects automatic integration. With
#'   \code{adaptive = TRUE}, quadrature order is selected internally and an
#'   explicit value is ignored with a warning. Full gas solutions calculate
#'   overlaps from normalized expansion coefficients and do not use this
#'   quadrature count. Fixed quadrature must resolve the angular functions.
#'   Automatic quadrature integrates the retained angular expansion products
#'   using the Gauss-Legendre polynomial exactness condition.
#'   Under this condition, orthogonality allows the integral to be evaluated
#'   algebraically from the complete expansion coefficients. Fixed integration
#'   counts continue to use their requested quadrature rule.
#'   See \link{gauss_legendre}.}
#'   \item{\code{sound_speed_sw}}{Seawater sound speed (\eqn{m~s^{-1}}).}
#'   \item{\code{density_sw}}{Seawater density (\eqn{kg~m^{-3}}).}
#' }
#'
#' @section Theory:
#' some relevant parameters are used to describe the geometry. A spheroid
#' surface is given by \eqn{\xi = \xi_0 = \mathrm{constant}}. Denoting the
#' major radius \eqn{a}, the minor radius \eqn{b}, and the semi-focal-length
#' \eqn{q}, then:
#'
#' \deqn{
#'   a = \xi_0 q \\
#'   \xi_0 = \left[1 - \left(\frac{b}{a}\right)^2\right]^{-1/2}
#' }
#'
#' Densities are expressed by \eqn{\rho}, sound speeds by \eqn{c}, and
#' wavenumbers by \eqn{k}, each followed by a subscript 0 for the surrounding
#' medium or 1 for the spheroidal body. For the spheroid, an alternative of the
#' reduced frequency \eqn{ka} is \eqn{h = kq}.
#'
#' The far-field scattering amplitude is given by:
#'
#' \deqn{
#'   f_\infty(\theta, \phi | \theta', \phi') =
#'   \frac{2}{j k_0} \sum\limits_{m=0}^\infty \sum\limits_{n=m}^\infty
#'   \frac{\epsilon_m}{N_{mn}(h_0)} S_{mn}(h_0, \cos\theta') \,
#'   A_{mn} S_{mn}(h_0, \cos\theta) \cos m(\phi - \phi')
#' }
#'
#' where \eqn{\epsilon_m = (-1)^{m/2}} for even m, and \eqn{N_{mn}(h_0)} is the
#' norm. \eqn{S_{mn}} is the angle spheroidal wave function of the first kind
#' of order \eqn{m} and degree \eqn{n}, and \eqn{A_{mn}} is the expansion
#' coefficient for the scattered wave, determined by the boundary conditions.
#' The parameters \eqn{\theta} and \eqn{\phi} denote the spherical angle
#' coordinates of an observed point along the body. The \eqn{\theta'} and
#' \eqn{\phi'} denote similar spherical angle coordinates of the incident
#' direction. Effectively, \eqn{\phi} and \eqn{\theta} correspond to the roll
#' and tilt angles (relative the incident sound wave), respectively. It is
#' assumed that \eqn{\theta'} and \eqn{\phi'} are perpendicular to the
#' incident wave such that:
#'
#' \deqn{
#'   \theta' = \pi - \theta \\
#'   \phi' = \pi + \phi
#' }
#'
#' For pressure release (or soft) and rigid spheroids:
#'
#' \deqn{
#'   A_{mn} = \frac{-\Delta R_{mn}^{(1)}(h_0, \xi_0)}{
#'   \Delta R_{mn}^{(3)}(h_0, \xi_0)
#'   }
#' }
#'
#' where \eqn{\Delta = 1} for soft and \eqn{\Delta = \partial / \partial \xi} f
#' or rigid spheroid, and \eqn{R_{mn}^{(i)}} is the radial spheroidal wave
#' function of the i-th kind.
#'
#' For the case of a fluid-filled spheroid, the following simultaneous equation
#' must be solved:
#'
#' \deqn{
#'   \sum\limits_{n'=m}^\infty K_{nn'}^{(3)} A_{mn'} + \sum\limits_{n'=m}^\infty
#'   K_{nn'}^{(1)} \alpha_{mn'}^{m} E_{n'}^{m(1)} = 0
#' }
#'
#' where
#'
#' \deqn{
#'   K_{nl}^{(i)} = \frac{1}{N_{mn}(h_0)}
#'   \int\limits_{-1}^{1} S_{mn}(h_0, \cos\theta)
#'   S_{ml}(h_1, \eta) d\eta
#' }
#'
#' and
#'
#' \deqn{
#'   \alpha_{mn}^{m} = \frac{1}{N_{mn}(h_1)}
#'   \int\limits_{-1}^{1} S_{mn}(h_0, \eta)
#'   S_{mn}(h_1, \eta) d\eta
#' }
#'
#' In practice, these kernel matrices are solved numerically to determine
#' \eqn{A_{mn}}. The implementation uses compiled dense linear algebra and
#' keeps the fluid-filled solve in the requested arithmetic so that the
#' linear-system stage does not become detached from the rest of the chosen
#' precision pathway.
#'
#' For full gas solutions, angular overlaps are dot products
#' of expansion coefficients in an orthonormal associated Legendre basis.
#' This follows from the spheroidal expansion and Legendre orthogonality
#' identities in \href{https://dlmf.nist.gov/30.8.E1}{DLMF 30.8.1} and
#' \href{https://dlmf.nist.gov/14.17.E6}{DLMF 14.17.6}. The coefficients and
#' overlap sums retain the requested arithmetic precision.
#'
#' In the case that \eqn{h_0 \approx h_1}, \eqn{\alpha_{mn}^{m} \approx 0} for
#' \eqn{n \neq l}, and a much simplified expression is derived:
#'
#' \deqn{
#'   A_{mn} = -\frac{E_{mn}^{m(1)}}{E_{mn}^{m(3)}}
#' }
#'
#' The maximum values of \eqn{m} and \eqn{n} can be estimated by:
#'
#' \deqn{
#'   m_{\max} = [2 k_0 b] \\
#'   n_{\max} = m_{\max} + [h_0 / 2]
#' }
#'
#' Full liquid and gas solutions also use the interior-wave and
#' incidence-dependent estimates accompanying Khodabandeloo et al. (2025).
#' These estimates define the retained modal ceilings. Adaptive calculations
#' can truncate negligible modal tails before reaching those ceilings.
#'
#' Full fluid solutions preserve pressure and normal-velocity continuity.
#' Radial mantissas and exponents remain separate during assembly to avoid
#' intermediate overflow.
#'
#' @section Implementation:
#'
#' **\code{C++}**
#'
#' This model is primarily implemented in \code{C++} since it leverages the
#' \code{Fortran} algorithm developed by Arnie Lee van Buren and Jeffery
#' Boisvert. Another reason this is primarily written in \code{C++} is due to
#' computational and performance reasons. As \eqn{k_0}, \eqn{b}, and \eqn{h}
#' increase, so do \eqn{m_{\max}} and \eqn{n_{\max}}. While this is not as much
#' of a concern for the pressure-release and fixed rigid cases, this results in
#' increasingly large and unwieldy kernel matrices required for solving the
#' fluid-filled expansion matrices (\eqn{A_{mn}}). Consequently, this can
#' result in the model taking an impractical amount of time to compute
#' \eqn{f_\infty(\theta, \phi | \theta', \phi')}. \code{C++} helps reduce this
#' burden compared to a pure \code{R} implementation. The current implementation
#' uses compiled matrix algebra throughout, applies backscatter parity to avoid
#' recomputing one of the two angular \eqn{S_{mn}} matrices, and solves the
#' fluid-filled kernel system natively in the requested arithmetic.
#'
#' **Precision**
#'
#' Another consideration is the floating point precision inherent to \code{R}.
#' \code{R} uses double-precision, which stores numeric values in a 64-bit
#' format. At low \eqn{m_{\max}} and \eqn{n_{\max}}, there is lower numerical
#' instability in the prolate spheroidal wave functions used in the model.
#' However, this is not the case at greater \eqn{m_{\max}} and \eqn{n_{\max}}
#' where the difference between double- (64-bit) and quadruple-precision
#' (128-bit) can contribute to differences in target strength of more than 1
#' dB. While \code{R}-packages like \code{Rmpfr} provide access to the
#' (\code{GNU}) \code{MPFR C} library, the (\code{GCC}) \code{libquadmath C++}
#' library. Quad precision nevertheless remains much more computationally
#' intensive than double precision because the spheroidal function evaluations,
#' overlap integrals, and kernel systems all grow rapidly with
#' \eqn{m_{\max}} and \eqn{n_{\max}}.
#'
#' @seealso
#' \code{\link{target_strength}}, \code{\link{FLS}}, \code{\link{GAS}},
#' \code{\link{ESS}}, \code{\link{ProlateSpheroid}},
#' \code{\link{prolate_spheroid}}, \code{\link{Smn}}, \code{\link{Rmn}}
#'
#' @references
#'
#' Furusawa, M. (1988). Prolate spheroidal models for predicting general trends
#' of fish target strength. Journal of the Acoustical Society of Japan, 9:
#' 13-24.
#'
#' Khodabandeloo, B., Heggelund, Y., Ystad, B., Marx, S. A. B., and Pedersen, G.
#' (2025). High-precision model and open-source software for acoustic
#' backscattering by liquid- and gas-filled prolate spheroids across a wide
#' frequency range and incident angles: Implications for fisheries acoustics.
#' Journal of Sound and Vibration, 616: 119227.
#' \doi{10.1016/j.jsv.2025.119227}.
#'
#' Sanderson, C., and Curtin, R. (2019). Practical sparse matrices in C++ with
#' hybrid storage and template-based expression optimisation. Mathematical and
#' Computational Applications, 24: 70.
#'
#' Spencer, R.D., and Granger, S. (1951). The scattering of sound from a
#' prolate spheroid. The Journal of the Acoustical Society of America, 23:
#' 701-706.
#'
#' Van Buren, A. L. and Boisvert, J. E. "Prolate Spheroidal Wave Functions."
#' GitHub repository:
#' \url{https://github.com/MathieuandSpheroidalWaveFunctions/Prolate_swf}
#'
#' @return No value; this help topic documents the PSMS model.
#' @examples
#' subset(available_models(), model == "psms")
#'
#' @name PSMS
#' @aliases psms PSMS
#' @keywords models acoustics internal
NULL

#' Validate that the current PSMS implementation is only used for prolates.
#' @noRd
.psms_validate_shape <- function(scatterer_shape) {
  # Restrict the modal-series initializer to prolate spheroids =================
  if (scatterer_shape$shape != "ProlateSpheroid") {
    stop(
      "The modal series solution for a prolate spheroid requires scatterer to ",
      "be shape-type 'ProlateSpheroid'. Input scatterer is shape-type ",
      paste0("'", scatterer_shape$shape, "'.")
    )
  }

  invisible(TRUE)
}

#' Validate the PSMS boundary label used to choose the Amn formulation.
#' @noRd
.psms_validate_boundary <- function(boundary) {
  # Confirm that the requested PSMS boundary is implemented ====================
  if (!(boundary %in% c(
    "liquid_filled", "gas_filled", "fixed_rigid", "pressure_release"
  ))) {
    stop(
      "Only the following values for 'method' are available in this ",
      "implementation of the prolate spheroid modal series solution: ",
      "'liquid_filled' (default), 'gas_filled', 'fixed_rigid', ",
      "'pressure_release'."
    )
  }

  boundary
}

#' Report whether this installation has the complete native quad backend.
#' @noRd
.quad_precision_available <- function() {
  isTRUE(quad_precision_available_cpp())
}

#' Reject requests for quad precision when configure could not build the
#' matching C++ and Fortran binary128 implementations.
#' @noRd
.validate_quad_precision_available <- function(precision) {
  if (identical(precision, "quad") && !.quad_precision_available()) {
    stop(
      "Native quad precision is unavailable with the C++ and Fortran ",
      "compilers used to build this installation; use precision = \"double\".",
      call. = FALSE
    )
  }

  precision
}

#' Validate the PSMS precision label.
#' @noRd
.psms_validate_precision <- function(precision) {
  # Restrict the precision choice to the supported backends ====================
  if (!precision %in% c("double", "quad")) {
    stop("'precision' must be either 'double' or 'quad'.")
  }

  .validate_quad_precision_available(precision)
}

#' Validate the PSMS adaptive-mode flag.
#' @noRd
.psms_validate_adaptive <- function(adaptive) {
  # Require a scalar logical adaptive flag =====================================
  if (!is.logical(adaptive) || length(adaptive) != 1 || is.na(adaptive)) {
    stop("'adaptive' must be either TRUE or FALSE.")
  }

  adaptive
}

#' Resolve the Amn formulation associated with one PSMS boundary condition.
#' @noRd
.psms_Amn_method <- function(boundary, simplify_Amn) {
  # Map the public PSMS boundary labels onto the internal kernel names =========
  switch(boundary,
    liquid_filled = ifelse(simplify_Amn, "Amn_fluid_simplify", "Amn_fluid"),
    gas_filled = ifelse(simplify_Amn, "Amn_fluid_simplify", "Amn_fluid_gas"),
    fixed_rigid = "Amn_fixed_rigid",
    pressure_release = "Amn_pressure_release"
  )
}

#' Identify full penetrable PSMS boundary solves.
#' @noRd
.psms_is_full_fluid_method <- function(Amn_method) {
  isTRUE(Amn_method %in% c("Amn_fluid", "Amn_fluid_gas", "Amn_fluid_coefficients"))
}

#' Resolve the hydrated PSMS body properties against the surrounding medium.
#' @noRd
.psms_body_state <- function(object, sound_speed_sw, density_sw) {
  # Hydrate contrasts into absolute body properties ============================
  body <- .hydrate_contrasts(
    extract(object, "body"),
    sound_speed_sw,
    density_sw
  )

  list(
    body = body,
    medium = .init_medium_params(sound_speed_sw, density_sw)
  )
}

#' Build the PSMS model-parameter list prior to geometry-specific bookkeeping.
#' @noRd
.psms_model_parameters <- function(frequency,
                                   sound_speed_sw,
                                   body_h,
                                   precision,
                                   n_integration,
                                   adaptive) {
  # Initialize the common acoustics table and scalar controls ==================
  list(
    acoustics = .init_acoustics_df(
      frequency,
      k_sw = sound_speed_sw,
      k_f = body_h * sound_speed_sw
    ),
    precision = precision,
    n_integration = n_integration,
    adaptive = adaptive
  )
}

#' Build the stored body metadata for one initialized PSMS object.
#' @noRd
.psms_body_parameters <- function(scatterer_shape,
                                  body,
                                  phi_body,
                                  sound_speed_sw,
                                  density_sw) {
  # Convert the stored prolate geometry into spheroidal coordinates ============
  geometry <- .prolate_spheroidal_geometry(scatterer_shape)
  body_params <- list(
    xi = geometry$xi,
    phi_body = phi_body,
    phi_scatter = phi_body + pi,
    theta_body = body$theta,
    theta_scatter = pi - body$theta,
    density = body$g * density_sw,
    sound_speed = body$h * sound_speed_sw
  )
  body_params$q <- geometry$q

  body_params
}

#' Resolve the quadrature order used by the PSMS kernels.
#' @noRd
.psms_n_integration_label <- function(Amn_method) {
  switch(Amn_method,
    Amn_fluid = "full liquid-filled PSMS solve",
    Amn_fluid_gas = "full gas-filled PSMS solve",
    "full penetrable PSMS solve"
  )
}

#' Resolve the quadrature order used by the PSMS kernels.
#' @noRd
.psms_n_integration <- function(n_integration, adaptive, Amn_method) {
  # Let the full liquid solve choose adaptive quadrature =======================
  use_adaptive_quadrature <- isTRUE(adaptive) && identical(Amn_method, "Amn_fluid")
  if (use_adaptive_quadrature) {
    if (!is.null(n_integration)) {
      warning(
        "'n_integration' is ignored when 'adaptive = TRUE' for the ",
        .psms_n_integration_label(Amn_method),
        "."
      )
    }
    return(NA_integer_)
  }

  # Integrate complete liquid expansions when no fixed count is supplied =======
  if (is.null(n_integration) && identical(Amn_method, "Amn_fluid")) {
    return(NA_integer_)
  }

  # Otherwise validate or fill in the explicit quadrature order ================
  if (is.null(n_integration)) {
    return(96L)
  }
  if (!is.numeric(n_integration) ||
    length(n_integration) != 1 ||
    !is.finite(n_integration) ||
    n_integration < 1 ||
    n_integration > .Machine$integer.max ||
    n_integration %% 1 != 0) {
    stop("'n_integration' must be a single positive integer within R's integer range.")
  }

  as.integer(n_integration)
}

#' Attach the reduced frequencies and modal truncation limits to the PSMS
#' acoustics table.
#' @noRd
.psms_complete_acoustics <- function(
  model_params, body_params,
  scatterer_shape
) {
  # Add reduced frequencies for the surrounding medium and interior ============
  model_params$acoustics$chi_sw <- model_params$acoustics$k_sw * body_params$q
  model_params$acoustics$chi_body <- model_params$acoustics$k_f * body_params$q
  # Initialize modal ceilings from the published exterior estimates ============
  geometry <- .prolate_spheroidal_geometry(scatterer_shape)
  model_params$acoustics$m_max <- ceiling(
    2 * model_params$acoustics$k_sw * geometry$semiminor_length
  )
  model_params$acoustics$n_max <- model_params$acoustics$m_max +
    ceiling(0.5 * model_params$acoustics$chi_sw)

  model_params
}

#' Promote the retained modal ceilings for full liquid and gas PSMS solves.
#' @noRd
.psms_promote_full_fluid_modal_limits <- function(
  model_params, body_params,
  boundary
) {
  full_gas <- identical(boundary, "gas_filled") &&
    identical(model_params$Amn_method, "Amn_fluid_gas")
  full_liquid <- identical(boundary, "liquid_filled") &&
    model_params$Amn_method %in% c("Amn_fluid", "Amn_fluid_coefficients")
  if (!full_gas && !full_liquid) {
    return(model_params)
  }

  broadside_weight <- abs(sin(body_params$theta_body))
  chi_body_abs <- abs(model_params$acoustics$chi_body)
  sound_speed_ratio <- ifelse(
    chi_body_abs > .Machine$double.eps,
    abs(model_params$acoustics$chi_sw / model_params$acoustics$chi_body),
    1
  )
  m_seed <- if (isTRUE(body_params$sound_speed < 500)) {
    0.5 * sound_speed_ratio * broadside_weight^3 * chi_body_abs + 8
  } else {
    0.3 * sound_speed_ratio * broadside_weight^3 * chi_body_abs + 6
  }
  # Retain the Khodabandeloo et al. (2025) software's high-accuracy estimate ===
  n_offset <- 12
  angle_discount <- 7
  n_seed <- 0.75 * chi_body_abs + n_offset - angle_discount * broadside_weight
  adaptive_m_max <- ceiling(m_seed)
  adaptive_n_max <- ceiling(n_seed)

  model_params$acoustics$m_max <- pmax(
    model_params$acoustics$m_max,
    adaptive_m_max
  )
  model_params$acoustics$n_max <- pmax(
    model_params$acoustics$n_max,
    adaptive_n_max,
    model_params$acoustics$m_max
  )

  model_params
}

#' Initialize object for the modal series solution for a prolate spheroid
#' @param object Scatterer-class object.
#' @param frequency Frequency vector (Hz).
#' @param phi_body Body rotation angle (radians).
#' @param boundary Boundary condition.
#' @param adaptive Apply adaptive modal support and quadrature
#' selection.
#' @param precision Numerical precision.
#' @param n_integration Number of Gauss-Legendre integration points for full
#' liquid-filled solutions. `NULL` selects automatic integration.
#' With `adaptive = TRUE`, the quadrature order is selected internally.
#' @param simplify_Amn Calculate the boundary expansion coefficient using a
#' simplified formulation.
#' @param sound_speed_sw Seawater sound speed (m/s).
#' @param density_sw Seawater density (kg/m³).
#' @noRd
psms_initialize <- function(object,
                            frequency,
                            phi_body = pi,
                            boundary = "liquid_filled",
                            adaptive = FALSE,
                            precision = "double",
                            n_integration = NULL,
                            simplify_Amn = FALSE,
                            sound_speed_sw = 1500,
                            density_sw = 1026) {
  # Detect object class ========================================================
  scatterer_shape <- acousticTS::extract(object, "shape_parameters")
  .psms_validate_shape(scatterer_shape)
  body_state <- .psms_body_state(object, sound_speed_sw, density_sw)
  body <- body_state$body
  medium_params <- body_state$medium
  # Define model parameters recipe =============================================
  boundary <- .psms_validate_boundary(boundary)
  precision <- .psms_validate_precision(precision)
  adaptive <- .psms_validate_adaptive(adaptive)
  model_params <- .psms_model_parameters(
    frequency = frequency,
    sound_speed_sw = sound_speed_sw,
    body_h = body$h,
    precision = precision,
    n_integration = n_integration,
    adaptive = adaptive
  )
  # Determine expansion coefficient Amn method =================================
  model_params$Amn_method <- .psms_Amn_method(boundary, simplify_Amn)
  # Compute body parameters ====================================================
  body_params <- .psms_body_parameters(
    scatterer_shape = scatterer_shape,
    body = body,
    phi_body = phi_body,
    sound_speed_sw = sound_speed_sw,
    density_sw = density_sw
  )
  model_params <- .psms_complete_acoustics(
    model_params = model_params,
    body_params = body_params,
    scatterer_shape = scatterer_shape
  )
  model_params <- .psms_promote_full_fluid_modal_limits(
    model_params = model_params,
    body_params = body_params,
    boundary = boundary
  )
  model_params$n_integration <- .psms_n_integration(
    n_integration = n_integration,
    adaptive = adaptive,
    Amn_method = model_params$Amn_method
  )
  .init_model_slots(
    object = object,
    model_name = "PSMS",
    frequency = frequency,
    model_parameters = list(
      parameters = model_params,
      medium = medium_params,
      body = body_params
    )
  )
}

#' Traditional kernel solver for the PSMS.
#' @keywords internal
#' @noRd
.prolate_spheroidal_kernels_fixed <- function(
  acoustics,
  body,
  medium,
  boundary_method,
  n_integration = 96,
  precision = "double",
  adaptive = FALSE
) {
  # Generate nodes and weights for quadrature ==================================
  quad_pts <- if (boundary_method %in% c("Amn_fluid_gas",
    "Amn_fluid_coefficients") ||
    (isTRUE(adaptive) && identical(boundary_method, "Amn_fluid"))) {
    list(nodes = numeric(), weights = numeric())
  } else {
    gauss_legendre(n = n_integration, a = -1, b = 1)
  }
  # Retain the requested count while refining Gauss points in native precision =
  attr(quad_pts, "acousticTS_gauss_legendre") <- TRUE
  # Calculate the linear scattering coefficient, fbs ===========================
  prolate_spheroid_fbs(
    acoustics, body, medium, quad_pts, precision, boundary_method, adaptive,
    FALSE
  )
}

#' Solve Expansion Coefficients for a Liquid-Filled Spheroidal Scatterer
#'
#' @description
#' Computes the modal expansion coefficients \eqn{A_{mn}} for acoustic
#' scattering from a liquid-filled prolate spheroidal body by matching pressure
#' and normal velocity in the requested arithmetic precision.
#'
#' @details
#' This function solves the linear system arising from matching boundary
#' conditions at the surface of a liquid-filled spheroidal scatterer. The
#' expansion coefficients  \eqn{A_{mn}} relate the scattered field to the
#' incident field through the spheroidal wave function expansion.
#' @keywords internal
#' @noRd
prolate_spheroidal_kernels <- function(
  acoustics,
  body,
  medium,
  boundary_method,
  n_integration = 96,
  precision = "double",
  adaptive = FALSE
) {
  # Resolve automatic liquid integrals through Ferrers orthogonality ===========
  if (!isTRUE(adaptive) && identical(boundary_method, "Amn_fluid") &&
    (is.null(n_integration) || (length(n_integration) == 1 &&
      is.na(n_integration)))) {
    boundary_method <- "Amn_fluid_coefficients"
  }
  # Normalize missing quadrature input =========================================
  if (is.null(n_integration) || (length(n_integration) == 1 &&
    is.na(n_integration))) {
    n_integration <- 96L
  }
  # Fall back to the fixed-order kernel solver =================================
  .prolate_spheroidal_kernels_fixed(
    acoustics = acoustics,
    body = body,
    medium = medium,
    boundary_method = boundary_method,
    n_integration = n_integration,
    precision = precision,
    adaptive = adaptive
  )
}

#' Calculates the theoretical target strength of a prolate spheroid using a
#' modal series solution
#' @noRd
PSMS <- function(object) {
  # Extract model parameters/inputs ============================================
  model_params <- extract(object, "model_parameters")$PSMS
  parameters <- model_params$parameters
  acoustics <- parameters$acoustics
  medium <- model_params$medium
  body <- model_params$body
  # Compute the linear scattering coefficient, f_bs ============================
  f_bs <- prolate_spheroidal_kernels(
    acoustics, body, medium, parameters$Amn_method, parameters$n_integration,
    parameters$precision, parameters$adaptive
  )
  # Calculate backscatter and return ===========================================
  # Compute sigma_bs +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  o_bs <- abs(-2i / acoustics$k_sw * f_bs)
  methods::slot(object, "model")$PSMS <- data.frame(
    frequency = acoustics$frequency,
    f_bs = f_bs,
    sigma_bs = o_bs,
    TS = 20 * log10(o_bs)
  )
  object
}
