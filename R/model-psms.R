################################################################################
# Prolate spheroid modal series solution
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
#'   \item{\code{boundary}}{Boundary condition at the cylinder surface.
#'     One of \code{"fixed_rigid"}, \code{"pressure_release"},
#'     \code{"liquid_filled"}, or \code{"gas_filled"}. See the
#'     \link[=boundary_conditions]{boundary conditions documentation} for more
#'     details on these different boundary conditions.}
#'   \item{\code{adaptive}}{A boolean argument controlling whether the PSMS
#'   backend is allowed to stop once the retained modal tail becomes
#'   numerically negligible and to choose an adaptive quadrature order for full
#'   fluid- or gas-filled overlap integrals. When \code{FALSE} (default), the
#'   implementation uses the published hard truncation rule and a fixed
#'   quadrature order unless the user overrides it directly. When
#'   \code{TRUE}, the hard \eqn{m_{\max}} and \eqn{n_{\max}} rules remain the
#'   upper bounds, but the backscatter implementation is allowed to stop early
#'   once the retained tail is both numerically small and gradient-flat.}
#'   \item{\code{simplify_Amn}}{A boolean argument that flags whether or not to
#'   use the simplified calculation for the exansion matrix, \eqn{A_{mn}}, for a
#'   fluid-filled scatterer. See Theory for the mathematical formulation and
#'   assumptions. **Note:** this argument **only** applies to when
#'   \code{boundary = "liquid_filled"} or \code{boundary = "gas_filled"}. It is
#'   otherwise left unused for all other boundary conditions.}
#'   \item{\code{precision}}{A literal argument that allows for levels of
#'   precision when calculating the expansion matrix, \eqn{A_{mn}}. There are
#'   two choices: \code{"double"} for double-precision (64-bit) and
#'   \code{"quad"} for quadruple-precision (128-bit). See Details for more
#'   double- and quadruple-precision usages are implemented.}
#'   \item{\code{n_integration}}{An integer argument that informs the model how
#'   many integration points will be used to calculate \eqn{\alpha_{mn}^{m}},
  #'   which is numerically calculated using Gauss-Legendre quadrature. When left
  #'   as \code{NULL}, the model uses 96 integration points unless
  #'   \code{adaptive = TRUE}, in which case a reduced-frequency-based
  #'   quadrature rule is selected internally for the full fluid- or gas-filled
  #'   solve. See
#'   \link{gauss_legendre} for a full description of how \code{n_integration} is
#'   used. **Note:** this argument **only** applies to when
#'   \code{boundary = "liquid_filled"} or \code{boundary = "gas_filled"}. It is
#'   otherwise left unused for all other boundary conditions.}
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
#' @section Implementation:
#'
#' <u>**\code{C++}**</u>
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
#' <u>**Precision**</u>
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
#' @name PSMS
#' @aliases psms PSMS
#' @docType data
#' @keywords models acoustics
NULL

#' Initialize object for the modal series solution for a prolate spheroid
#' @param object Scatterer-class object.
#' @param frequency Frequency vector (Hz).
#' @param phi_body Body rotation angle (radians).
#' @param boundary Boundary condition.
#' @param adaptive Apply adaptive modal-tail truncation and quadrature selection.
#' @param precision Numerical precision.
  #' @param n_integration Number of integration points. When left as `NULL`, the
  #' model uses 96 integration points unless `adaptive = TRUE`, in which case a
  #' reduced-frequency-based rule is used internally for full fluid- or gas-filled runs.
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
  scatterer_type <- class(object)
  # Detect object shape ========================================================
  scatterer_shape <- acousticTS::extract(object, "shape_parameters")
  # Validate shape +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  if (scatterer_shape$shape != "ProlateSpheroid") {
    stop(
      "The modal series solution for a prolate spheroid requires scatterer to ",
      "be shape-type 'ProlateSpheroid'. Input scatterer is shape-type ",
      paste0("'", scatterer_shape, "'.")
    )
  }
  # Parse body =================================================================
  body <- .hydrate_contrasts(extract(object, "body"), sound_speed_sw, density_sw)
  body_h <- body$h
  body_g <- body$g
  medium_params <- .init_medium_params(sound_speed_sw, density_sw)
  # Define model parameters recipe =============================================
  model_params <- list(
    acoustics = .init_acoustics_df(
      frequency,
      k_sw = sound_speed_sw,
      k_f = body_h * sound_speed_sw
    ),
    precision = precision,
    n_integration = n_integration,
    adaptive = adaptive
  )
  # Determine expansion coefficient Amn method =================================
  # Validate method ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  if (
    !(boundary %in% c(
      "liquid_filled", "gas_filled", "fixed_rigid", "pressure_release"
    )
    )) {
    stop(
      "Only the following values for 'method' are available in this ",
      "implementation of the prolate spheroid modal series solution: ",
      "'liquid_filled' (default), 'gas_filled', 'fixed_rigid',
      'pressure_release'."
    )
  }
  if (! precision %in% c("double", "quad")) {
    stop(
      "'precision' must be either 'double' or 'quad'."
    )
  }
  if (!is.logical(adaptive) || length(adaptive) != 1 || is.na(adaptive)) {
    stop("'adaptive' must be either TRUE or FALSE.")
  }
  # Assign method ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  model_params$Amn_method <- switch(
    boundary,
    liquid_filled = ifelse(simplify_Amn, "Amn_fluid_simplify", "Amn_fluid"),
    gas_filled = ifelse(simplify_Amn, "Amn_fluid_simplify", "Amn_fluid"),
    fixed_rigid = "Amn_fixed_rigid",
    pressure_release = "Amn_pressure_release"
  )
  # Compute body parameters ====================================================
  body_params <- list(
    # Prolate spheroidal coordinate 'xi'++++++++++++++++++++++++++++++++++++++++
    xi = 1 / sqrt(
      1 - (scatterer_shape$radius / (scatterer_shape$length / 2))^2
    ),
    # Roll angle 'phi' (incident direction) ++++++++++++++++++++++++++++++++++++
    phi_body = phi_body,
    # Roll angle 'phi' (scattering direction) ++++++++++++++++++++++++++++++++++
    phi_scatter = phi_body + pi,
    # Theta angle 'theta' (incident direction) +++++++++++++++++++++++++++++++++
    theta_body = body$theta,
    # Theta angle 'theta' (scattering direction) +++++++++++++++++++++++++++++++
    theta_scatter = pi - body$theta,
    # Density ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    density = body_g * density_sw,
    # Sound speed ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    sound_speed = body_h * sound_speed_sw
  )
  # Focal length 'q' +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  body_params$q <- scatterer_shape$length / 2 / body_params$xi
  # Compute the reduced frequency for the surrounding medium ===================
  model_params$acoustics$chi_sw <- model_params$acoustics$k_sw * body_params$q
  # Compute the reduced frequency for the scatterer body =======================
  model_params$acoustics$chi_body <- model_params$acoustics$k_f * body_params$q
  # Determine quadrature order =================================================
  use_adaptive_quadrature <- isTRUE(adaptive) && identical(model_params$Amn_method, "Amn_fluid")
  if (use_adaptive_quadrature) {
    if (!is.null(n_integration)) {
      warning(
        "'n_integration' is ignored when 'adaptive = TRUE' for the full fluid- or gas-filled PSMS solve."
      )
    }
    n_integration <- NULL
  } else if (is.null(n_integration)) {
    n_integration <- 96L
  }
  if (!is.null(n_integration)) {
    if (
      !is.numeric(n_integration) ||
      length(n_integration) != 1 ||
      is.na(n_integration) ||
      n_integration < 1 ||
      n_integration %% 1 != 0
    ) {
      stop("'n_integration' must be a single positive integer.")
    }
    n_integration <- as.integer(n_integration)
  } else {
    n_integration <- NA_integer_
  }
  model_params$n_integration <- n_integration
  # Define limits for 'm' and 'n' iterators ====================================
  model_params$acoustics$m_max <- ceiling(
    2 * model_params$acoustics$k_sw * scatterer_shape$radius
  )
  model_params$acoustics$n_max <- model_params$acoustics$m_max +
    ceiling(0.5 * model_params$acoustics$chi_sw)
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
