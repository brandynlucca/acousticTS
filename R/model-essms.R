################################################################################
# Goodman and Stern (1962) modal series solution for elastic-shelled spheres
################################################################################
#' Elastic-shelled sphere modal series (ESSMS) solution
#'
#' @description
#' Calculates the far-field scattering amplitude and related quantities
#' elastic shelled sphere using the modal series solution, as described by
#' Goodman and Stern (1962).
#'
#' @section Usage:
#' This model is accessed via:
#' \preformatted{
#' target_strength(
#'   ...,
#'   model="essms",
#'   sound_speed_sw,
#'   density_sw,
#'   m_limit
#' )
#' }
#'
#' @section Arguments:
#' \describe{
#'    \item{\code{sound_speed_sw}}{Seawater sound speed (\eqn{m~s^{-1}}).}
#'    \item{\code{density_sw}}{Seawater density (\eqn{kg~m^{-3}}).}
#'    \item{\code{m_limit}}{Optional model truncation limit used to cap the
#'   number of modes in the numerical calculation. By default, this is set to
#'   the maximum value of \eqn{ka + 10}.}
#' }
#'
#' @section Theory:
#' The elastic-shelled sphere model solves the acoustic scattering problem
#' by expanding the incident and scattered fields in terms of spherical Bessel
#' and Hankel functions and Legendre polynomials. The displacement field in the
#' shell is represented as the sum of a scalar potential and a vector potential:
#'
#' \deqn{
#'   \mathbf{u} = \nabla \phi + \nabla \times \mathbf{\Psi}
#' }
#'
#' where \eqn{\phi} and \eqn{\mathbf{\Psi}} satisfy the Helmholtz equations
#' with velocities determined by the elastic properties of the shell:
#'
#' \deqn{
#'   (1/c_L^2) \frac{\partial^2 \phi}{\partial t^2} = \nabla^2 \phi, \qquad
#'   (1/c_T^2) \frac{\partial^2 \mathbf{\Psi}}{\partial t^2} = \nabla^2
#'   \mathbf{\Psi}
#' }
#'
#' with \eqn{c_L = \sqrt{(\lambda + 2\mu)/\rho}} (longitudinal) and
#' \eqn{c_T = \sqrt{\mu/\rho}} (transverse) wave speeds, where \eqn{\lambda}
#' and \eqn{\mu} are the Lamé parameters and \eqn{\rho} is the shell density.
#'
#' The scattered field is expanded as:
#'
#' \deqn{
#'   f_{bs} = -\frac{i}{k} \sum_{m=0}^{\infty} (2m+1) b_m P_m(\cos\theta)
#' }
#'
#' where \eqn{k} is the wavenumber in the surrounding fluid, \eqn{P_m} is the
#' Legendre polynomial of order \eqn{m}, and \eqn{b_m} are the modal
#' coefficients determined by the shell and fluid properties.
#'
#' The modal coefficients \eqn{b_m} are determined by enforcing boundary
#' conditions at the shell-fluid interfaces, resulting in a system of six
#' equations for each mode. These are solved using Cramer's rule as the ratio
#' of two 6 \eqn{\times} 6 determinants.
#'
#' At mode 0:
#'
#' \deqn{
#'   b_m(0) = \frac{
#'     \left|
#'     \begin{array}{ccccc}
#'       a_1 & \alpha_{12} & \alpha_{14} & 0 & 0 \\
#'       a_2 & \alpha_{22} & \alpha_{24} & 0 & 0 \\
#'       0   & \alpha_{42} & \alpha_{44} & \alpha_{46} & 0 \\
#'       0   & \alpha_{52} & \alpha_{54} & \alpha_{56} & 0 \\
#'       0   & \alpha_{62} & \alpha_{64} & \alpha_{66} & 0
#'     \end{array}
#'     \right|
#'   }{
#'     \left|
#'     \begin{array}{ccccc}
#'       \alpha_{11} & \alpha_{12} & \alpha_{14} & 0 & 0 \\
#'       \alpha_{21} & \alpha_{22} & \alpha_{24} & 0 & 0 \\
#'       0   & \alpha_{42} & \alpha_{44} & \alpha_{46} & 0 \\
#'       0   & \alpha_{52} & \alpha_{54} & \alpha_{56} & 0 \\
#'       0   & \alpha_{62} & \alpha_{64} & \alpha_{66} & 0
#'     \end{array}
#'     \right|
#'   }
#' }
#'
#' At modes greater than 0:
#'
#' \deqn{
#'   b_m(m) = -i^m (2m+1) \frac{
#'     \left|
#'     \begin{array}{cccccc}
#'       a_1 & \alpha_{12} & \alpha_{13} & \alpha_{14} & \alpha_{15} & 0 \\
#'       a_2 & \alpha_{22} & \alpha_{23} & \alpha_{24} & \alpha_{25} & 0 \\
#'       0   & \alpha_{32} & \alpha_{33} & \alpha_{34} & \alpha_{35} & 0 \\
#'       0   & \alpha_{42} & \alpha_{43} & \alpha_{44} & \alpha_{45} &
#'       \alpha_{46} \\
#'       0   & \alpha_{52} & \alpha_{53} & \alpha_{54} & \alpha_{55} &
#'       \alpha_{56} \\
#'       0   & \alpha_{62} & \alpha_{63} & \alpha_{64} & \alpha_{65} &
#'       \alpha_{66}
#'     \end{array}
#'     \right|
#'   }{
#'     \left|
#'     \begin{array}{cccccc}
#'       \alpha_{11} & \alpha_{12} & \alpha_{13} & \alpha_{14} & \alpha_{15} &
#'       0 \\
#'       \alpha_{21} & \alpha_{22} & \alpha_{23} & \alpha_{24} & \alpha_{25} &
#'       0 \\
#'       0   & \alpha_{32} & \alpha_{33} & \alpha_{34} & \alpha_{35} & 0 \\
#'       0   & \alpha_{42} & \alpha_{43} & \alpha_{44} & \alpha_{45} &
#'       \alpha_{46} \\
#'       0   & \alpha_{52} & \alpha_{53} & \alpha_{54} & \alpha_{55} &
#'       \alpha_{56} \\
#'       0   & \alpha_{62} & \alpha_{63} & \alpha_{64} & \alpha_{65} &
#'       \alpha_{66}
#'     \end{array}
#'     \right|
#'   }
#' }
#'
#' The elements of these matrices depend on the shell's elastic moduli,
#' thickness, densities, and the acoustic properties of the interior and
#' exterior fluids. The exact values for each element can be calcaulted using
#' equations provided by Goodman and Stern (1962) and Stanton (1990).
#' Since \eqn{\theta = \pi} in the backscattering case, the equation for
#' \eqn{f_{bs}} becomes:
#'
#' \deqn{
#'  f_{bs} = -\frac{i}{k} \sum_{m=0}^{\infty} (2m+1) b_m P_m(\cos\theta) \\
#'  \phantom{f_{bs}} = -\frac{i}{k} \sum_{m=0}^{\infty} (2m+1) b_m P_m(-1) \\
#'  \phantom{f_{bs}} = -\frac{i}{k} \sum_{m=0}^{\infty} (2m+1) b_m (-1)^m
#' }
#'
#' @section Implementation:
#'
#' <u>**\code{C++}**</u>
#'
#' The computation for \eqn{b_m} was done in \code{C++} due to relatively large
#' computational costs with increasing \eqn{ka} (and subsequently larger limits
#' for \eqn{m}) and the efficiency in solving the systems of linear equations
#' containing complex values. Since each matrix is a square, this
#' implementation specifically utilizes lower-upper (LU) decomposition to break
#' down the matrices into the products of the resulting lower and upper
#' triangular matrices. However, this means that the determinant computations
#' are sensitive to ill-conditioned matrices that can amplify numerical errors.
#'
#' There are guards in-place that partially address singularity issues when the
#' denominator is 0. The algorithm does not handle near-singular matrices
#' directly, but it will raise a warning when a matrix is ill-conditioned. This
#' is determined based on the pivot ratio from the calculated condition number.
#' Partial pivoting and row-scaling are also incorporated to improve numerical
#' stability and reduce the effect of high-leverage values in a matrix.
#'
#' <u>**Modal Truncation**</u>
#'
#' The maximum number of terms for \eqn{n} is chosen as \eqn{k_w a_{shell} +
#' 10} (rounded to the nearest integer), which is sufficient for convergence in
#' most practical cases.
#'
#' @seealso
#' See the boundary conditions documentation for
#' more details on how elastic-shelled scattering relates to other boundary
#' conditions,
#' \code{\link{target_strength}}, \code{\link{ESS}}, \code{\link{Sphere}},
#' \code{\link{sphere}}
#'
#' @references
#'
#' Anderson, V.C. (1950). Sound scattering from a fluid sphere. The Journal of
#' The Acoustical Society of America, 22: 426–431.
#'
#' Gaunaurd, G.C., and Wertman, W. (1991). Transient acoustic scattering by
#' fluid-loaded elastic shells. International Journal of Solids and Structures,
#' 27: 699-811.
#'
#' Stanton, T.K. (1990). Sound scattering by spherical and elongated shelled
#' bodies. The Journal of the Acoustical Society of America, 88: 1619-1633.
#'
#' @name ESSMS
#' @aliases essms ESSMS
#' @docType data
#' @keywords models acoustics internal
NULL

#' Initialize ESS-object for modal series solution.
#' @param object ESS-class object.
#' @param frequency Frequency vector (Hz).
#' @param sound_speed_sw Seawater sound speed (m/s).
#' @param density_sw Seawater density (kg/m³).
#' @param m_limit Modal series limit (i.e. max "m").
#' @noRd
essms_initialize <- function(object,
                             frequency,
                             sound_speed_sw = .SEAWATER_SOUND_SPEED_DEFAULT,
                             density_sw = .SEAWATER_DENSITY_DEFAULT,
                             m_limit = NULL) {
  # Detect object class ========================================================
  scatterer_type <- class(object)
  # Parse shell ================================================================
  shell <- extract(object, "shell")
  # Parse fluid ================================================================
  fluid <- extract(object, "fluid")
  # Calculate wavenumbers ======================================================
  k1 <- wavenumber(frequency, sound_speed_sw)
  # Define the maximum Modal series iterator limit =============================
  m_limit <- if (is.null(m_limit)) {
    round(k1 * shell$radius) + 10
  } else {
    m_limit
  }
  # Define model parameters recipe =============================================
  model_params <- list(
    acoustics = transform(
      .init_acoustics_df(frequency, k_sw = sound_speed_sw),
      m_limit = m_limit
    )
  )
  # Initialize lists for fluid and shell =======================================
  # Establish the correct material properties ==================================
  shell_params <- .extract_material_props(shell, sound_speed_sw, density_sw)
  fluid_params <- .extract_material_props(fluid, sound_speed_sw, density_sw)
  # Add morphometrics ==========================================================
  # Shell ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  shell_params$radius <- shell$radius
  shell_params$shell_thickness <- shell$shell_thickness
  # Fluid ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  fluid_params$radius <- fluid$radius
  .init_model_slots(
    object = object,
    model_name = "ESSMS",
    frequency = frequency,
    model_parameters = list(
      parameters = model_params,
      shell = shell_params,
      fluid = fluid_params,
      medium = .init_medium_params(sound_speed_sw, density_sw)
    )
  )
}

#' Calculate ka matrix for Goodman and Stern (1962) model
#' @param frequency Frequency vector
#' @param sound_speed_sw Seawater sound speed
#' @param sound_speed_fluid Fluid sound speed
#' @param sound_speed_longitudinal Longitudinal sound speed in shell
#' @param sound_speed_transversal Transversal sound speed in shell
#' @param radius_shell Shell radius
#' @param radius_fluid Fluid radius
#' @return Matrix of ka values
#' @keywords internal
#' @noRd
.calculate_ka_matrix <- function(frequency, sound_speed_sw,
                                 sound_speed_fluid, sound_speed_longitudinal,
                                 sound_speed_transversal,
                                 radius_shell, radius_fluid) {
  # Calculate layer-specific wave numbers ======================================
  k1 <- wavenumber(frequency, sound_speed_sw)
  k3 <- wavenumber(frequency, sound_speed_fluid)
  kL <- wavenumber(frequency, sound_speed_longitudinal)
  kT <- wavenumber(frequency, sound_speed_transversal)
  # Bind into collective matrix ================================================
  ka_matrix <- rbind(
    k1a_shell = k1 * radius_shell,
    kLa_shell = kL * radius_shell,
    kTa_shell = kT * radius_shell,
    k1a_fluid = k1 * radius_fluid,
    kTa_fluid = kT * radius_fluid,
    kLa_fluid = kL * radius_fluid,
    k3a_fluid = k3 * radius_fluid
  )
  ka_matrix
}

#' Calculates the theoretical TS of an elastic-shelled sphere using the modal
#' series solution from Goodman and Stern (1962).
#'
#' @param object ESS-class object.
#' @noRd
ESSMS <- function(object) {
  # Extract model parameters/inputs ============================================
  model <- extract(object, "model_parameters")$ESSMS
  shell <- model$shell
  fluid <- model$fluid
  medium <- model$medium
  acoustics <- model$parameters$acoustics
  # Get requisite elastic properties ===========================================
  G <- shell$G
  lambda <- shell$lambda
  density_shell <- shell$density
  # Extract the required morphometrics =========================================
  radius_shell <- shell$radius
  radius_fluid <- fluid$radius
  # Get the internal fluid density =============================================
  density_fluid <- fluid$density
  # Calculate shell sound speeds in the longitudinal and transverse directions =
  sound_speed_longitudinal <- sqrt((lambda + 2 * G) / density_shell)
  sound_speed_transversal <- sqrt(G / density_shell)
  # Calculate the associated wavenumbers =======================================
  kL <- wavenumber(acoustics$frequency, sound_speed_longitudinal)
  kT <- wavenumber(acoustics$frequency, sound_speed_transversal)
  # Calculate the reindexed ka values ==========================================
  ka_matrix <- .calculate_ka_matrix(
    acoustics$frequency,
    model$medium$sound_speed,
    fluid$sound_speed,
    sound_speed_longitudinal,
    sound_speed_transversal,
    radius_shell,
    radius_fluid
  )
  # Resolve the boundary conditions ============================================
  Am <- elastic_shell_boundary_conditions(
    ka_matrix,
    acoustics$m_limit,
    lambda,
    G,
    medium$density / density_shell,
    density_fluid / density_shell
  )
  # Compute the linear scattering coefficient ==================================
  m_vec <- 0:max(acoustics$m_limit)
  fbs <- -1i / acoustics$k_sw * rowSums(
    sweep(Am, 2, (-1)^m_vec * (2 * m_vec + 1), FUN = "*"),
    na.rm = TRUE
  )
  # Convert to sigma_bs ========================================================
  sigma_bs <- abs(fbs)^2
  # Define MSS slot for ESS-type scatterer =====================================
  methods::slot(object, "model")$ESSMS <- list(
    frequency = acoustics$frequency,
    ka_shell = acoustics$k_sw * radius_shell,
    ka_fluid = acoustics$k_sw * radius_fluid,
    f_bs = fbs,
    sigma_bs = sigma_bs,
    TS = db(sigma_bs)
  )
  # Return object ==============================================================
  object
}
