################################################################################
# Kirchoff-Ray Mode approximation
################################################################################
#' Kirchhoff Ray Mode (KRM) scattering model
#'
#' @description
#' Computes the far-field acoustic scattering amplitude and derived quantities
#' for fish and similar elongated scatterers using the Kirchhoff–Ray Mode (KRM)
#' approximation described by Clay and Horne (1994). The KRM model is widely
#' used in fisheries acoustics for estimating target strength, particularly
#' for fish with gas-filled swimbladders.
#'
#' The model represents the fish body and swimbladder as a series of contiguous
#' cylindrical elements and computes the coherent sum of their scattered fields.
#' Depending on frequency and element size, scattering from each segment is
#' evaluated using either a low-frequency modal solution or a high-frequency
#' Kirchhoff (ray-based) approximation.
#'
#' @section Usage:
#'
#' This model is accessed via:
#' \preformatted{
#' target_strength(
#'   ...,
#'   model="krm",
#'   sound_speed_sw,
#'   density_sw
#' )
#' }
#'
#' @section Arguments:
#'
#' \describe{
#'   \item{\code{sound_speed_sw}}{Seawater sound speed (\eqn{m~s^{-1}}).}
#'   \item{\code{density_sw}}{Seawater density (\eqn{kg~m^{-3}}).}
#' }
#'
#' @section Scatterer representation:
#' The scatterer must provide body and swimbladder geometry discretized into
#' axial segments. Each segment is described by its longitudinal position and
#' cross-sectional dimensions. Typical geometric quantities include
#' \eqn{x(j)}, \eqn{z_U(j)}, \eqn{z_L(j)}, and radius \eqn{a(j)} for both the
#' body and the swimbladder.
#'
#' The incident plane wave is defined by angle \eqn{\theta}, and coordinates
#' are transformed into a rotated system according to:
#'
#'  \deqn{
#'    u(j) = x(j) \sin \theta - z(j) \cos \theta,
#'  }
#'  \deqn{
#'    v(j) = x(j) \cos \theta + z(j) \sin \theta,
#'  }#'
#'  \deqn{
#'    \Delta u_j = [x(j+1) - x(j)] \sin \theta,
#'  }
#'  \deqn{
#'    a_j = [w(j) + w(j+1)] / 4.
#'  }
#'
#' Segment lengths and effective radii are computed from adjacent points and
#' used to evaluate the scattering contribution from each element.
#'
#' @section Theory:
#'
#' The total backscattering amplitude is computed as the coherent sum of all
#' scatterer segments:
#'
#'  \deqn{
#'    \mathcal{L} =
#'      \sum\limits_{b=1}^{N_b} \mathcal{L}_{B,b} +
#'      \sum\\limits_{s=1}^{N_s} \mathcal{L}_{SB,s},
#'  }
#'
#' where \eqn{\mathcal{S}_{B,b}} is the contribution from the \eqn{b}-th body
#' segment and \eqn{\mathcal{L}_{SB,s}} is from the \eqn{s}-th swimbladder
#' segment. For body segments:
#'
#'  \deqn{
#'    \mathcal{L}_{B,b} \approx -i \frac{(ka_B \sin \theta)^{1/2}}{2\sqrt{\pi}}
#'    \Delta u_b [\mathcal{R}_{wb} e^{-i2kz_U} -
#'    \mathcal{T}_{wb}\mathcal{T}_{bw}
#'    e^{-i2k z_U + i 2 k_B (w_b - v_b - v_L) + i \psi_B}],
#'  }
#'
#' where \eqn{\mathcal{R}_{wb}} and \eqn{\mathcal{T}_{wb}} are reflection and
#' transmission coefficients at the water-body interface, \eqn{k} is the
#' wavenumber in water, and \eqn{k_B} in the body. The upper/lower surface
#' coordinates are \eqn{z_U} and \eqn{z_L}, and \eqn{\psi_B} is an empirical
#' phase correction.
#'
#' The scattering from the swimbladder depends on the dimensionless frequency
#' parameter \eqn{ka_e}, \eqn{a_e} is the equivalent swimbladder radius. For low
#' frequencies (\eqn{ka < 0.15}), the swimbladder is treated as a finite-length
#' gas-filled cylinder, and the scattering is dominated by the first
#' cylindrical mode (breathing mode):
#'
#'  \deqn{
#'    L_M(ka)|_{m=0} = e^{i(\chi - \pi/4)} \frac{L_e}{\pi}
#'    \frac{\sin \Delta}{\Delta} b_0,
#'  }
#'  \deqn{
#'    b_0 = -\frac{1}{1 + i C_0},
#'  }
#'
#' where \eqn{C_0} is a mode coefficient determined by the material properties
#' and  boundary conditions of the swimbladder. For higher frequencies
#' \eqn{ka \ge 0.15}), the Kirchhoff-ray approximation is used:
#'
#' \deqn{
#'  \mathcal{L}_{SB,s} \approx -i \frac{A_{SB,s}
#'  [(k_B a_s + 1)\sin\theta]^{1/2}}{2\sqrt{\pi}} \Delta u_s
#'  \mathcal{R}_{bc} \mathcal{T}_{wb} e^{-i (2 k_B v_s + \psi_p)}.
#' }
#'
#' Here, \eqn{a_s} and \eqn{v_s} are the averaged radius and longitudinal
#' position for  the segment, \eqn{\mathcal{R}_{bc}} is the reflection
#' coefficient at the body-cylinder interface, and \eqn{A_{SB,s}} and
#' \eqn{\psi_p)} are empirical amplitude and phase adjustments, respectively.
#'
#' This hybrid approach ensures that the KRM model correctly captures the
#' low-frequency  resonance of the swimbladder while maintaining high-frequency
#' accuracy via the Kirchhoff-ray mode approximation.
#'
#' \strong{Assumptions and limitations}
#'
#' The KRM formulation assumes that the scatterer is smooth, elongated, and
#' approximately axisymmetric, such that both the Kirchhoff (ray) approximation
#' and low-order cylindrical mode solutions are valid representations of the
#' local scattering physics. The body andswimbladder are discretized into short
#' axial segments, each with its own local height and width. While the segments
#' are treated as locally cylindrical for the purpose of computing the
#' scattered field, their dimensions may vary independently along the body
#' axis. The total scattered field is obtained as a coherent sum of the
#' contributions from these segments, which implicitly assumes that geometric
#' and material properties vary gradually relative to the acoustic wavelength.
#'
#' The model further assumes that scattering is dominated by first-order
#' interactions between the incident field and each local element. Multiple
#' internal reflections and higher-order multiple scattering between
#' different segments are neglected. For the fish body, the material is
#' treated as fluid-like, with no explicit treatment of elastic shear waves
#' or bending modes. As a result, elastic shell effects, bone resonance, and
#' complex internal structural scattering are not represented explicitly,
#' except insofar as they may be approximated through effective reflection
#' or transmission coefficients.
#'
#' The Kirchhoff approximation employed in the high-frequency regime assumes
#' that the local radius of curvature of each segment is large compared to
#' the acoustic wavelength, and that the incident wave interacts primarily
#' with the illuminated surface. Accuracy degrades for end-on incidence,
#' sharp geometric discontinuities, or highly concave regions, where shadowing
#' and diffraction effects become significant. The low-frequency mode
#' solution used for small values of \eqn{ka} is limited to the lowest-order
#' cylindrical modes and does not capture higher-order resonances that may
#' arise for more complex internal geometries.
#'
#' Orientation dependence is treated deterministically through the incident
#' angle \eqn{\theta}, and the model does not account for stochastic body
#' deformations, posture changes, or dynamic swimbladder shape variations.
#' Surface roughness may be included through empirical attenuation of the
#' reflection coefficient, but this treatment is phenomenological and does
#' not represent scattering from discrete roughness elements. Consequently,
#' the KRM is best suited for predicting mean or orientation-specific target
#' strength for smooth-bodied organisms, and its accuracy decreases when
#' applied to organisms with highly irregular shapes, strong elastic
#' contrasts, or complex internal skeletal structures.
#'
#' @section Fluid-only and alternative scattering configurations:
#'
#' Although the KRM formulation is most commonly applied to fish with
#' gas-filled swimbladders, the model can also be used to compute target
#' strength using only the fluid-like scattering contribution of the body.
#'
#' In this configuration, the swimbladder scattering term is omitted and the
#' total scattered field is obtained solely from the fluid contrast between
#' the animal body and the surrounding medium. This mode of operation is
#' appropriate for organisms that lack swimbladders (e.g. many invertebrates,
#' elasmobranchs, or larval fish), or when swimbladder geometry is unknown or
#' intentionally excluded.
#'
#' More generally, the swimbladder component in the KRM framework may be
#' replaced or supplemented by alternative internal scattering features,
#' provided they can be represented geometrically and assigned appropriate
#' acoustic boundary conditions. Examples include rigid or elastic skeletal
#' structures (e.g. vertebral columns) or other localized impedance contrasts
#' within an otherwise fluid-like body.
#'
#' In all cases, the total target strength is computed from the coherent sum
#' of the selected scattering components. Users are responsible for ensuring
#' that the assumed boundary conditions and material contrasts are consistent
#' with the biological structure being modeled.
#'
#' @section Implementation:
#'
#' The implementation extracts geometric and acoustic parameters from the
#' input object, transforms coordinates to the rotated reference frame,
#' evaluates the appropriate modal or ray-based scattering expression for each
#' segment, and coherently sums the contributions to obtain the total
#' scattering amplitude and target strength.
#'
#' @seealso
#' See the \link[=boundary_conditions]{boundary conditions documentation} for
#' more details on fluid-like scattering assumptions,
#' \code{\link{target_strength}}, \code{\link{SBF}}, \code{\link{FLS}}
#'
#' @references
#' Clay, C.S. (1991). Low resolution acoustic scattering models: Fluid-filled
#' cylinders and fish with swimbladders. The Journal of the Acoustical Society
#' of America, 89: 2168-2179.
#'
#' Clay, C.S. (1992). Composite ray-mode approximations for backscattered sound
#' from gas-filled cylinders and swimbladders.The Journal of the Acoustical
#' Society of America, 92: 2173-2180.
#'
#' Clay, C.S., and Horne, J.K. (1994). Acoustic models of fish: The Atlantic cod
#' (<i>Gadus morhua</i>). The Journal of the Acoustical Society of America, 96:
#' 1661-1668.
#'
#' @name KRM
#' @aliases krm KRM
#' @docType data
#' @keywords models acoustics
#' @export
NULL

#' Initialize SBF-class object for KRM calculations.
#' @param object SBF-class object
#' @param frequency Frequency (Hz).
#' @param sound_speed_sw Seawater sound speed (m/s).
#' @param density_sw Seawater density (kg/m³).
#' @noRd
krm_initialize <- function(object,
                           frequency,
                           sound_speed_sw = 1500,
                           density_sw = 1026) {
  # Detect scatterer object class ==============================================
  scatterer_type <- class(object)
  # Parse body =================================================================
  body <- acousticTS::extract(object, "body")
  # Parse shape ================================================================
  shape <- acousticTS::extract(object, "shape_parameters")
  # Helper to complete contrast/absolute pairs =================================
  fill_props <- function(comp, medium_sound_speed, medium_density) {
    # derive contrasts if missing
    if (is.null(comp$g) && !is.null(comp$density)) {
      comp$g <- comp$density / medium_density
    }
    if (is.null(comp$h) && !is.null(comp$sound_speed)) {
      comp$h <- comp$sound_speed / medium_sound_speed
    }
    # derive absolutes if missing
    if (is.null(comp$density) && !is.null(comp$g)) {
      comp$density <- comp$g * medium_density
    }
    if (is.null(comp$sound_speed) && !is.null(comp$h)) {
      comp$sound_speed <- comp$h * medium_sound_speed
    }
    comp
  }
  # Derive contrasts from absolute properties when needed =====================
  contrasts <- .derive_contrasts(body, sound_speed_sw, density_sw)
  body$h <- if (!is.null(body$h)) body$h else contrasts$h
  body$g <- if (!is.null(body$g)) body$g else contrasts$g
  body <- fill_props(body, sound_speed_sw, density_sw)
  # Define medium parameters ===================================================
  medium_params <- data.frame(
    sound_speed = sound_speed_sw,
    density = density_sw
  )
  # Define model parameters recipe =============================================
  model_params <- list(
    acoustics = data.frame(
      # Frequency (Hz) =========================================================
      frequency = frequency,
      # Wavenumber (medium) ====================================================
      k_sw = acousticTS::wavenumber(
        frequency,
        sound_speed_sw
      )
    )
  )
  # Class-specific model parameters ============================================
  if (scatterer_type == "FLS") {
    # Define body parameters recipe ============================================
    body_params <- data.frame(
      length = shape$length,
      theta = body$theta,
      density = density_sw * body$g,
      sound_speed = sound_speed_sw * body$h
    )
    # Body wave number =========================================================
    model_params$acoustics$k_b <- acousticTS::wavenumber(
      frequency,
      sound_speed_sw * body$h
    )
    # Define body segment count ================================================
    model_params$ns_b <- ncol(body$rpos)
    # Slot new model parameters input for KRM ==================================
    methods::slot(object, "model_parameters")$KRM <- list(
      parameters = model_params,
      medium = medium_params,
      body = body_params
    )
  } else if (scatterer_type == "SBF") {
    # Define body parameters recipe ============================================
    body_params <- data.frame(
      length = shape$body$length,
      theta = body$theta,
      density = body$density,
      sound_speed = body$sound_speed
    )
    # Pull bladder information =================================================
    bladder <- acousticTS::extract(object, "bladder")
    bladder <- fill_props(bladder, sound_speed_sw, density_sw)
    body <- fill_props(body, sound_speed_sw, density_sw)
    # Define bladder parameters recipe =========================================
    bladder_params <- data.frame(
      length = shape$bladder$length,
      theta = bladder$theta,
      density = bladder$density,
      sound_speed = bladder$sound_speed
    )
    # Body wave number =========================================================
    model_params$acoustics$k_b <- acousticTS::wavenumber(
      frequency,
      body$sound_speed
    )
    # Define body segment count ================================================
    model_params$ns_b <- shape$body$n_segments
    # Define bladder segment count =============================================
    model_params$ns_sb <- shape$bladder$n_segments
    # Slot new model parameters input for KRM ==================================
    methods::slot(object, "model_parameters")$KRM <- list(
      parameters = model_params,
      medium = medium_params,
      body = body_params,
      bladder = bladder_params
    )
  }
  # Slot new model output for KRM ==============================================
  methods::slot(object, "model")$KRM <- data.frame(frequency,
                                                   sigma_bs = rep(
                                                     NA,
                                                     length(frequency)
                                                   )
  )
  # Output =====================================================================
  return(object)
}

#' Calculates the theoretical TS using Kirchoff-ray Mode approximation.
#' @noRd
KRM <- function(object) {
  # Detect object class ========================================================
  scatterer_type <- class(object)
  # Extract model parameter inputs =============================================
  model <- extract(object, "model_parameters")$KRM
  # Extract body parameters ====================================================
  body <- extract(object, "body")
  # Calculate reflection coefficient for medium-body interface =================
  R12 <- reflection_coefficient(model$medium, model$body)
  # Calculate transmission coefficient and its reverse =========================
  T12T21 <- 1 - R12 * R12
  # Sum across body position vector ============================================
  rpos <- switch(scatterer_type,
                 FLS = rbind(
                   x = body$rpos[1, ],
                   w = c(
                     body$radius[2],
                     body$radius[2:(length(body$radius) - 1)],
                     body$radius[(length(body$radius)) - 1]
                   ) * 2,
                   zU = c(
                     body$radius[2],
                     body$radius[2:(length(body$radius) - 1)],
                     body$radius[(length(body$radius)) - 1]
                   ),
                   zL = -c(
                     body$radius[2],
                     body$radius[2:(length(body$radius) - 1)],
                     body$radius[(length(body$radius)) - 1]
                   )
                 ),
                 SBF = body$rpos
  )
  body_rpos_sum <- along_sum(rpos, model$parameters$ns_b)
  # Approximate radius of body cylinders =======================================
  a_body <- switch(scatterer_type,
                   FLS = body_rpos_sum[2, ] / 4,
                   SBF = body_rpos_sum[2, ] / 4
  )
  # Combine wavenumber (k) and radii to calculate "ka" =========================
  ka_body <- matrix(
    data = rep(a_body,
               each = length(model$parameters$acoustics$k_sw)
    ),
    ncol = length(a_body),
    nrow = length(model$parameters$acoustics$k_b)
  ) * model$parameters$acoustics$k_b
  # Convert c-z coordinates to required u-v rotated coordinates ================
  uv_body <- acousticTS::body_rotation(
    body_rpos_sum,
    body$rpos,
    body$theta,
    length(model$parameters$acoustics$k_sw)
  )
  # Calculate body empirical phase shift function ==============================
  body_dorsal_sum <- switch(scatterer_type,
                            FLS = matrix(
                              data = rep(
                                body_rpos_sum[3, ],
                                each = length(model$parameters$acoustics$k_sw)
                              ),
                              ncol = length(body_rpos_sum[3, ]),
                              nrow = length(model$parameters$acoustics$k_sw)
                            ) / 2,
                            SBF = matrix(
                              data = rep(
                                body_rpos_sum[3, ],
                                each = length(model$parameters$acoustics$k_sw)
                              ),
                              ncol = length(body_rpos_sum[3, ]),
                              nrow = length(model$parameters$acoustics$k_sw)
                            ) / 2
  )
  Psi_b <- -pi * model$parameters$acoustics$k_b * body_dorsal_sum /
    (2 * (model$parameters$acoustics$k_b * body_dorsal_sum + 0.4))
  # Estimate natural log function (phase, etc.) ================================
  exp_body <- exp(-2i * model$parameters$acoustics$k_sw * uv_body$vbU) -
    T12T21 * exp(-2i * model$parameters$acoustics$k_sw * uv_body$vbU +
                   2i * model$parameters$acoustics$k_b *
                   (uv_body$vbU - uv_body$vbL) + 1i * Psi_b)
  # Resolve summation term =====================================================
  body_summation <- sqrt(ka_body) * uv_body$delta_u
  # Calculate linear scattering length (m) =====================================
  f_body <- rowSums(-((1i * (R12 / (2 * sqrt(pi)))) *
                        body_summation * exp_body))
  if (scatterer_type == "FLS") {
    # Define KRM slot for FLS-type scatterer ===================================
    methods::slot(object, "model")$KRM <- data.frame(
      frequency = model$parameters$acoustics$frequency,
      ka = model$parameters$acoustics$k_sw *
        stats::median(a_body, na.rm = TRUE),
      f_bs = f_body,
      sigma_bs = abs(f_body) * abs(f_body),
      TS = 20 * log10(abs(f_body))
    )
  } else if (scatterer_type == "SBF") {
    #### Repeat process for bladder ============================================
    # Extract bladder parameters ===============================================
    bladder <- acousticTS::extract(object, "bladder")
    # Calculate reflection coefficient for bladder =============================
    R23 <- acousticTS::reflection_coefficient(
      body,
      bladder
    )
    # Sum across body/swimbladder position vectors =============================
    bladder_rpos_sum <- acousticTS::along_sum(
      bladder$rpos,
      model$parameters$ns_sb
    )
    # Approximate radii of bladder discrete cylinders ==========================
    a_bladder <- bladder_rpos_sum[2, ] / 4
    # Combine wavenumber (k) and radii to calculate "ka" for bladder ===========
    ka_bladder <- matrix(
      data = rep(a_bladder, each = length(model$parameters$acoustics$k_sw)),
      ncol = length(a_bladder),
      nrow = length(model$parameters$acoustics$k_sw)
    ) * model$parameters$acoustics$k_sw
    # Calculate Kirchoff approximation empirical factor, A_sb ==================
    A_sb <- ka_bladder / (ka_bladder + 0.083)
    # Calculate empirical phase shift for a fluid cylinder, Psi_p ==============
    Psi_p <- ka_bladder / (40 + ka_bladder) - 1.05
    # Convert x-z coordinates to requisite u-v rotated coordinates =============
    uv_bladder <- acousticTS::bladder_rotation(
      bladder_rpos_sum,
      bladder$rpos,
      bladder$theta,
      length(model$parameters$acoustics$k_sw)
    )
    # Estimate natural log functions  ==========================================
    exp_bladder <- exp(-1i * (2 * model$parameters$acoustics$k_b *
                                uv_bladder$v + Psi_p)) * uv_bladder$delta_u
    # Calculate the summation term =============================================
    bladder_summation <- A_sb * sqrt((ka_bladder + 1) * sin(bladder$theta))
    # Estimate backscattering length, f_fluid/f_soft ===========================
    f_bladder <- rowSums(-1i * (R23 * T12T21) / (2 * sqrt(pi)) *
                           bladder_summation * exp_bladder)
    # Estimate total backscattering length, f_bs ===============================
    f_bs <- f_body + f_bladder
    # Define KRM slot for FLS-type scatterer ===================================
    methods::slot(object, "model")$KRM <- data.frame(
      frequency = model$parameters$acoustics$frequency,
      ka = model$parameters$acoustics$k_sw *
        stats::median(a_body, na.rm = TRUE),
      f_body = f_body,
      f_bladder = f_bladder,
      f_bs = f_bs,
      sigma_bs = abs(f_bs) * abs(f_bs),
      TS = 20 * log10(abs(f_bs))
    )
  }
  # Return object ==============================================================
  object
}
