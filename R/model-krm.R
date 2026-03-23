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
#'   density_sw,
#'   krm_variant = "lowcontrast"
#' )
#' }
#'
#' @section Arguments:
#'
#' \describe{
#'   \item{\code{sound_speed_sw}}{Seawater sound speed (\eqn{m~s^{-1}}).}
#'   \item{\code{density_sw}}{Seawater density (\eqn{kg~m^{-3}}).}
#'   \item{\code{krm_variant}}{Swimbladder-medium convention for combined
#'   body-plus-bladder targets. Use \code{"lowcontrast"} for the low-contrast
#'   approximation \eqn{k_B \approx k} in both swimbladder regimes,
#'   \code{"body_embedded"} for the literal body-embedded interpretation of
#'   Clay and Horne (1994), or \code{"mixed"} for the intermediate convention
#'   in which the high-\eqn{ka} swimbladder term uses the body medium but the
#'   low-\eqn{ka} breathing-mode term uses the low-contrast approximation.
#'   This argument is ignored for fluid-only body targets.}
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
#'  }
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
#'      \sum\limits_{s=1}^{N_s} \mathcal{L}_{SB,s},
#'  }
#'
#' where \eqn{\mathcal{L}_{B,b}} is the contribution from the \eqn{b}-th body
#' segment and \eqn{\mathcal{L}_{SB,s}} is from the \eqn{s}-th swimbladder
#' segment. For body segments:
#'
#'  \deqn{
#'    \mathcal{L}_{B,b} \approx -i \frac{\mathcal{R}_{wb}}{2\sqrt{\pi}}
#'    (k a_b)^{1/2} \Delta u_b \left[
#'      e^{-i2k v_{U,b}} -
#'      \mathcal{T}_{wb}\mathcal{T}_{bw}
#'      e^{-i2k v_{U,b} + i 2 k_B (v_{U,b} - v_{L,b}) + i \psi_{B,b}}
#'    \right],
#'  }
#'
#' where \eqn{\mathcal{R}_{wb}} and \eqn{\mathcal{T}_{wb}} are reflection and
#' transmission coefficients at the water-body interface, \eqn{k} is the
#' wavenumber in water, and \eqn{k_B} in the body. The quantities
#' \eqn{v_{U,b}} and \eqn{v_{L,b}} are the rotated upper- and lower-surface
#' coordinates of the \eqn{b}-th segment, \eqn{\Delta u_b = [x(b+1)-x(b)]
#' \sin\theta}, and \eqn{\psi_{B,b}} is an empirical phase correction.
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
#'  \mathcal{L}_{SB,s} \approx -i \frac{\mathcal{R}_{bc}
#'  \mathcal{T}_{wb}\mathcal{T}_{bw}}{2\sqrt{\pi}} A_{SB,s}
#'  [(k a_s + 1)\sin\theta]^{1/2} \Delta u_s
#'  e^{-i (2 k_B v_s + \psi_{p,s})}.
#' }
#'
#' Here, \eqn{a_s} and \eqn{v_s} are the averaged radius and longitudinal
#' position for  the segment, \eqn{\mathcal{R}_{bc}} is the reflection
#' coefficient at the body-cylinder interface, and \eqn{A_{SB,s}} and
#' \eqn{\psi_{p,s}} are empirical amplitude and phase adjustments,
#' respectively.
#'
#' For body-plus-swimbladder targets, the default coherent KRM combines the
#' two complex scattering lengths as:
#'
#' \deqn{
#'   f_{bs}^{(\mathrm{coh})} = f_{body} + f_{bladder},
#' }
#'
#' which gives:
#'
#' \deqn{
#'   \sigma_{bs}^{(\mathrm{coh})} =
#'   \left|f_{body} + f_{bladder}\right|^2.
#' }
#'
#' @section Swimbladder medium conventions:
#' Clay and Horne (1994) derive the swimbladder term for a gas-filled inclusion
#' embedded in body tissue, but also note that the low body-water contrast can
#' justify the approximation \eqn{k_B \approx k}. For combined body-plus-
#' bladder targets, `acousticTS` exposes that choice through
#' \code{krm_variant}:
#'
#' \describe{
#'   \item{\code{"body_embedded"}}{Use the body medium for both the
#'   high-\eqn{ka} Kirchhoff swimbladder term and the low-\eqn{ka}
#'   breathing-mode term. This is the most literal interpretation of the
#'   body-embedded swimbladder geometry in Clay and Horne (1994).}
#'   \item{\code{"mixed"}}{Use the body medium for the high-\eqn{ka}
#'   swimbladder term but the water approximation for the low-\eqn{ka}
#'   breathing-mode term. This follows the later mixed convention discussed in
#'   the fisheries-acoustics literature.}
#'   \item{\code{"lowcontrast"}}{Use the low-contrast approximation
#'   \eqn{k_B \approx k} for both swimbladder regimes, i.e. evaluate both the
#'   high-\eqn{ka} and low-\eqn{ka} swimbladder terms with the external-medium
#'   wavenumber rather than the body-medium wavenumber.}
#' }
#'
#' These named variants make the swimbladder-medium assumption explicit while
#' keeping the public API tied to the scientific interpretation rather than to
#' implementation-specific knob names.
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
#' @param krm_variant Swimbladder-medium convention for combined body-plus-
#' bladder targets. One of \code{"lowcontrast"}, \code{"mixed"}, or
#' \code{"body_embedded"}.
#' @noRd
krm_initialize <- function(object,
                           frequency,
                           sound_speed_sw = 1500,
                           density_sw = 1026,
                           krm_variant = c("lowcontrast",
                                           "mixed",
                                           "body_embedded")) {
  krm_variant <- .krm_variant_settings(krm_variant)
  # Detect scatterer object class ==============================================
  scatterer_type <- class(object)
  # Parse body =================================================================
  body <- acousticTS::extract(object, "body")
  # Parse shape ================================================================
  shape <- acousticTS::extract(object, "shape_parameters")
  # Derive contrasts from absolute properties when needed =====================
  contrasts <- .derive_contrasts(body, sound_speed_sw, density_sw)
  body$h <- if (!is.null(body$h)) body$h else contrasts$h
  body$g <- if (!is.null(body$g)) body$g else contrasts$g
  body <- .complete_material_props(body, sound_speed_sw, density_sw)
  # Define medium parameters ===================================================
  medium_params <- .init_medium_params(sound_speed_sw, density_sw)
  # Define model parameters recipe =============================================
  model_params <- list(
    acoustics = .init_acoustics_df(
      frequency,
      k_sw = sound_speed_sw
    )
  )
  # Class-specific model parameters ============================================
  if (scatterer_type == "FLS") {
    model_parameters <- .krm_initialize_fls(
      body = body,
      shape = shape,
      model_params = model_params,
      medium_params = medium_params,
      sound_speed_sw = sound_speed_sw,
      density_sw = density_sw,
      krm_variant = krm_variant
    )
  } else if (scatterer_type == "SBF") {
    model_parameters <- .krm_initialize_sbf(
      object = object,
      body = body,
      shape = shape,
      model_params = model_params,
      medium_params = medium_params,
      sound_speed_sw = sound_speed_sw,
      density_sw = density_sw,
      krm_variant = krm_variant
    )
  }

  .init_model_slots(
    object = object,
    model_name = "KRM",
    frequency = frequency,
    model_parameters = model_parameters
  )
}
################################################################################
.krm_initialize_fls <- function(body,
                                shape,
                                model_params,
                                medium_params,
                                sound_speed_sw,
                                density_sw,
                                krm_variant) {
  body_params <- data.frame(
    length = shape$length,
    theta = body$theta,
    density = density_sw * body$g,
    sound_speed = sound_speed_sw * body$h
  )

  model_params$acoustics$k_b <- wavenumber(frequency = model_params$acoustics$frequency,
                                           sound_speed = sound_speed_sw * body$h)
  model_params$ns_b <- ncol(body$rpos)

  list(
    parameters = model_params,
    medium = medium_params,
    body = body_params,
    variant = krm_variant
  )
}

################################################################################
.krm_initialize_sbf <- function(object,
                                body,
                                shape,
                                model_params,
                                medium_params,
                                sound_speed_sw,
                                density_sw,
                                krm_variant) {
  bladder <- acousticTS::extract(object, "bladder")
  bladder <- .complete_material_props(bladder, sound_speed_sw, density_sw)
  body <- .complete_material_props(body, sound_speed_sw, density_sw)

  body_params <- data.frame(
    length = shape$body$length,
    theta = body$theta,
    density = body$density,
    sound_speed = body$sound_speed
  )
  bladder_params <- data.frame(
    length = shape$bladder$length,
    theta = bladder$theta,
    density = bladder$density,
    sound_speed = bladder$sound_speed
  )

  model_params$acoustics$k_b <- wavenumber(
    frequency = model_params$acoustics$frequency,
    sound_speed = body$sound_speed
  )
  model_params$acoustics$k_sb <- wavenumber(
    frequency = model_params$acoustics$frequency,
    sound_speed = bladder$sound_speed
  )
  model_params$ns_b <- ncol(body$rpos)
  model_params$ns_sb <- ncol(bladder$rpos)

  list(
    parameters = model_params,
    medium = medium_params,
    body = body_params,
    bladder = bladder_params,
    variant = krm_variant
  )
}

################################################################################
.krm_variant_settings <- function(krm_variant = c("lowcontrast",
                                                  "mixed",
                                                  "body_embedded")) {
  krm_variant <- match.arg(krm_variant)

  switch(
    krm_variant,
    lowcontrast = list(
      name = "lowcontrast",
      high_ka_medium = "water",
      low_ka_medium = "water",
      description = paste(
        "Low-contrast approximation for both swimbladder regimes,",
        "i.e. use the external-medium wavenumber in both swimbladder terms."
      )
    ),
    mixed = list(
      name = "mixed",
      high_ka_medium = "body",
      low_ka_medium = "water",
      description = paste(
        "Body medium for the high-ka swimbladder term and water",
        "approximation for the low-ka breathing-mode term."
      )
    ),
    body_embedded = list(
      name = "body_embedded",
      high_ka_medium = "body",
      low_ka_medium = "body",
      description = paste(
        "Literal body-embedded swimbladder convention after",
        "Clay and Horne (1994)."
      )
    )
  )
}
################################################################################
.krm_bladder_geometry <- function(rpos, theta) {
  # Parse bladder coordinate vectors ==========================================
  x <- rpos[1, ]
  w <- rpos[2, ]
  z_u <- rpos[3, ]
  z_l <- rpos[4, ]
  # Resolve local axial spacing ===============================================
  dx <- diff(x)
  # Summarize adjacent-cylinder geometry ======================================
  sum_rpos <- along_sum(rpos, ncol(rpos))
  center_z <- (sum_rpos[3, ] + sum_rpos[4, ]) / 4
  v_mid <- sum_rpos[1, ] * cos(theta) / 2 + center_z * sin(theta)
  # Approximate local semi-axes at each segment end ===========================
  half_width0 <- w[-length(w)] / 2
  half_width1 <- w[-1] / 2
  half_height0 <- (z_u[-length(z_u)] - z_l[-length(z_l)]) / 2
  half_height1 <- (z_u[-1] - z_l[-1]) / 2
  # Estimate linear taper across each segment =================================
  dx_safe <- ifelse(abs(dx) < sqrt(.Machine$double.eps), NA_real_, dx)
  width_slope <- (half_width1 - half_width0) / dx_safe
  height_slope <- (half_height1 - half_height0) / dx_safe
  # Integrate local elliptical-cylinder volume =================================
  dv <- pi * abs(
    half_height0 * half_width0 * dx +
      (height_slope * half_width0 + width_slope * half_height0) * dx^2 / 2 +
      height_slope * width_slope * dx^3 / 3
  )
  # Convert local volume to equivalent cylinder radius ========================
  a_eq <- sqrt(dv / (pi * abs(dx)))
  a_eq_fallback <- sqrt(pmax(
    ((half_width0 + half_width1) / 2) * ((half_height0 + half_height1) / 2),
    0
  ))
  # Replace degenerate segments with averaged cross-sectional radii ============
  bad_idx <- !is.finite(a_eq)
  if (any(bad_idx)) {
    a_eq[bad_idx] <- a_eq_fallback[bad_idx]
  }
  # Resolve whole-bladder equivalent-cylinder dimensions ======================
  length_eq <- max(x, na.rm = TRUE) - min(x, na.rm = TRUE)
  a_eq_total <- sqrt(sum(dv, na.rm = TRUE) / (pi * length_eq))
  # Return geometry bookkeeping ===============================================
  list(
    sum_rpos = sum_rpos,
    delta_x = dx,
    v_mid = v_mid,
    a_eq = a_eq,
    length_eq = length_eq,
    a_eq_total = a_eq_total
  )
}
################################################################################
.krm_low_mode_b0 <- function(k_medium, k_bladder, a_eq, g13, h13) {
  # Combine acoustic size terms for exterior and bladder media ================
  k2a <- outer(k_medium, a_eq)
  k3a <- outer(k_bladder, a_eq)
  # Evaluate cylindrical Bessel terms entering the m = 0 mode coefficient =====
  J0_2 <- jc(0, k2a)
  Y0_2 <- yc(0, k2a)
  J0p_2 <- jcd(0, k2a)
  Y0p_2 <- ycd(0, k2a)
  J0_3 <- jc(0, k3a)
  J0p_3 <- jcd(0, k3a)
  # Resolve the modal boundary-condition coefficient C0 =======================
  ratio_3 <- J0p_3 / J0_3
  C0 <- (
    ratio_3 * Y0_2 - g13 * h13 * Y0p_2
  ) / (
    ratio_3 * J0_2 - g13 * h13 * J0p_2
  )
  # Convert C0 to the breathing-mode scattering coefficient ===================
  -1 / (1 + 1i * C0)
}

################################################################################
.krm_body_rpos <- function(body, scatterer_type) {
  switch(scatterer_type,
    FLS = if ("w" %in% rownames(body$rpos)) {
      body$rpos
    } else {
      radius_profile <- .dwba_body_radius(body)
      rbind(
        x = body$rpos[1, ],
        w = radius_profile * 2,
        zU = radius_profile,
        zL = -radius_profile
      )
    },
    SBF = body$rpos
  )
}

################################################################################
.krm_body_term <- function(model, body, scatterer_type) {
  rpos <- .krm_body_rpos(body, scatterer_type)
  body_rpos_sum <- along_sum(rpos, model$parameters$ns_b)
  a_body <- body_rpos_sum[2, ] / 4

  ka_body <- matrix(
    data = rep(a_body, each = length(model$parameters$acoustics$k_sw)),
    ncol = length(a_body),
    nrow = length(model$parameters$acoustics$k_sw)
  ) * model$parameters$acoustics$k_sw

  uv_body <- body_rotation(
    body_rpos_sum,
    body$rpos,
    body$theta,
    length(model$parameters$acoustics$k_sw)
  )

  body_dorsal_sum <- matrix(
    data = rep(body_rpos_sum[3, ], each = length(model$parameters$acoustics$k_sw)),
    ncol = length(body_rpos_sum[3, ]),
    nrow = length(model$parameters$acoustics$k_sw)
  ) / 2
  Psi_b <- -pi * model$parameters$acoustics$k_b * body_dorsal_sum /
    (2 * (model$parameters$acoustics$k_b * body_dorsal_sum + 0.4))

  R12 <- reflection_coefficient(model$medium, model$body)
  T12T21 <- 1 - R12 * R12

  exp_body <- exp(-2i * model$parameters$acoustics$k_sw * uv_body$vbU) -
    T12T21 * exp(
      -2i * model$parameters$acoustics$k_sw * uv_body$vbU +
        2i * model$parameters$acoustics$k_b * (uv_body$vbU - uv_body$vbL) +
        1i * Psi_b
    )

  delta_u_body <- matrix(
    data = rep(uv_body$delta_u, each = length(model$parameters$acoustics$k_sw)),
    ncol = length(uv_body$delta_u),
    nrow = length(model$parameters$acoustics$k_sw)
  )
  body_summation <- sqrt(ka_body) * delta_u_body

  list(
    f_body = rowSums(
      -((1i * (R12 / (2 * sqrt(pi)))) * body_summation * exp_body)
    ),
    a_body = a_body,
    R12 = R12,
    T12T21 = T12T21
  )
}

################################################################################
.krm_variant_media <- function(model, krm_variant) {
  if (is.null(krm_variant)) {
    krm_variant <- .krm_variant_settings("lowcontrast")
  }

  list(
    k_bladder_high = if (identical(krm_variant$high_ka_medium, "body")) {
      model$parameters$acoustics$k_b
    } else {
      model$parameters$acoustics$k_sw
    },
    k_bladder_low = if (identical(krm_variant$low_ka_medium, "body")) {
      model$parameters$acoustics$k_b
    } else {
      model$parameters$acoustics$k_sw
    },
    low_density = if (identical(krm_variant$low_ka_medium, "body")) {
      model$body$density
    } else {
      model$medium$density
    },
    low_sound_speed = if (identical(krm_variant$low_ka_medium, "body")) {
      model$body$sound_speed
    } else {
      model$medium$sound_speed
    }
  )
}

################################################################################
.krm_swimbladder_high_term <- function(model,
                                       bladder,
                                       bladder_geom,
                                       R23,
                                       T12T21,
                                       k_bladder_high) {
  bladder_rpos_sum <- bladder_geom$sum_rpos
  a_bladder <- bladder_rpos_sum[2, ] / 4

  ka_bladder_empirical <- matrix(
    data = rep(a_bladder, each = length(model$parameters$acoustics$k_sw)),
    ncol = length(a_bladder),
    nrow = length(model$parameters$acoustics$k_sw)
  ) * model$parameters$acoustics$k_sw
  ka_bladder <- matrix(
    data = rep(a_bladder, each = length(model$parameters$acoustics$k_sw)),
    ncol = length(a_bladder),
    nrow = length(model$parameters$acoustics$k_sw)
  ) * k_bladder_high

  A_sb <- ka_bladder_empirical / (ka_bladder_empirical + 0.083)
  Psi_p <- ka_bladder_empirical / (40 + ka_bladder_empirical) - 1.05

  uv_bladder <- bladder_rotation(
    bladder_rpos_sum,
    bladder$rpos,
    bladder$theta,
    length(model$parameters$acoustics$k_sw)
  )
  delta_u_bladder <- matrix(
    data = rep(uv_bladder$delta_u, each = length(model$parameters$acoustics$k_sw)),
    ncol = length(uv_bladder$delta_u),
    nrow = length(model$parameters$acoustics$k_sw)
  )
  exp_bladder <- exp(-1i * (2 * k_bladder_high * uv_bladder$v + Psi_p)) *
    delta_u_bladder
  bladder_summation <- A_sb * sqrt((ka_bladder + 1) * sin(bladder$theta))

  list(
    f = rowSums(
      -1i * (R23 * T12T21) / (2 * sqrt(pi)) *
        bladder_summation * exp_bladder
    ),
    ka_bladder_eq = NULL
  )
}

################################################################################
.krm_swimbladder_low_term <- function(model,
                                      bladder,
                                      bladder_geom,
                                      k_bladder_low,
                                      low_density,
                                      low_sound_speed) {
  ka_bladder_eq <- k_bladder_low * bladder_geom$a_eq_total

  g13 <- bladder$density / low_density
  h13 <- bladder$sound_speed / low_sound_speed
  b0 <- .krm_low_mode_b0(
    k_medium = k_bladder_low,
    k_bladder = model$parameters$acoustics$k_sb,
    a_eq = bladder_geom$a_eq_total,
    g13 = g13,
    h13 = h13
  )

  dx_matrix <- matrix(
    data = rep(bladder_geom$delta_x, each = length(model$parameters$acoustics$k_sw)),
    ncol = length(bladder_geom$delta_x),
    nrow = length(model$parameters$acoustics$k_sw)
  )
  phase_low <- rowSums(
    exp(2i * outer(k_bladder_low, bladder_geom$v_mid)) * dx_matrix,
    na.rm = TRUE
  ) / bladder_geom$length_eq
  Delta_low <- k_bladder_low * bladder_geom$length_eq * cos(bladder$theta)
  sinc_low <- rep(1, length(Delta_low))
  nz_delta <- abs(Delta_low) > sqrt(.Machine$double.eps)
  sinc_low[nz_delta] <- sin(Delta_low[nz_delta]) / Delta_low[nz_delta]

  list(
    f = -1i / pi * bladder_geom$length_eq * sinc_low * b0 * phase_low,
    low_ka_mask = is.finite(ka_bladder_eq) & ka_bladder_eq <= 0.15
  )
}

################################################################################
.krm_swimbladder_term <- function(model, body, bladder, krm_variant) {
  media <- .krm_variant_media(model, krm_variant)
  R23 <- reflection_coefficient(body, bladder)
  bladder_geom <- .krm_bladder_geometry(bladder$rpos, bladder$theta)

  high_term <- .krm_swimbladder_high_term(
    model = model,
    bladder = bladder,
    bladder_geom = bladder_geom,
    R23 = R23,
    T12T21 = 1 - reflection_coefficient(model$medium, model$body)^2,
    k_bladder_high = media$k_bladder_high
  )
  low_term <- .krm_swimbladder_low_term(
    model = model,
    bladder = bladder,
    bladder_geom = bladder_geom,
    k_bladder_low = media$k_bladder_low,
    low_density = media$low_density,
    low_sound_speed = media$low_sound_speed
  )

  f_bladder <- high_term$f
  f_bladder[low_term$low_ka_mask] <- low_term$f[low_term$low_ka_mask]
  f_bladder
}

################################################################################
.krm_combine_terms <- function(f_body, f_bladder = NULL) {
  if (is.null(f_bladder)) {
    sigma_bs <- abs(f_body) * abs(f_body)
    return(list(
      f_bs = f_body,
      sigma_bs = sigma_bs,
      TS = 20 * log10(abs(f_body))
    ))
  }

  f_bs <- f_body + f_bladder
  sigma_body <- abs(f_body) * abs(f_body)
  sigma_bladder <- abs(f_bladder) * abs(f_bladder)
  sigma_bs <- abs(f_bs) * abs(f_bs)

  list(
    f_bs = f_bs,
    sigma_body = sigma_body,
    sigma_bladder = sigma_bladder,
    sigma_bs = sigma_bs,
    TS = 10 * log10(sigma_bs)
  )
}

#' Calculates the theoretical TS using Kirchoff-ray Mode approximation.
#' @noRd
KRM <- function(object) {
  # Detect object class ========================================================
  scatterer_type <- class(object)
  if (inherits(object, "FLS")) {
    object <- .as_krm_profile(object)
  }
  # Extract model parameter inputs =============================================
  model <- extract(object, "model_parameters")$KRM
  # Extract body parameters ====================================================
  body <- extract(object, "body")
  body_terms <- .krm_body_term(model, body, scatterer_type)
  if (scatterer_type == "FLS") {
    totals <- .krm_combine_terms(body_terms$f_body)
    # Define KRM slot for FLS-type scatterer ===================================
    methods::slot(object, "model")$KRM <- data.frame(
      frequency = model$parameters$acoustics$frequency,
      ka = model$parameters$acoustics$k_sw *
        stats::median(body_terms$a_body, na.rm = TRUE),
      f_bs = totals$f_bs,
      sigma_bs = totals$sigma_bs,
      TS = totals$TS
    )
  } else if (scatterer_type == "SBF") {
    bladder <- acousticTS::extract(object, "bladder")
    f_bladder <- .krm_swimbladder_term(
      model = model,
      body = body,
      bladder = bladder,
      krm_variant = model$variant
    )
    totals <- .krm_combine_terms(body_terms$f_body, f_bladder)
    # Define KRM slot for FLS-type scatterer ===================================
    methods::slot(object, "model")$KRM <- data.frame(
      frequency = model$parameters$acoustics$frequency,
      ka = model$parameters$acoustics$k_sw *
        stats::median(body_terms$a_body, na.rm = TRUE),
      f_body = body_terms$f_body,
      f_bladder = f_bladder,
      sigma_body = totals$sigma_body,
      sigma_bladder = totals$sigma_bladder,
      f_bs = totals$f_bs,
      sigma_bs = totals$sigma_bs,
      TS = totals$TS
    )
  }
  # Return object ==============================================================
  object
}
