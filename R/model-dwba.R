################################################################################
# Distorted wave Born approximation (DWBA)
################################################################################
#' Distorted wave Born approximation (DWBA) for weak scatterers
#'
#' @description
#' Calculates backscatter for fluid-like weak scatterers using the distorted
#' wave Born approximation (DWBA). Frequencies in Hz; sound speed in m/s;
#' density in kg/m^3. Material properties may be provided as contrasts or
#' absolute values (contrasts derived relative to seawater).
#'
#' @section Usage:
#' This model is accessed via:
#' \preformatted{
#' target_strength(
#'   ...,
#'   model="dwba",
#'   sound_speed_sw,
#'   density_sw
#' )
#' }
#'
#' @section Arguments:
#' \describe{
#'    \item{\code{sound_speed_sw}}{Seawater sound speed (\eqn{m~s^{-1}}).}
#'    \item{\code{density_sw}}{Seawater density (\eqn{kg~m^{-3}}).}
#' }
#'
#' @section Theory:
#' The DWBA approach is derived under the weak scattering assumption, meaning
#' that the differences in compressibility (\eqn{\kappa}) and density
#' (\eqn{\rho}) between the scatterer and surrounding fluid are sufficiently
#' small to linearize the acoustic scattering problem. This linearization allows
#' the scattered field to be expressed as an integral over the scatterer volume.
#'
#' The far-field backscattering amplitude is given by:
#'
#'  \deqn{
#'    f_{bs} = \frac{k_1}{4\pi} \int \int \int\limits_v
#'    (\gamma_\kappa - \gamma_\rho)
#'    e^{2i k_2 \cdot r_v} dv
#'  }
#'
#' where
#'
#'  \deqn{
#'    \gamma_\kappa = \frac{\kappa_2 - \kappa_1}{\kappa_1}, \quad
#'    \gamma_\rho = \frac{\rho_2 - \rho_1}{\rho_2},
#'  }
#'
#' and
#'
#'  \deqn{
#'    \kappa = (\rho c^2)^{-1}.
#'  }
#'
#' where \eqn{c} is the sound speed (m s<sup>-1</sup>). For elongated bodies,
#' the volume integral reduces to a line intergral along the body axis:
#'
#'  \deqn{
#'    f_{bs} = \frac{k_1}{4} \int_{r_{pos}} (\gamma_\kappa - \gamma_\rho)
#'    e^{2i k_2 \cdot r_{pos}}
#'    \frac{\text{J}_1(2 k_2 a \cos \beta_{tilt})}{\cos \beta_{tilt}}
#'    |dr_{pos}|,
#'  }
#'
#' where \eqn{a} is the local radius and \eqn{\beta_{tilt}} the local tilt
#' angle. The cylindrical Bessel function of the first kind of order 1,
#' \eqn{\text{J}_1}, accounts for the circular cross-section of each segment.
#' The wavenumber \eqn{k_2} is evaluated inside the scatterer.
#'
#' \strong{Assumptions}
#'
#' The DWBA assumes that the scatterer is weakly inhomogenous where the
#' material property contrasts for the interior (\eqn{c_2}, \eqn{\rho_2}) and
#' surrounding fluid (\eqn{c_1}, \eqn{\rho_1}) where:
#'
#'  \deqn{
#'    g = \frac{\rho_2}{\rho_1} \approx 1, \quad
#'    h = \frac{c_2}{c_1} \approx 1.
#'  }
#'
#' In practice, \eqn{c_2} and \eqn{\rho_2} within 5% of \eqn{c_1} and
#' \eqn{\rho_1}, respectively, can be considered to be sufficent for the
#' weak scattering assumption whereby:
#'
#'  \deqn{
#'    |h - 1| \le 0.05, \quad |g - 1| \le 0.05.
#'  }
#'
#' This model also assumes that the scatterer's body has no sharp edges or
#' irregularities, and its cross-section is symmetric around the longitudinal
#' axis. This allows the scattering integral to be reduced to a line integral
#' along the body axis in the first place, simplifying the computation of phase
#' contributions from different segments. Moreover, this enables the body to be
#' discretized into along-axis segments that can approximate arbitrary body
#' shapes. Since the body is axisymmetric and smooth, this further allows the
#' DWBA to be applied for arbitrary orientation angles without additional
#' correction terms.
#'
#'
#' The DWBA provides a first-order approximation that neglects multiple
#' scattering within the body whereby secondary interactions between different
#' parts of the body are ignored. This is valid for weakly scattering objects
#' where the amplitude of scattered waves is small and is consistent with the
#' Born approximation in wave physics. This model also disregards any
#' elastic or shelled effects, treating scatterers are being purely fluid-like.
#' Consequently, the lack of internal elasticity means there is no support for
#' shear waves or resonances due to solid boundaries. Mathematically, this means
#' that only the contrasts in compressibility and density between the scatterer
#' and the surrounding medium contribute to the scattered field.
#'
#' @section Implementation:
#' The model extracts geometric and acoustic parameters from the input object,
#' constructs rotation matrices and wavenumber matrices, and evaluates the
#' integral numerically for each frequency of interest. The DWBA is
#' computationally efficient and handles elongated, weakly scattering targets
#' such as zooplankton and small fish.
#'
#' @seealso
#' See the boundary conditions documentation for
#' more details on how weak scattering relates to other boundary conditions,
#' \code{\link{target_strength}}, \code{\link{FLS}}
#'
#' @references
#' Morse, P.M., and Ingard, K.U. (1968). Theoretical Acoustics. Princeton
#' University Press.
#'
#' Stanton, T.K., Chu, D., and Wiebe, P.H. (1998). Sound scattering by several
#' zooplankton groups. II. Scattering models. The Journal of the Acoustical
#' Society of America, 103, 236-253.
#'
#' @name DWBA
#' @aliases dwba DWBA dwba_curved DWBA_CURVED
#' @docType data
#' @keywords models acoustics internal
NULL

#' Initialize FLS-class object for TS modeling.
#' @param object FLS-class object.
#' @param frequency Frequency vector (Hz).
#' @param sound_speed_sw Seawater sound speed (m/s).
#' @param density_sw Seawater density (kg/m³).
#' @noRd
dwba_initialize <- function(object,
                            frequency,
                            sound_speed_sw = 1500,
                            density_sw = 1026) {
  object_profiled <- .as_dwba_profile(object)
  # Parse shape ================================================================
  shape <- acousticTS::extract(object_profiled, "shape_parameters")
  # Parse body =================================================================
  body <- .hydrate_contrasts(
    acousticTS::extract(object_profiled, "body"),
    sound_speed_sw,
    density_sw
  )
  # Define model parameters recipe =============================================
  model_params <- list(
    acoustics = .init_acoustics_df(
      frequency,
      k_sw = sound_speed_sw,
      k_f = body$h * sound_speed_sw
    ),
    ncyl_b = shape$n_segments
  )
  # Define body parameters recipe ==============================================
  body_params <- list(
    rpos = body$rpos,
    radius = body$radius,
    h = body$h,
    g = body$g,
    theta = body$theta
  )
  .init_model_slots(
    object = object,
    model_name = "DWBA",
    frequency = frequency,
    model_parameters = list(
      parameters = model_params,
      medium = .init_medium_params(sound_speed_sw, density_sw),
      body = body_params
    )
  )
}

#' Calculates the theoretical TS of a fluid-like scatterer at a given frequency
#' using the distorted Born wave approximation (DWBA) model.
#'
#' @param object FLS-class scatterer.
#' @noRd
DWBA <- function(object) {
  # Extract model parameters/inputs ============================================
  model <- acousticTS::extract(object, "model_parameters")$DWBA
  model_body <- acousticTS::extract(model, "body")
  theta <- model_body$theta
  r0 <- model_body$rpos[1:3, , drop = FALSE]
  # Material properties calculation ============================================
  g <- mean(model_body$g)
  h <- mean(model_body$h)
  R <- 1 / (g * h * h) + 1 / g - 2
  # Update position matrices  ==================================================
  rpos <- rbind(r0, a = model$body$radius)
  # Calculate linear scatter response ==========================================
  f_bs <- dwba_fbs_cpp(
    rpos = rpos,
    k_sw = model$parameters$acoustics$k_sw,
    theta = theta,
    h = h,
    R = R
  )
  # Update scatterer object ====================================================
  methods::slot(object, "model")$DWBA <- data.frame(
    frequency = model$parameters$acoustics$frequency,
    ka = model$parameters$acoustics$k_sw * stats::median(model$body$radius),
    f_bs = f_bs,
    sigma_bs = abs(f_bs) * abs(f_bs),
    TS = 10 * log10(abs(f_bs) * abs(f_bs))
  )
  # Return object ==============================================================
  object
}

#' Initialize FLS-class object for TS modeling.
#' @param object FLS-class object.
#' @param frequency Frequency vector (Hz).
#' @param sound_speed_sw Seawater sound speed (m/s).
#' @param density_sw Seawater density (kg/m³).
#' @param radius_curvature_ratio Radius of curvature ratio (length-to-curvature
#' ratio).
#' @param theta Angle of incident soundwave (pi / 2 is broadside).
#' @noRd
dwba_curved_initialize <- function(object,
                                   frequency,
                                   sound_speed_sw = 1500,
                                   density_sw = 1026,
                                   radius_curvature_ratio = NULL,
                                   theta = pi / 2) {
  .Deprecated(
    msg = paste(
      "`DWBA_curved` is deprecated.",
      "Apply `brake()` to the scatterer first, then run `DWBA`",
      "or `target_strength(..., model = \"dwba\")`."
    )
  )
  # Parse shape ================================================================
  shape <- extract(object, "shape_parameters")
  # Parse body =================================================================
  body <- extract(object, "body")
  # Bend body ==================================================================
  body <- brake(body,
                radius_curvature = ifelse(!is.null(radius_curvature_ratio),
                                          radius_curvature_ratio,
                                          body$radius_curvature_ratio
                ),
                mode = "ratio"
  )
  # Define medium parameters ===================================================
  medium_params <- data.frame(
    sound_speed = sound_speed_sw,
    density = density_sw
  )
  # Define model parameters recipe =============================================
  model_params <- list(
    acoustics = data.frame(
      frequency = frequency,
      # Wavenumber (medium) ====================================================
      k_sw = acousticTS::wavenumber(
        frequency,
        sound_speed_sw
      ),
      # Wavenumber (fluid) =====================================================
      k_f = acousticTS::wavenumber(
        frequency,
        body$h * sound_speed_sw
      )
    ),
    ncyl_b = shape$n_segments
  )
  # Define body parameters recipe ==============================================
  body_params <- list(
    rpos = body$rpos,
    radius = body$radius,
    h = body$h,
    g = body$g,
    theta = body$theta,
    radius_curvature_ratio = body$radius_curvature_ratio,
    arc_length = body$arc_length
  )
  # Update object to reflect curvature =========================================

  # Add model parameters slot to scattering object =============================
  methods::slot(
    object,
    "model_parameters"
  )$DWBA_curved <- list(
    parameters = model_params,
    medium = medium_params,
    body = body_params
  )
  # Add model results slot to scattering object ================================
  methods::slot(
    object,
    "model"
  )$DWBA_curved <- data.frame(
    frequency = frequency,
    sigma_bs = rep(
      NA,
      length(frequency)
    )
  )
  # Output =====================================================================
  return(object)
}

#' Calculates the theoretical TS of a fluid-like scatterer at a given frequency
#' using the distorted Born wave approximation (DWBA) model.
#'
#' @param object FLS-class scatterer.
#' @noRd
DWBA_curved <- function(object) {
  # Extract model parameters/inputs ============================================
  model <- acousticTS::extract(object, "model_parameters")$DWBA_curved
  body <- acousticTS::extract(object, "body")
  # Parse position matrices ====================================================
  r0 <- model$body$rpos[1:3, ]
  theta <- model$body$theta
  # Material properties calculation ============================================
  g <- body$g
  h <- body$h
  R <- 1 / (g * h * h) + 1 / g - 2
  # Update position matrices  ==================================================
  rpos <- rbind(r0, a = model$body$radius)
  # Calculate linear scatter response ==========================================
  f_bs <- dwba_fbs_cpp(
    rpos = rpos,
    k_sw = model$parameters$acoustics$k_sw,
    theta = theta,
    h = h,
    R = R
  )
  # Update scatterer object ====================================================
  methods::slot(object, "model")$DWBA_curved <- data.frame(
    frequency = model$parameters$acoustics$frequency,
    ka = model$parameters$acoustics$k_sw * stats::median(model$body$radius),
    f_bs = f_bs,
    sigma_bs = abs(f_bs) * abs(f_bs),
    TS = 10 * log10(abs(f_bs) * abs(f_bs))
  )
  # Return object ==============================================================
  object
}
