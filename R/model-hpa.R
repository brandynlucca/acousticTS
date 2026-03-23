################################################################################
# Ray-based high pass approximation
################################################################################
#' High-pass approximation (HPA) scattering model
#'
#' @description
#' Computes the far-field scattering amplitude and related quantities spherical
#' and elongated scatterers using a high-pass approximation model. The
#' high-pass model provides a simplified analytical expression for the
#' backscattering cross-section that is valid for all values of \eqn{ka},
#' where \eqn{k} is the acoustic wavenumber and \eqn{a} is a characteristic
#' dimension of the scatterer. The model is named after the analogy to a
#' two-pole high-pass filter in electrical circuit theory, where the frequency
#' response resembles the backscattering behavior as a function of \eqn{ka}.
#' The high-pass model is computationally efficient and captures the essential
#' physics of acoustic scattering from weakly scattering bodies,
#' including spheres, prolate spheroids, and straight or bent cylinders.
#' Two implementations are available: the Johnson (1977) formulation for fluid
#' spheres, and the Stanton (1989) generalization for spheres,
#' prolate spheroids, and cylinders.
#'
#' @section Usage:
#' This model is accessed via:
#' \preformatted{
#' target_strength(
#'   ...,
#'   model = "HPA",
#'   method,
#'   deviation_fun,
#'   null_fun,
#'   sound_speed_sw,
#'   density_sw
#' )
#' }
#'
#' @section Arguments:
#' \describe{
#'   \item{\code{method}}{Method for computing the high-pass model. One of
#'   \code{"johnson"} (Johnson 1977 formulation for fluid spheres) or
#'   \code{"stanton"} (Stanton 1989 formulation for spheres, prolate spheroids,
#'   and cylinders).}
#'   \item{\code{deviation_fun}}{Function or scalar specifying the expected
#'   deviation in the phase of the scattered wave due to shape complexity,
#'   flexure, or other factors. Expressed as a function of \eqn{ka}, or as a
#'   scalar constant. Default is 1.}
#'   \item{\code{null_fun}}{Function or scalar specifying the expected
#'   reduction in scattering amplitude due to nulls in the scattering pattern.
#'   Expressed as a function of \eqn{ka}, or as a scalar constant. Default is
#'   1.}
#'   \item{\code{sound_speed_sw}}{Seawater sound speed (\eqn{m~s^{-1}}).}
#'   \item{\code{density_sw}}{Seawater density (\eqn{kg~m^{-3}}).}
#' }
#'
#' @section Theory:
#' The high-pass model is based on the observation that the backscattering
#' cross-section of a scatterer can be approximated by a simple analytical
#' expression that captures the low-frequency Rayleigh scattering regime
#' (where \eqn{\sigma_{bs} \propto (ka)^4}) and the high-frequency geometric
#' scattering regime (where \eqn{\sigma_{bs}} approaches a constant). The
#' transition between these regimes is governed by an auxiliary material
#' property parameter, \eqn{\alpha}, which depends on the density contrast
#' \eqn{g} and sound speed contrast \eqn{h} between the scatterer and the
#' surrounding medium. For spheres and prolate spheroids, the auxiliary
#' parameter is given by
#'
#'  \deqn{
#'    \alpha_{\pi s} = \frac{1 - g h^2}{3 g h^2} + \frac{1 - g}{1 + 2g}
#'  }
#'
#' and for cylinders and prolate spheroids (in the cylindrical limit), the
#' auxiliary parameter is
#'
#'  \deqn{
#'    \alpha_{\pi c} = \frac{1 - g h^2}{2 g h^2} + \frac{1 - g}{1 + g}
#'  }
#'
#' where \eqn{g} is the density contrast (target to medium) and \eqn{h} is the
#' sound speed contrast (target to medium). The reflection coefficient at the
#' scatterer surface is given by
#'
#'  \deqn{
#'    \mathcal{R} = \frac{gh - 1}{gh + 1}
#'  }
#'
#' The Johnson (1977) high-pass model for a fluid sphere is expressed as
#'
#'  \deqn{
#'    \sigma_{bs} = \frac{a^2 (ka)^4 \alpha_{\pi s}^2}{1 + \frac{3}{2} (ka)^4}
#'  }
#'
#' where \eqn{a} is the spherical radius. This formulation is valid for all
#' \eqn{ka} and provides a good approximation for fluid spheres.
#' The Stanton (1989) high-pass model generalizes this approach to spheres,
#' prolate spheroids, and cylinders, and incorporates empirical terms to
#' account for nulls in the scattering pattern and deviations due to shape
#' complexity. For a sphere, the Stanton formulation is
#'
#'  \deqn{
#'    \sigma_{bs} = \frac{
#'      a^2 (ka)^4 \alpha_{\pi s}^2 \mathcal{G}
#'    }{
#'      1 + \frac{4(ka)^4 \alpha_{\pi s}^2}{\mathcal{R}^2 \mathcal{F}}
#'    }
#'  }
#'
#' where \eqn{\mathcal{G}} is a null function that accounts for reductions in
#' scattering amplitude at certain frequencies, and \eqn{\mathcal{F}} is a
#' deviation function that accounts for phase variability. For a prolate
#' spheroid, the Stanton formulation is
#'
#'  \deqn{
#'    \sigma_{bs} = \frac{
#'      \frac{1}{9} L^2 (ka)^4 \alpha_{\pi c}^2 \mathcal{G}
#'    }{
#'      1 + \frac{\frac{16}{9}(ka)^4 \alpha_{\pi c}^2}
#'      {\mathcal{R}^2 \mathcal{F}}
#'    }
#'  }
#'
#' where \eqn{L} is the length of the prolate spheroid. For a straight cylinder
#' at angle \eqn{\theta}, the Stanton formulation is
#'
#'  \deqn{
#'    \sigma_{bs} = \frac{
#'      \frac{1}{4} L^2 (Ka)^4 \alpha_{\pi c}^2 s^2 \mathcal{G}
#'    }{
#'      1 + \frac{\pi (Ka)^3 \alpha_{\pi c}^2}{\mathcal{R}^2 \mathcal{F}}
#'    }
#'  }
#'
#' where \eqn{K = k \sin \theta}, \eqn{a} is the cylindrical radius, and
#' \eqn{s = \sin(kL \cos \theta) / (kL \cos \theta)} accounts for the finite
#' length and orientation of the cylinder. For a bent cylinder with radius of
#' curvature \eqn{\rho_c} relative to the cylinder length \eqn{L}, the Stanton
#' formulation is
#'
#'  \deqn{
#'  \sigma_{bs} = \frac{
#'      \frac{1}{4} L^2 (ka)^4 \alpha_{\pi c}^2 \mathcal{H}^2 \mathcal{G}
#'    }{
#'      1 + \frac{
#'        L^2 (ka)^4 \alpha_{\pi c}^2 \mathcal{H}^2
#'      }{
#'        \rho_c a \mathcal{R}^2 \mathcal{F}
#'      }
#'    }
#'  }
#'
#' where \eqn{\mathcal{H} = \frac{1}{2} + \frac{1}{2} (\rho_c / L)
#' \sin(L / \rho_c)} is an effective length factor that accounts for the
#' curvature. The empirical functions \eqn{\mathcal{F}} and \eqn{\mathcal{G}}
#' are specified by the user as functions of \eqn{ka}, and can be used to fit
#' the model to measured or numerically simulated scattering data. In all
#' cases, the wavenumber satisfies \eqn{kr \gg 1} and
#' \eqn{L \ll 2\sqrt{r \lambda}}, where \eqn{r} is the  distance from the
#' scatterer to the receiver and \eqn{\lambda} is the acoustic wavelength,
#' ensuring that the far-field and elongated object assumptions are valid.
#'
#' @references
#' Johnson, R.K. (1977). Sound scattering from a fluid sphere revisited. The
#' Journal of the Acoustical Society of America, 61: 375-377.
#'
#' Johnson, R.K. (1978). Erratum: Sound scattering from a fluid sphere
#' revisited. The Journal of the Acoustical Society of America, 63: 626.
#'
#' Stanton, T.K. (1989). Simple approximate formulas for backscattering of
#' sound by spherical and elongated objects. The Journal of the Acoustical
#' Society of America, 86: 1499-1510.
#'
#' @name HP
#' @aliases hpa HPA high_pass HIGH_PASS
#' @docType data
#' @keywords models acoustics
NULL

#' Initialize Scatterer-class object for the high-pass approximation model
#'
#' @param object Scatterer-class object.
#' @param frequency Frequency vector (Hz).
#' @param method Method for computing the high-pass model.
#' @param deviation_fun Deviation function.
#' @param null_fun Null positioning function.
#' @param sound_speed_sw Seawater sound speed (m/s).
#' @param density_sw Seawater density (kg/m³).
#' @noRd
hpa_initialize <- function(object,
                           frequency,
                           method = "stanton",
                           deviation_fun = 1,
                           null_fun = 1,
                           sound_speed_sw = 1500,
                           density_sw = 1026) {
  # Validate function arguments ================================================
  if (! method %in% c("stanton", "johnson")) {
    stop("'method' must be one of the following: 'stanton', 'johnson'.")
  }
  if (is.function(deviation_fun)) {
    if (length(formals(deviation_fun)) > 1) {
      stop(
        "'deviation_fun' as a callable function must only comprise a single ",
        "argument for ka."
      )
    }
  } else if (length(deviation_fun) > 1) {
    stop(
      "'deviation_fun' as a numeric must be scalar."
    )
  }
  if (is.function(null_fun)) {
    if (length(formals(null_fun)) > 1) {
      stop(
        "'null_fun' as a callable function must only comprise a single ",
        "argument for ka."
      )
    }
  } else if (length(null_fun) > 1) {
    stop(
      "'null_fun' as a numeric must be scalar."
    )
  }
  # Extract scatterer shape information ========================================
  shape_params <- acousticTS::extract(object, "shape_parameters")
  if (methods::is(object, "ESS")) {
    scatterer_shape <- shape_params$shape
    body <- acousticTS::extract(object, "shell")
  } else if ("fluid" %in% names(shape_params)) {
    scatterer_shape <- shape_params$body$shape
    body <- acousticTS::extract(object, "fluid")
  } else {
    scatterer_shape <- shape_params$shape
    body <- acousticTS::extract(object, "body")
  }
  if (method == "stanton") {
    if (! scatterer_shape %in% c("Sphere", "ProlateSpheroid", "Cylinder")) {
      stop(
        "The high-pass approximation requires scatterer to one of the ",
        "following shape-types when 'method' = 'stanton': 'Cylinder', ",
        "'ProlateSpheroid', 'Sphere'. Input scatterer is shape-type ",
        paste0("'", scatterer_shape, "'.")
      )
    }
  } else if (scatterer_shape != "Sphere") {
    stop(
      "The high-pass approximation requires scatterer to one of shape-type ",
      "'Sphere'. Input scatterer is shape-type ",
      paste0("'", scatterer_shape, "'.")
    )
  }
  medium_params <- .init_medium_params(sound_speed_sw, density_sw)
  body <- .hydrate_contrasts(body, medium_params$sound_speed, medium_params$density)
  # Define model parameters recipe =============================================
  model_params <- list(
    acoustics = .init_acoustics_df(
      frequency,
      k_sw = sound_speed_sw,
      k_f = body$h * sound_speed_sw
    ),
    deviation_fun = deviation_fun,
    null_fun = null_fun,
    method = method
  )
  # Define body parameters recipe ==============================================
  body_params <- c(
    body,
    list(
      theta = body$theta,
      shape = shape_params$shape,
      length = if ("body" %in% names(shape_params)) {
        shape_params$body$length
      } else {
        shape_params$length
      }
    )
  )
  # Fill out material properties ===============================================
  body_params$h <- body$h
  body_params$g <- body$g
  .init_model_slots(
    object = object,
    model_name = "HPA",
    frequency = frequency,
    model_parameters = list(
      parameters = model_params,
      medium = medium_params,
      body = body_params
    )
  )
}

#' Auxiliary variable term, alpha, used in the high-pass models proposed by
#' Johnson (1977, 1978) and Stanton (1989).
#' @noRd
.alpha_pi <- function(body) {
  # Compute the shape specific auxiliary variable ==============================
  if (body$shape == "Sphere") {
    return(
      (1 - body$g * body$h^2) /
        (3 * body$g * body$h^2) + (1 - body$g) / (1 + 2 * body$g)
    )
  } else {
    (1 - body$g * body$h^2) / (2 * body$g * body$h^2) +
      (1 - body$g) / (1 + body$g)
  }
}

#' Johnson (1977) high-pass model for a fluid sphere and the accompanying
#' erratum from Johnson (1978).
#' @noRd
.johnson_hp <- function(acoustics, body, alpha) {
  # Calculate ka ===============================================================
  if (length(body$radius) > 1) {
    a <- max(body$radius)
  } else {
    a <- body$radius
  }
  ka <- acoustics$k_sw * a
  # Calculate sigma_bs =========================================================
  (body$radius^2 * ka^4 * alpha^2) / (1 + 1.5 * ka^4)
}

#' Stanton (1989) high-pass model for spheres.
#' @noRd
.stanton_sphere <- function(a, ka, r, alpha, gnull, fdevs) {
  # Compute sigma_bs ===========================================================
  (a^2 * ka^4 * alpha^2 * gnull) / (1 + (4 * ka^4 * alpha^2) / (r^2 * fdevs))
}

#' Stanton (1989) high-pass model for prolate spheroids.
#' @noRd
.stanton_prolate_spheroid <- function(ka, l, r, alpha, gnull, fdevs) {
  # Compute sigma_bs ===========================================================
  ((1 / 9) * l^2 * ka^4 * alpha^2 * gnull) / (1 + ((16 / 9) * ka^4 * alpha^2) /
                                                (r^2 * fdevs))
}

#' Stanton (1989) high-pass model for cylinders.
#' @noRd
.stanton_cylinder <- function(k, l, a, theta, rho_c, r, alpha, gnull, fdevs) {
  # Compute ka =================================================================
  ka <- k * a
  # Special case: bent cylinder ================================================
  if (length(rho_c) > 0 && !all(is.na(rho_c))) {
      # Compute the literal radius of curvature ++++++++++++++++++++++++++++++++
      r_c <- rho_c * l
      # Compute effective length of bent cylinder ++++++++++++++++++++++++++++++
      H <- 0.5 + 0.5 * rho_c * sin(1 / rho_c)
      # Compute sigma_bs +++++++++++++++++++++++++++++++++++++++++++++++++++++++
      return(
        (0.25 * l^2 * ka^4 * alpha^2 * H^2 * gnull) /
          (1 + (l^2 * ka^4 * alpha^2 * H^2) / (r_c * a * r^2 * fdevs))
      )
  }
  # Compute the effective length based on tilt angle ===========================
  s <- sin(k * l * cos(theta)) / (k * l * cos(theta))
  # Compute the wavenumber*radius based on the tilt angle ======================
  Ka <- ka * sin(theta)
  # Compute sigma_bs ===========================================================
  (0.25 * l^2 * Ka^4 * alpha^2 * s^2 * gnull) /
    (1 + (pi * Ka^3 * alpha^2) / (r^2 * fdevs))
}

#' Stanton (1989) high-pass model for spheres, prolate spheroids, and (bent)
#' cylinders.
#' @noRd
.stanton_hp <- function(acoustics, body, parameters, alpha) {
  # Calculate ka ===============================================================
  if (length(body$radius) > 1) {
    a <- max(body$radius)
  } else {
    a <- body$radius
  }
  ka <- acoustics$k_sw * a
  # Resolve empirical terms, G (nulls) and F (deviation) =======================
  # G ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  if (is.function(parameters$null_fun)) {
    gnull <- parameters$null_fun(ka)
  } else {
    gnull <- parameters$null_fun
  }
  # F ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  if (is.function(parameters$deviation_fun)) {
    fdevs <- parameters$deviation_fun(ka)
  } else {
    fdevs <- parameters$deviation_fun
  }
  # Calculate corrected reflection coefficient, R ==============================
  r <- (body$g * body$h - 1) / (body$g * body$h + 1)
  # Compute sigma_bs ===========================================================
  switch(
    body$shape,
    Sphere = .stanton_sphere(a, ka, r, alpha, gnull, fdevs),
    Cylinder = .stanton_cylinder(
      acoustics$k_sw, body$length, a, body$theta,
      body$radius_curvature_ratio, r, alpha, gnull, fdevs
      ),
    ProlateSpheroid = .stanton_prolate_spheroid(
      ka, body$length, r, alpha, gnull, fdevs
    )
  )
}

#' High-pass approximation (HPA) model
#'
#' Calculates the far-field scattering amplitude and related quantities for a
#' finite cylinder using the high-pass approximation model.
#'
#' @param object Scatterer-object with a Cylinder-class shape.
#' @noRd
HPA <- function(object) {
  # Extract model parameters/inputs ============================================
  model <- acousticTS::extract(object, "model_parameters")$HPA
  parameters <- model$parameters
  acoustics <- parameters$acoustics
  body <- model$body
  # Compute the alpha auxiliary variable =======================================
  alpha <- .alpha_pi(body)
  # Calculate the linear backscattering coefficient, sigma_bs ==================
  obs <- switch(
    parameters$method,
    johnson = .johnson_hp(acoustics, body, alpha),
    stanton = .stanton_hp(acoustics, body, parameters, alpha)
  )
  # Update and return modeled TS ===============================================
  methods::slot(object, "model")$HPA <- data.frame(
    frequency = acoustics$frequency,
    sigma_bs = obs,
    TS = 10 * log10(obs)
  )
  object
}
