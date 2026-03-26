################################################################################
# Sphere modal series solution
################################################################################
#' Spherical modal series (SPHMS) solution
#'
#' @description
#' Spherical modal series (SPHMS) backscatter model. Supports rigid,
#' pressure-release, fluid-filled, shelled (pressure-release/liquid/gas).
#' Frequencies in Hz; sound speed in m/s; density in kg/m^3. Material
#' properties may be contrasts or absolute (contrasts derived relative to
#' seawater). Full equations and boundary-condition details are documented in
#' the SPHMS vignette.
#'
#' @section Usage:
#' This model is accessed via:
#' \preformatted{
#' target_strength(
#'   ...,
#'   model="sphms",
#'   boundary,
#'   sound_speed_sw,
#'   density_sw,
#'   m_limit
#' )
#' }
#'
#' @section Arguments:
#' \describe{
#'   \item{\code{boundary}}{Boundary condition at a spherical surface.
#'     One of \code{"fixed_rigid"}, \code{"pressure_release"},
#'     \code{"liquid_filled"}, \code{"gas_filled"},
#'     \code{"shelled_pressure_release"}, \code{"shelled_liquid"}, or
#'     \code{"shelled_gas"}. See the
#'     \link[=boundary_conditions]{boundary conditions documentation} for more
#'     details on these different boundary conditions.}
#'   \item{\code{sound_speed_sw}}{Seawater sound speed (\eqn{m~s^{-1}}).}
#'   \item{\code{density_sw}}{Seawater density (\eqn{kg~m^{-3}}).}
#'   \item{\code{m_limit}}{Optional model truncation limit used to cap the
#'   number of modes in the numerical calculation.}
#' }
#'
#' @section Theory:
#' The modal series solution for a sphere expands backscatter from a spherical
#' scatterer in spherical Bessel/Hankel modes. The form function is:
#'
#' \eqn{
#'  f_{bs} = -\frac{i}{k_w} \sum_{n=0}^{\infty} (-1)^n (2n+1) A_n,
#' }
#'
#' @seealso
#' \code{\link{target_strength}}, \code{\link{FLS}}, \code{\link{GAS}},
#' \code{\link{ESS}}, \code{\link{Sphere}}, \code{\link{sphere}}
#'
#' @references
#' Anderson, V.C. (1950). Sound scattering from a fluid sphere. The Journal of
#' The Acoustical Society of America, 22: 426–431.
#'
#' @name SPHMS
#' @aliases sphms SPHMS
#' @docType data
#' @keywords models acoustics internal
NULL

#' Initialize object for the modal series solution for a sphere
#'
#' @param object Scatterer-class object that must be a sphere.
#' @param frequency Acoustic frequency (Hz).
#' @param boundary The sphere's boundary-type.
#' @param sound_speed_sw Seawater sound speed (\eqn{m s^{-1}}).
#' @param density_sw Seawater density (\eqn{kg m^{-3}}).
#' @param m_limit Optional argument for truncating the number of modes
#' used by the modal series solution.
#' @keywords internal
#' @noRd
sphms_initialize <- function(object,
                             frequency,
                             boundary = NULL,
                             sound_speed_sw = .SEAWATER_SOUND_SPEED_DEFAULT,
                             density_sw = .SEAWATER_DENSITY_DEFAULT,
                             m_limit = NULL) {

  # Detect object class ========================================================
  scatterer_type <- class(object)
  # Detect object shape ========================================================
  scatterer_shape <- acousticTS::extract(object, "shape_parameters")
  # Validate shape +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  if (scatterer_shape$shape != "Sphere") {
    stop(
      "The modal series solution for a sphere requires scatterer to be ",
      "shape-type 'Sphere'. Input scatterer is shape-type ",
      paste0("'", scatterer_shape, "'.")
    )
  }
  medium_params <- .init_medium_params(sound_speed_sw, density_sw)
  # Determine expansion coefficient Bm method ==================================
  # Validate method ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  if (is.null(boundary)) {
    # Set default boundary +++++++++++++++++++++++++++++++++++++++++++++++++++
    boundary <- if (methods::is(object, "CAL")) {
      stop(
        "Use 'model='calibration'' when modeling the TS of a solid sphere."
      )
    } else if (methods::is(object, "ESS")) {
      "shelled_liquid"
    } else if (methods::is(object, "FLS")) {
      "liquid_filled"
    } else if (methods::is(object, "GAS") || methods::is(object, "SBF")) {
      "gas_filled"
    } else {
      stop(
        "Could not infer a default SPHMS boundary for scatterer class '",
        class(object)[1], "'.",
        call. = FALSE
      )
    }
  }
  if (!(boundary %in% c(
    "liquid_filled", "fixed_rigid", "pressure_release", "gas_filled",
    "shelled_pressure_release", "shelled_liquid", "shelled_gas"
  ))) {
    # Check boundary compatibility with scattering class +++++++++++++++++++++++
    # ESS ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if (
      methods::is(object, "ESS") &
      !(boundary %in%
        c("shelled_pressure_release", "shelled_liquid", "shelled_gas")
      )
    ) {
      stop(
        "Only the following values for 'boundary' are available in this ",
        "implementation of the sphere modal series solution for the ",
        "'ESS'-class: ",
        "'shelled_pressure_release', 'shelled_liquid', 'shelled_gas'. Input ",
        "boundary is ", paste0("'", boundary, "'.")
      )
    } else if (
      # All others +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      boundary %in%
      c(
        "shelled_pressure_release", "shelled_liquid", "shelled_gas"
      )
    ) {
      stop(
        "Only the following values for 'boundary' are available in this ",
        "implementation of the sphere modal series solution for the ",
        paste0("'", class(object)[1], "'-class:" ),
        "'fixed_rigid', 'pressure_release', 'liquid_filled', 'gas_filled'."
      )
    }
  }
  # Assign method ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  model_params <- list()
  model_params$Bm_method <- switch(
    boundary,
    liquid_filled = "Bm_fluid",
    gas_filled = "Bm_fluid",
    fixed_rigid = "Bm_rigid",
    pressure_release = "Bm_pressure_release",
    shelled_pressure_release = "Bm_shelled_pressure_release",
    shelled_liquid = "Bm_shelled_fluid",
    shelled_gas = "Bm_shelled_fluid"
  )
  # Compute sphere material properties =========================================
  exterior <- if ("shell" %in% methods::slotNames(object)) {
    acousticTS::extract(object, "shell")
  } else {
    acousticTS::extract(object, "body")
  }
  g_exterior <- exterior$g %||%
    if (!is.null(exterior$density)) exterior$density / density_sw else NA_real_
  h_exterior <- exterior$h %||%
    if (!is.null(exterior$sound_speed)) exterior$sound_speed / sound_speed_sw else NA_real_
  # Material properties: body/shell-to-seawater interface ++++++++++++++++++++++
  if (methods::is(object, "ESS")) {
    g21 <- if (!is.null(exterior$density)) exterior$density / density_sw else
      exterior$g %||% NA_real_
    h21 <- if (!is.null(exterior$sound_speed)) exterior$sound_speed / sound_speed_sw else
      exterior$h %||% NA_real_
  } else {
    g21 <- g_exterior
    h21 <- h_exterior
  }
  # Material properties: fluid-to-shell interface ++++++++++++++++++++++++++++++
  if (methods::is(object, "ESS")) {
    fluid <- acousticTS::extract(object, "fluid")
    g32 <- if (!is.null(fluid$density)) fluid$density / (g21 * density_sw) else
      if (!is.null(fluid$g)) fluid$g / g21 else NA_real_
    h32 <- if (!is.null(fluid$sound_speed)) fluid$sound_speed / (h21 * sound_speed_sw) else
      if (!is.null(fluid$h)) fluid$h / h21 else NA_real_
  } else {
    g32 <- NA_real_
    h32 <- NA_real_
  }
  # Material properties: fluid-to-seawater interface +++++++++++++++++++++++++++
  if (methods::is(object, "ESS")) {
    g31 <- if (!is.null(fluid$density)) fluid$density / density_sw else
      fluid$g %||% NA_real_
    h31 <- if (!is.null(fluid$sound_speed)) fluid$sound_speed / sound_speed_sw else
      fluid$h %||% NA_real_
  } else {
    g31 <- g21
    h31 <- h21
  }
  # Outer radius +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  radius_shell <- switch(
    class(object),
    ESS = exterior$radius,
    NA
  )
  radius_fluid <- switch(
    class(object),
    ESS = fluid$radius,
    exterior$radius
  )
  # Add to list ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  body_params <- list(
    g21 = g21, g32 = g32, g31 = g31,
    h21 = h21, h32 = h32, h31 = h31,
    radius_shell = radius_shell,
    radius_fluid = radius_fluid
  )
  # Define model parameters recipe =============================================
  acoustics <- .init_acoustics_df(frequency, k_sw = sound_speed_sw)
  # Wavenumber at shell-to-seawater interface ++++++++++++++++++++++++++++++++++
  acoustics$k_shell <- acoustics$k_sw / h21
  # Wavenumber at fluid-to-shell interface +++++++++++++++++++++++++++++++++++++
  acoustics$k_fluid <- if( !is.na(h32)) {
    acoustics$k_shell / h32
  } else {
    acoustics$k_sw / h31
  }
  model_params$acoustics <- acoustics
  # Define limits for 'm' modal series iterator ================================
  if (!is.null(m_limit)) {
    model_params$acoustics$m_limit <- m_limit
  } else if(!is.na(radius_shell)){
    model_params$acoustics$m_limit <- round(
      model_params$acoustics$k_sw * radius_shell + 20
    )
  } else {
    model_params$acoustics$m_limit <- round(
      model_params$acoustics$k_sw * radius_fluid + 20
    )
  }
  .init_model_slots(
    object = object,
    model_name = "SPHMS",
    frequency = frequency,
    model_parameters = list(
      parameters = model_params,
      medium = medium_params,
      body = body_params
    )
  )
}


#' Calculate the TS of a sphere using the modal series solution
#' @noRd
SPHMS <- function(object) {
  # Extract model parameters/inputs ============================================
  model_params <- extract(object, "model_parameters")$SPHMS
  parameters <- model_params$parameters
  acoustics <- parameters$acoustics
  medium <- model_params$medium
  body <- model_params$body
  # Pull out maximum 'm' specifically ==========================================
  m_limit <- acoustics$m_limit
  # Compute Delta (shell thickness) ============================================
  b_delta <- body$radius_shell - body$radius_fluid
  # Compute ka and kb ==========================================================
  k1a <- if (!is.na(body$radius_shell)) {
    acoustics$k_sw * body$radius_shell
  } else {
    acoustics$k_sw * body$radius_fluid
  }
  k2a <- acoustics$k_shell * body$radius_shell
  k2b <- acoustics$k_shell * body$radius_fluid
  k3a <- acoustics$k_fluid * body$radius_fluid
  k3b <- acoustics$k_fluid * body$radius_fluid
  # Define the boundary expansion coefficient method ===========================
  sph_bm <- switch(
    parameters$Bm_method,
    Bm_rigid = .sphms_bm_rigid(k1a, m_limit),
    Bm_pressure_release = .sphms_bm_prelease(k1a, m_limit) ,
    Bm_fluid = .sphms_bm_fluid(k1a, k3a, body$g31, body$h31, m_limit),
    Bm_shelled_pressure_release = .sphms_bm_shelled_prelease(
        k1a, k2a, k2b, body$g21, body$h21, m_limit
      ),
    Bm_shelled_fluid = .sphms_bm_shelled_fluid(
        k1a, k2a, k2b, k3b,
        body$g21, body$g31, body$g32,
        body$h21, body$h31, body$h32,
        m_limit
      )
  )
  # Compute the linear scattering coefficient, f_bs ============================
  f_bs <- -1i / acoustics$k_sw * sph_bm
  # Calculate backscatter and return ===========================================
  # Compute sigma_bs +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  o_bs <- abs(f_bs)
  methods::slot(object, "model")$SPHMS <- data.frame(
    frequency = acoustics$frequency,
    f_bs = f_bs,
    sigma_bs = o_bs,
    TS = 20 * log10(o_bs)
  )
  object
}

#' Helper function for fixed rigid sphere
#' @keywords internal
#' @noRd
.sphms_bm_rigid <- function(k1a, m_limit) {
  m_max <- max(m_limit)
  Am <- .modal_series_apply(
    m_limit = m_limit,
    FUN = function(k1a, ml) {
      -(jsd(0:ml, k1a) / hsd(0:ml, k1a))
    },
    k1a = k1a
  )
  .modal_weighted_sum(Am, (2 * (0:m_max) + 1) * (-1)^(0:m_max))
}

#' Helper function for pressure release sphere
#' @keywords internal
#' @noRd
.sphms_bm_prelease <- function(k1a, m_limit) {
  m_max <- max(m_limit)
  Am <- .modal_series_apply(
    m_limit = m_limit,
    FUN = function(k1a, ml) {
      -(js(0:ml, k1a) / hs(0:ml, k1a))
    },
    k1a = k1a
  )
  .modal_weighted_sum(Am, (2 * (0:m_max) + 1) * (-1)^(0:m_max))
}

#' Helper function for fluid sphere
#' @keywords internal
#' @noRd
.sphms_bm_fluid <- function(k1a, k3a, g31, h31, m_limit) {
  # Get material properties product ============================================
  gh <- g31 * h31
  m_max <- max(m_limit)
  cm_num <- .modal_series_apply(
    m_limit = m_limit,
    FUN = function(k1a, k3a, ml) {
      m <- 0:ml
      (jsd(m, k3a) * ys(m, k1a)) / (js(m, k3a) * jsd(m, k1a)) -
        (gh * ysd(m, k1a) / jsd(m, k1a))
    },
    k1a = k1a,
    k3a = k3a
  )
  cm_denom <- .modal_series_apply(
    m_limit = m_limit,
    FUN = function(k1a, k3a, ml) {
      m <- 0:ml
      (jsd(m, k3a) * js(m, k1a)) / (js(m, k3a) * jsd(m, k1a)) - gh
    },
    k1a = k1a,
    k3a = k3a
  )
  cm <- cm_num / cm_denom
  Am <- (-1 / (1 + 1i * cm))
  .modal_weighted_sum(Am, (2 * (0:m_max) + 1) * (-1)^(0:m_max))
}

#' Helper function for shelled pressure release sphere
#' @keywords internal
#' @noRd
.sphms_bm_shelled_prelease <- function(k1a, k2a, k2b, g21, h21, m_limit) {
  # Get maxmimum m_limit =======================================================
  m_max <- max(m_limit)
  # Calculate the expansion coefficient, A_m ===================================
  Am <- mapply(FUN = function(k1a, k2a, k2b, g21, h21, ml) {
    # Expand modal series iterator =============================================
    m <- 0:ml
    # Compute elements for boundary matrix +++++++++++++++++++++++++++++++++++++
    a11 <- -hs(m, k1a)
    a12 <- js(m, k1a)
    a21 <- -g21 * h21 * hsd(m, k1a)
    b1 <- js(m, k1a)
    b2 <- jsd(m, k1a) * g21 * h21
    d1 <- js(m, k2a) * ys(m, k2b) - js(m, k2b) * ys(m, k2a)
    d2 <- jsd(m, k2a) * ys(m, k2b) - js(m, k2b) * ysd(m, k2a)
    # Calculate coefficient ++++++++++++++++++++++++++++++++++++++++++++++++++++
    Am <- (b1 * d2 - d1 * b2) / (a11 * d2 - d1 * a21)
    if (length(Am) < (m_max + 1)) {
      difference <- (m_max + 1) - length(Am)
      Am <- c(Am, rep(NA, difference))
    }
    Am
  }, k1a, k2a, k2b, g21, h21, m_limit)
  # Return the entire boundary modal term ======================================
  colSums((2 * (0:m_max) + 1) * (-1)^(0:m_max) * Am, na.rm = TRUE)
}

#' Helper function for shelled fluid sphere
#' @keywords internal
#' @noRd
.sphms_bm_shelled_fluid <- function(
    k1a, k2a, k2b, k3b, g21, g31, g32, h21, h31, h32, m_limit
) {
  # Get maximum m_limit ========================================================
  m_max <- max(m_limit)
  # Compute boundary conditions ================================================
  Am <- mapply(FUN = function(k1a, k2a, k2b, k3b, g21, g31, g32, h21, h31, h32,
                              ml) {
    # Expand modal series iterator =============================================
    m <- 0:ml
    # Get material properties product ==========================================
    gh <- g31 * h31
    # Calculate the expansion coefficient, A_m =================================
    # Compute elements for boundary matrix +++++++++++++++++++++++++++++++++++++
    a11 <- -hs(m, k1a)
    a12 <- js(m, k2a)
    a13 <- ys(m, k2a)
    a21 <- -g21 * h21 * hsd(m, k1a)
    a22 <- jsd(m, k2a)
    a23 <- ysd(m, k2a)
    a32 <- js(m, k2b) * jsd(m, k3b) - g32 * h32 * jsd(m, k2b) *
      js(m, k3b)
    a33 <- ys(m, k2b) * jsd(m, k3b) - g32 * h32 * ysd(m, k2b) *
      js(m, k3b)
    b1 <- js(m, k1a)
    b2 <- jsd(m, k1a) * g21 * h21
    # Calculate coefficient ++++++++++++++++++++++++++++++++++++++++++++++++++++
    Am <- (b1*a22*a33 + a13*b2*a32 - a12*b2*a33 - b1*a23*a32) /
      (a11*a22*a33 + a13*a21*a32 - a12*a21*a33 - a11*a23*a32)
    if (length(Am) < (m_max + 1)) {
        difference <- (m_max + 1) - length(Am)
        Am <- c(Am, rep(NA, difference))
      }
    Am
    }, k1a, k2a, k2b, k3b, g21, g31, g32, h21, h31, h32, m_limit)

  # Return the entire boundary modal term ======================================
  colSums((2 * (0:m_max) + 1) * (-1)^(0:m_max) * Am, na.rm = TRUE)
}
