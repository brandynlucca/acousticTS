################################################################################
# ACOUSTIC AND SIGNAL PROCESSING VARIABLE CALCULATIONS
################################################################################
################################################################################
# GENERIC ACOUSTIC VARIABLES
################################################################################
################################################################################
#' Calculate the acoustic wavenumber (\eqn{k}) for a given frequency and sound 
#' speed.
#' 
#' @description
#' Calculates the acoustic wavenumber (\eqn{k}) for a given frequency and sound 
#' speed in water. The wavenumber is defined as:
#' \deqn{k = \frac{2\pi f}{c}}
#' where \eqn{f} is the frequency (Hz) and \eqn{c} is the sound speed 
#' (\eqn{ms^{-1}}).
#' The wavenumber describes the spatial frequency of a sound wave and is 
#' fundamental in acoustic calculations.
#'
#' @param sound_speed Sound speed (c, \eqn{m~s^{-1}})
#' @param frequency Frequency (f, Hz)
#' @return
#' Calculates the acoustic wavenumber (k) based on the sound speed of water.
#' @rdname wavenumber
#' @export
wavenumber <- function(frequency, sound_speed) 2 * pi * frequency / sound_speed
################################################################################
#' Calculates the linear backscattering coefficient (\eqn{\sigma_\text{bs}}) 
#' from the linear scattering length/coefficient, \eqn{f_\text{bs}}.
#' @param f_bs Linear scattering length (m), or related expression
#' @return
#' Returns the linear backscattering coefficient that can be then converted
#' into TS.
#' @keywords internal
#' @noRd
.sigma_bs <- function(f_bs) abs(f_bs) * abs(f_bs)
################################################################################
#' Convert between logarithmic (dB) and linear domains for backscatter values.
#'
#' @description
#' The `linear` function converts a value from the logarithmic (dB) domain to 
#' the linear domain, while the `db` function converts a value from the linear 
#' domain to the logarithmic (dB) domain. These are commonly used for target 
#' strength (TS) and backscattering coefficient (\eqn{\sigma_{bs}}) conversions.
#'
#' The conversions are defined as:
#' \deqn{\text{linear}(x) = c^{x / c}}
#' \deqn{\text{db}(x) = c \log_c(x)}
#' where \eqn{c} is the coefficient (default 10).
#'
#' @param value Numeric value to convert. For `linear`, this is a logarithmic 
#' value (e.g., dB TS); for `db`, this is a linear value (e.g., 
#' \eqn{\sigma_{bs}}).
#' @param coefficient Optional. Numeric coefficient (base) for the logarithm. 
#' Default is 10.
#' @return
#' For `linear`, returns the value converted to the linear domain. For `db`, 
#' returns the value converted to the logarithmic (dB) domain.
#' @rdname linear
#' @export
linear <- function(value, coefficient = 10) {
  coefficient^(value / coefficient)
}
#' @rdname linear
#' @export
db <- function(value, coefficient = 10) {
  coefficient * log(value, base = coefficient)
}
################################################################################
################################################################################
# REFLECTIVITY VARIABLES, COEFFICIENTS, AND EQUATIONS
################################################################################
################################################################################
#' Plane wave/plane interface reflection coefficient
#' @param interface1 Dataframe object containing density (kg/m^3) and sound
#' speed (m/s) values for a boundary/interface (1)
#' @param interface2 Dataframe object containing density (kg/m^3) and sound
#' speed (m/s) values for a boundary/interface (2)
#' @param mode Two options: coefficient calculation for "DWBA" and "KRM"
#' @noRd
reflection_coefficient <- function(interface1, interface2, mode = "DWBA") {
  # Calculate acoustic impedance of the first interface ========================
  Z1 <- interface1$density * interface1$sound_speed
  # Calculate acoustic impedance of the second interface =======================
  Z2 <- interface2$density * interface2$sound_speed
  # Calculate the reflection coefficient =======================================
  R <- (Z2 - Z1) / (Z2 + Z1)
  # Output =====================================================================
  return(R)
}
################################################################################
#' Plane wave/plane interface transmission coefficient
#' @param interface1 Dataframe object containing density (kg/m^3) and sound
#' speed (m/s) values for a boundary/interface (1)
#' @param interface2 Dataframe object containing density (kg/m^3) and sound
#' speed (m/s) values for a boundary/interface (2)
#' @param mode Two options: coefficient calculation for "DWBA" and "KRM"
#' @return Pressure-amplitude transmission coefficient at normal incidence.
#' @export
transmission_coefficient <- function(interface1, interface2, mode = "DWBA") {
  Z1 <- interface1$density * interface1$sound_speed
  Z2 <- interface2$density * interface2$sound_speed
  2 * Z2 / (Z1 + Z2)
}
################################################################################
#' Calculate she compressibility (\eqn{\kappa}) of a scattering 
#' boundary/interface.
#' 
#' @description
#' Calculates the compressibility contrast (\eqn{\kappa}) between a scattering 
#' interface and the surrounding medium. Compressibility is defined as:

#' \deqn{
#'  K = \frac{1}{\rho c^2}
#' }
#' where \eqn{\rho} is density (\eqn{kg~m^{-3}}) and \eqn{c} is sound speed 
#' (\eqn{m~s^{-1}}).
#'
#' The compressibility contrast is then:
#' \deqn{
#'  \kappa = \frac{K_2 - K_1}{K_1}
#' }
#'
#' where \eqn{K_1} is the compressibility of the medium and \eqn{K_2} is that 
#' of the target interface.
#'
#' @param medium Dataframe object containing density (\eqn{kgm^{-3}}) and 
#' sound speed (\eqn{ms^{-1}}) values for a fluid medium external to a 
#' scattering interface (e.g., seawater).
#' @param target Dataframe object containing density (\eqn{kgm^{-3}}) and 
#' sound speed (\eqn{ms^{-1}}) values for a target boundary.
#' 
#' @return Compressibility contrast (\eqn{\kappa}), dimensionless.
#' @export
compressibility <- function(medium, target) {
  # Calculate acoustic compressibility of the first interface ==================
  K1 <- (medium$density * medium$sound_speed^2)^(-1)
  # Calculate acoustic compressibility of the second interface =================
  K2 <- (target$density * target$sound_speed^2)^(-1)
  # Calculate the reflection coefficient =======================================
  gamma_kappa <- (K2 - K1) / K1
  # Output =====================================================================
  gamma_kappa
}
################################################################################
#' Calculate the density contrast (\eqn{\rho}) of a scattering boundary
#'
#' @param medium Dataframe object containing density (\eqn{kgm^{-3}}) and
#' sound speed (\eqn{ms^{-1}}) values for the external medium.
#' @param target Dataframe object containing density (\eqn{kgm^{-3}}) and
#' sound speed (\eqn{ms^{-1}}) values for the target boundary.
#' @return Density contrast, defined as
#'   \eqn{(\rho_{target} - \rho_{medium})/\rho_{target}}.
#' @export
rho <- function(medium, target) {
  (target$density - medium$density) / target$density
}
################################################################################
################################################################################
# ELASTICITY CALCULATIONS AND EQUATIONS
################################################################################
################################################################################
#' @title Calculate the Poisson's ratio (\eqn{\nu})
#' 
#' @description
#' Calculates Poisson's ratio (\eqn{\nu}) from two of the three other elastic 
#' moduli: bulk modulus (K), Young's modulus (E), or shear modulus (G). 
#' Assumes 3D material properties.
#'
#' The relationships used are:
#' \deqn{\nu = \frac{E}{2G} - 1}
#' \deqn{\nu = \frac{3K - 2G}{2(3K + G)}}
#' \deqn{\nu = \frac{3K - E}{6K}}
#' 
#' @param K Bulk modulus (K, Pa).
#' @param E Young's modulus (E, Pa).
#' @param G Shear modulus (Pa).
#' @return Poisson's ratio (\eqn{\nu}), dimensionless.
#' 
#' @keywords elastic
#' @rdname pois
#' @encoding UTF-8
#' @export
pois <- function(K = NULL, E = NULL, G = NULL) {
  # Check inputs ===============================================================
  inputs <- list(K = K, E = E, G = G)
  provided <- !vapply(inputs, is.null, logical(1))
  # Checksum ===================================================================
  if (sum(provided) < 2) {
    stop(
      paste0(
        "At least two elasticity moduli values are required to ",
        "calculate Poisson's ratio."
      )
    )
  }
  # Use first available combination ============================================
  if (provided["E"] && provided["G"]) {
    return(E / (2 * G) - 1)
  } else if (provided["K"] && provided["G"]) {
    return(3 * K - 2 * G) / (2 * (3 * K + G))
  } else if (provided["K"] && provided["E"]) {
    return((3 * K - E) / (6 * K))
  }
}
################################################################################
#' Calculate the bulk modulus (K).
#' 
#' @description
#' Calculates the bulk modulus (K) from two of the three other elastic moduli:
#' Young's modulus (E), shear modulus (G), or Poisson's ratio (\eqn{\nu}).
#' Assumes 3D material properties.
#'
#' The relationships used are:
#' \deqn{K = \frac{E G}{3(3G - E)}}
#' \deqn{K = \frac{2G(1 + \nu)}{3(1 - 2\nu)}}
#' \deqn{K = \frac{E}{3(1 - 2\nu)}}
#' 
#' @param E Young's modulus (Pa).
#' @param G Shear modulus (Pa).
#' @param nu Poisson's ratio (Dimensionless).
#' 
#' @return Bulk modulus (K, Pa).
#' 
#' @encoding UTF-8
#' @keywords elastic
#' @rdname bulk
#' @export
bulk <- function(E = NULL, G = NULL, nu = NULL) {
  # Check inputs ===============================================================
  inputs <- list(E = E, G = G, nu = nu)
  provided <- !vapply(inputs, is.null, logical(1))
  # Checksum ===================================================================
  if (sum(provided) < 2) {
    stop(
      paste0(
        "At least two elasticity moduli values are required to ",
        "calculate the bulk modulus."
      )
    )
  }
  # Use first available combination ============================================
  if (provided["E"] && provided["G"]) {
    return(E * G / (3 * (3 * G - E)))
  } else if (provided["G"] && provided["nu"]) {
    return(2 * G * (1 + nu) / (3 * (1 - 2 * nu)))
  } else if (provided["E"] && provided["nu"]) {
    return(E / (3 * (1 - 2 * nu)))
  }
}
################################################################################
#' Calculate Young's modulus (E).
#' 
#' @description
#' Calculates Young's modulus (E) from two of the three other elastic moduli:
#' bulk modulus (K), shear modulus (G), or Poisson's ratio (\eqn{\nu}).
#' Assumes 3D material properties.
#'
#' The relationships used are:
#' \deqn{E = \frac{9KG}{3K + G}}
#' \deqn{E = 3K(1 - 2\nu)}
#' \deqn{E = 2G(1 + \nu)}
#' 
#' @param K Bulk modulus (Pa).
#' @param G Shear modulus (Pa).
#' @param nu Poisson's ratio (Dimensionless).
#' 
#' @return Young's modulus (E, Pa).
#' 
#' @encoding UTF-8
#' @keywords elastic
#' @rdname young
#' @export
young <- function(K = NULL, G = NULL, nu = NULL) {
  # Check inputs ===============================================================
  inputs <- list(K = K, G = G, nu = nu)
  provided <- !vapply(inputs, is.null, logical(1))
  # Checksum ===================================================================
  if (sum(provided) < 2) {
    stop(
      paste0(
        "At least two elasticity moduli values are required to ",
        "calculate Young's modulus."
      )
    )
  }
  # Use first available combination ============================================
  if (provided["K"] && provided["G"]) {
    return(9 * K * G / (3 * K + G))
  } else if (provided["K"] && provided["nu"]) {
    return(3 * K * (1 - 2 * nu))
  } else if (provided["G"] && provided["nu"]) {
    return(2 * G * (1 + nu))
  }
}
################################################################################
#' Calculate the shear modulus (G)
#' 
#' @description
#' #' Calculates the shear modulus (G) from two of the three other elastic 
#' moduli: bulk modulus (K), Young's modulus (E), or Poisson's ratio 
#' (\eqn{\nu}). Assumes 3D material properties.
#'
#' The relationships used are:
#' \deqn{G = \frac{3KE}{9K - E}}
#' \deqn{G = \frac{3K(1 - 2\nu)}{2(1 + \nu)}}
#' \deqn{G = \frac{E}{2(1 + \nu)}}
#' 
#' @param K Bulk modulus (Pa).
#' @param E Young's modulus (Pa).
#' @param nu Poisson's ratio (Dimensionless).
#' 
#' @return Shear modulus (G, Pa).
#' 
#' @encoding UTF-8
#' @keywords elastic
#' @rdname shear
#' @export
shear <- function(K = NULL, E = NULL, nu = NULL) {
  # Check inputs ===============================================================
  inputs <- list(K = K, E = E, nu = nu)
  provided <- !vapply(inputs, is.null, logical(1))
  # Checksum ===================================================================
  if (sum(provided) < 2) {
    stop(
      paste0(
        "At least two elasticity moduli values are required to ",
        "calculate the shear modulus."
      )
    )
  }
  # Use first available combination ============================================
  if (provided["K"] && provided["E"]) {
    return(3 * K * E / (9 * K - E))
  } else if (provided["K"] && provided["nu"]) {
    return(3 * K * (1 - 2 * nu) / (2 * (1 + nu)))
  } else if (provided["E"] && provided["nu"]) {
    return(E / (2 * (1 + nu)))
  }
}
################################################################################
#' @encoding UTF-8
#' @title Calculate Lam&eacute;'s first parameter (\eqn{\lambda})
#' 
#' @description
#' Calculates Lam&eacute;'s first parameter (\eqn{\lambda}) from two of the 
#' four other elastic moduli: bulk modulus (K), Young's modulus (E), shear 
#' modulus (G), or Poisson's ratio (\eqn{\nu}). Assumes 3D material properties.
#'
#' The relationships used are:
#' \deqn{\lambda = K - \frac{2G}{3}}
#' \deqn{\lambda = \frac{E\nu}{(1 + \nu)(1 - 2\nu)}}
#' \deqn{\lambda = \frac{2G\nu}{1 - 2\nu}}
#' \deqn{\lambda = \frac{3K\nu}{1 + \nu}}
#' \deqn{\lambda = \frac{3K(3K - E)}{9K - E}}
#' \deqn{\lambda = \frac{G(E - 2G)}{3G - E}}
#' 
#' @param K Bulk modulus (Pa).
#' @param E Young's modulus (Pa).
#' @param G Shear modulus (Pa).
#' @param nu Poisson's ratio (Dimensionless).
#' @return Lam&eacute;'s first parameter (\eqn{\lambda}, Pa).
#'
#' @keywords elastic
#' @rdname lame
#' @export
lame <- function(K = NULL, E = NULL, G = NULL, nu = NULL) {
  # Check inputs ===============================================================
  inputs <- list(K = K, E = E, G = G, nu = nu)
  provided <- !vapply(inputs, is.null, logical(1))
  # Checksum ===================================================================
  if (sum(provided) < 2) {
    stop(
      paste0(
        "At least two elasticity moduli values are required to ",
        "calculate the Lam\U00E9 parameter."
      )
    )
  }
  # Use first available combination ============================================
  if (provided["K"] && provided["G"]) {
    return(K - 2 * G / 3)
  } else if (provided["E"] && provided["nu"]) {
    return(E * nu / ((1 + nu) * (1 - 2 * nu)))
  } else if (provided["G"] && provided["nu"]) {
    return(2 * G * nu / (1 - 2 * nu))
  } else if (provided["K"] && provided["nu"]) {
    return(3 * K * nu / (1 + nu))
  } else if (provided["K"] && provided["E"]) {
    return(3 * K * (3 * K - E) / (9 * K - E))
  } else if (provided["E"] && provided["G"]) {
    return(G * (E - 2 * G) / (3 * G - E))
  }
}
################################################################################
#' Wrapper function to model acoustic target strength
#' @param object Scatterer-class object.
#' @param frequency Frequency (Hz).
#' @param model Model name. Available models currently include
#'   \code{"dwba"} (\code{\link{DWBA}}),
#'   \code{"bbfm"} (\code{\link{BBFM}}),
#'   \code{"pcdwba"} (\code{\link{PCDWBA}}),
#'   \code{"sdwba"} (\code{\link{SDWBA}}),
#'   \code{"fcms"} (\code{\link{FCMS}}),
#'   \code{"bcms"} (\code{\link{BCMS}}),
#'   \code{"ecms"} (\code{\link{ECMS}}),
#'   \code{"hpa"} (\code{\link{HPA}}),
#'   \code{"krm"} (\code{\link{KRM}}),
#'   \code{"psms"} (\code{\link{PSMS}}),
#'   \code{"sphms"} (\code{\link{SPHMS}}),
#'   \code{"essms"} (\code{\link{ESSMS}}),
#'   \code{"vesms"} (\code{\link{VESMS}}),
#'   \code{"trcm"} (\code{\link{TRCM}}), and
#'   \code{"calibration"} / \code{"soems"} (\code{\link{SOEMS}}).
#' @param verbose Prints current procedural step occurring from model
#' initialization to calculating TS. Defaults to FALSE.
#' @param ... Additional optional model inputs/parameters.
#' @details
#' This is the main high-level entry point for running target-strength models in
#' \code{acousticTS}. The supplied scatterer object is checked, the requested
#' model or models are initialized, and the resulting outputs are stored back on
#' the same object.
#'
#' The available model families span exact modal-series solutions and
#' approximation-based solutions. Readers should consult the model-specific help
#' topics for the physical assumptions, valid object types, boundary-condition
#' options, and model-specific arguments used by each implementation:
#'
#' \itemize{
#'   \item \code{\link{DWBA}} for the distorted-wave Born approximation
#'   applied to weakly scattering fluid-like bodies.
#'   \item \code{\link{BBFM}} for a composite flesh-plus-backbone model that
#'   combines a DWBA flesh term with an elastic-cylinder backbone term.
#'   \item \code{\link{PCDWBA}} for the phase-compensated distorted-wave Born
#'   approximation applied to elongated fluid-like bodies.
#'   \item \code{\link{SDWBA}} for the stochastic distorted-wave Born
#'   approximation.
#'   \item \code{\link{FCMS}} for the finite cylinder modal series solution.
#'   \item \code{\link{BCMS}} for the bent-cylinder modal series solution.
#'   \item \code{\link{ECMS}} for the elastic-cylinder modal series solution.
#'   \item \code{\link{HPA}} for the high-pass approximation family.
#'   \item \code{\link{KRM}} for the Kirchhoff-Ray Mode model.
#'   \item \code{\link{PSMS}} for the prolate spheroidal modal series
#'   solution.
#'   \item \code{\link{SPHMS}} for the spherical modal series solution.
#'   \item \code{\link{ESSMS}} for the elastic-shelled spherical modal series
#'   solution.
#'   \item \code{\link{VESMS}} for the viscous-elastic spherical scattering
#'   model applied to gas-filled elastic shells with an external viscous layer.
#'   \item \code{\link{TRCM}} for the two-ray cylindrical model.
#'   \item \code{\link{SOEMS}} for the solid elastic calibration-sphere model,
#'   accessed through \code{"calibration"} or \code{"soems"}.
#' }
#'
#' Model-specific inputs are passed through \code{...}. For example, some models
#' require a \code{boundary} argument, \code{HPA} uses a \code{method}
#' argument, and several models expose additional numerical controls. The legacy
#' curved-entry wrappers \code{"dwba_curved"} and \code{"sdwba_curved"} are
#' deprecated; apply \code{\link{brake}} to the scatterer first, then run
#' \code{"dwba"} or \code{"sdwba"} on the curved object.
#' @seealso
#' \code{\link{DWBA}}, \code{\link{BBFM}}, \code{\link{PCDWBA}}, \code{\link{SDWBA}}, \code{\link{FCMS}},
#' \code{\link{BCMS}}, \code{\link{ECMS}}, \code{\link{HPA}}, \code{\link{KRM}}, \code{\link{PSMS}},
#' \code{\link{SPHMS}}, \code{\link{ESSMS}}, \code{\link{VESMS}}, \code{\link{TRCM}},
#' \code{\link{SOEMS}}
#' @export
target_strength <- function(object, frequency, model, verbose = FALSE, ...) {
  # Validate inputs ============================================================
  if (missing(object)) stop("Scattering object ('object') is required")
  if (missing(frequency)) stop("Frequency (Hz) ('frequency') is required")
  if (missing(model)) stop("Target strength model ('model') is required")

  # Store the object in a variable that won't conflict with model internals ====
  target_object <- object

  # Capture all arguments including ... ========================================
  arg_pull <- list(object = target_object, frequency = frequency, ...)

  # Handle model names (convert to uppercase for consistency) ==================
  model <- tolower(model)
  ts_model <- gsub("(_.*)", "\\L\\1", paste0(toupper(model)), perl = TRUE)
  ts_model <- ifelse(ts_model %in% c(
    "CALIBRATION",
    "HIGH_pass_stanton"
  ),
  tolower(ts_model),
  ts_model
  )

  # Initialize objects to input model parameters ==============================
  idx <- 1
  repeat {
    if (idx > length(model)) {
      break
    }

    # Pull correct formal arguments ==========================================
    model_name <- paste0(model[idx], "_initialize")

    # Check if initialization function exists ================================
    if (!exists(model_name)) {
      stop(
        "Initialization function ",
        model_name,
        " not found for model ",
        model[idx]
      )
    }

    # Filter out inappropriate parameters ====================================
    true_args <- .filter_shape_args(model_name, arg_pull)

    # Initialize ==============================================================
    object_copy <- do.call(model_name, true_args)

    # Store model parameters and results =====================================
    methods::slot(
      target_object,
      "model_parameters"
    )[ts_model[idx]] <- extract(
      object_copy, "model_parameters"
    )[ts_model[idx]]
    methods::slot(
      target_object,
      "model"
    )[ts_model[idx]] <- extract(
      object_copy, "model"
    )[ts_model[idx]]

    if (verbose) {
      cat(
        toupper(model[idx]), "model for", paste0(
          class(target_object), "-object: ",
          extract(target_object, "metadata")$ID
        ),
        "initialized.\n\n"
      )
    }

    idx <- idx + 1
  }

  # Run the models =============================================================
  idx <- 1
  repeat {
    if (idx > length(model)) {
      break
    }

    if (verbose) {
      cat(
        "Beginning TS modeling via", toupper(model[idx]),
        "model for", paste0(
          class(target_object), "-object: ",
          extract(target_object, "metadata")$ID
        ), "\n"
      )
    }

    # Check if model function exists ===========================================
    if (!exists(ts_model[idx])) {
      stop("Model function ", ts_model[idx], " not found")
    }

    # Calculate modeled TS using do.call instead of eval(parse(...)) ===========
    target_object <- do.call(ts_model[idx], list(object = target_object))

    if (verbose) {
      cat(toupper(model[idx]), "TS model predictions for", paste0(
        class(target_object), "-object: ",
        extract(target_object, "metadata")$ID
      ), "complete.\n\n")
    }

    idx <- idx + 1
  }

  # Output object ==============================================================
  return(target_object)
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
  k1 <- wavenumber(frequency, sound_speed_sw)
  k3 <- wavenumber(frequency, sound_speed_fluid)
  kL <- wavenumber(frequency, sound_speed_longitudinal)
  kT <- wavenumber(frequency, sound_speed_transversal)

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
################################################################################
#' Solve Expansion Coefficients for a Liquid-Filled Spheroidal Scatterer
#'
#' @description
#' Computes the modal expansion coefficients \eqn{A_{mn}} for acoustic 
#' scattering from a liquid-filled prolate spheroidal body using a truncated 
#' singular value decomposition (SVD) pseudoinverse approach.
#'
#' @details
#' This function solves the linear system arising from matching boundary
#' conditions at the surface of a liquid-filled spheroidal scatterer. The
#' expansion coefficients  \eqn{A_{mn}} relate the scattered field to the
#' incident field through the spheroidal wave function expansion.
#' @keywords internal
#' @noRd
.psms_adaptive_n_integration <- function(
  chi_sw,
  chi_body,
  m_max,
  n_max,
  precision
) {
  chi_max <- pmax(abs(chi_sw), abs(chi_body))
  step_n <- 8L
  floor_n <- if (identical(precision, "quad")) 32L else 24L
  cap_n <- 96L

  # The overlap quadrature must resolve the same angular content that drives the
  # retained PSMS modal system. Use the hard truncation ceilings as the main
  # difficulty scale, with a smaller reduced-frequency bonus for the sharper
  # high-chi oscillations that appear in the angular products.
  modal_difficulty <- 0.75 * n_max + 0.25 * m_max
  chi_bonus <- if (identical(precision, "quad")) {
    0.10 * sqrt(pmax(chi_max, 1))
  } else {
    0.05 * sqrt(pmax(chi_max, 1))
  }

  n_integration <- step_n * ceiling((modal_difficulty + chi_bonus) / step_n)
  n_integration <- pmax(floor_n, pmin(cap_n, n_integration))

  as.integer(n_integration)
}

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
  quad_pts <- gauss_legendre(n = n_integration, a = -1, b = 1)
  # Calculate the linear scattering coefficient, fbs ===========================
  prolate_spheroid_fbs(
    acoustics, body, medium, quad_pts, precision, boundary_method, adaptive,
    FALSE
  )
}

.prolate_spheroidal_kernels_adaptive <- function(
  acoustics,
  body,
  medium,
  boundary_method,
  precision = "double",
  adaptive = TRUE
) {
  n_by_freq <- .psms_adaptive_n_integration(
    chi_sw = acoustics$chi_sw,
    chi_body = acoustics$chi_body,
    m_max = acoustics$m_max,
    n_max = acoustics$n_max,
    precision = precision
  )
  f_bs <- complex(length = nrow(acoustics))

  for (n_integration_i in sort(unique(n_by_freq))) {
    idx <- which(n_by_freq == n_integration_i)
    f_bs[idx] <- .prolate_spheroidal_kernels_fixed(
      acoustics = acoustics[idx, , drop = FALSE],
      body = body,
      medium = medium,
      boundary_method = boundary_method,
      n_integration = n_integration_i,
      precision = precision,
      adaptive = adaptive
    )
  }

  attr(f_bs, "n_integration") <- n_by_freq
  f_bs
}

prolate_spheroidal_kernels <- function(
  acoustics,
  body,
  medium,
  boundary_method,
  n_integration = 96,
  precision = "double",
  adaptive = FALSE
) {
  if (isTRUE(adaptive) && identical(boundary_method, "Amn_fluid")) {
    return(
      .prolate_spheroidal_kernels_adaptive(
        acoustics = acoustics,
        body = body,
        medium = medium,
        boundary_method = boundary_method,
        precision = precision,
        adaptive = adaptive
      )
    )
  }

  if (is.null(n_integration) || (length(n_integration) == 1 && is.na(n_integration))) {
    n_integration <- 96L
  }

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
