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
  # Validate the available elastic moduli before selecting a formula ==========
  .validate_elastic_inputs(
    K = K,
    E = E,
    G = G,
    nu = nu,
    param_name = "the Lam\U00E9 parameter"
  )
  # Resolve the first supported parameter combination ==========================
  combo <- .lame_input_combo(K = K, E = E, G = G, nu = nu)

  switch(combo,
    K_G = K - 2 * G / 3,
    E_nu = E * nu / ((1 + nu) * (1 - 2 * nu)),
    G_nu = 2 * G * nu / (1 - 2 * nu),
    K_nu = 3 * K * nu / (1 + nu),
    K_E = 3 * K * (3 * K - E) / (9 * K - E),
    E_G = G * (E - 2 * G) / (3 * G - E)
  )
}

#' Resolve the precedence-ordered parameter pair used by lame()
#' @keywords internal
#' @noRd
.lame_input_combo <- function(K = NULL, E = NULL, G = NULL, nu = NULL) {
  # Record which elastic moduli are available ==================================
  provided <- !vapply(
    list(K = K, E = E, G = G, nu = nu),
    is.null,
    logical(1)
  )

  # Apply the historical precedence order used by the public helper ===========
  if (provided["K"] && provided["G"]) {
    return("K_G")
  }
  if (provided["E"] && provided["nu"]) {
    return("E_nu")
  }
  if (provided["G"] && provided["nu"]) {
    return("G_nu")
  }
  if (provided["K"] && provided["nu"]) {
    return("K_nu")
  }
  if (provided["K"] && provided["E"]) {
    return("K_E")
  }

  "E_G"
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
#' @param model_args Optional named list of per-model argument bundles. Each
#'   list name should match one of the requested model names
#'   case-insensitively, and each value should be either a named list or a
#'   named atomic vector of arguments to apply only to that model. When the
#'   same argument is supplied both through `...` and through
#'   `model_args[[model_name]]`, the model-specific entry takes precedence.
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
#'   \item \code{\link{TMM}} for the single-target transition-matrix family.
#'   \item \code{\link{TRCM}} for the two-ray cylindrical model.
#'   \item \code{\link{SOEMS}} for the solid elastic calibration-sphere model,
#'   accessed through \code{"calibration"} or \code{"soems"}.
#' }
#'
#' Model-specific inputs are usually passed through \code{...}. For example,
#' some models require a \code{boundary} argument, \code{HPA} uses a
#' \code{method} argument, and several models expose additional numerical
#' controls. When several models are requested together, shared arguments may be
#' supplied through \code{...} and per-model overrides may be supplied through
#' \code{model_args}. This is useful when different models should share the same
#' seawater properties but only one of them needs an extra stochastic or
#' numerical control:
#'
#' \preformatted{
#' target_strength(
#'   object,
#'   frequency,
#'   model = c("dwba", "sdwba"),
#'   density_sw = 1026,
#'   sound_speed_sw = 1478,
#'   model_args = list(
#'     sdwba = list(phase_sd_init = 0.77)
#'   )
#' )
#' }
#'
#' The legacy curved-entry wrappers \code{"dwba_curved"} and
#' \code{"sdwba_curved"} are deprecated; apply \code{\link{brake}} to the
#' scatterer first, then run \code{"dwba"} or \code{"sdwba"} on the curved
#' object. Model names are normalized internally, so case-insensitive inputs
#' such as \code{"DWBA"} and \code{"dwba"} resolve to the same family.
#' @return The input scatterer object with requested model parameters, model
#'   outputs, and target strength results stored in its model slots.
#' @seealso
#' \code{\link{DWBA}}, \code{\link{BBFM}}, \code{\link{PCDWBA}},
#' \code{\link{SDWBA}}, \code{\link{FCMS}}, \code{\link{BCMS}},
#' \code{\link{ECMS}}, \code{\link{HPA}}, \code{\link{KRM}}, \code{\link{PSMS}},
#' \code{\link{SPHMS}}, \code{\link{ESSMS}}, \code{\link{VESMS}},
#' \code{\link{TMM}}, \code{\link{TRCM}}, \code{\link{SOEMS}}
#' @export
target_strength <- function(object,
                            frequency,
                            model,
                            verbose = FALSE,
                            model_args = NULL,
                            ...) {
  # Validate inputs ============================================================
  if (missing(object)) stop("Scattering object ('object') is required")
  if (missing(frequency)) stop("Frequency (Hz) ('frequency') is required")
  if (missing(model)) stop("Target strength model ('model') is required")
  # Store the object in a variable that won't conflict with model internals ====
  target_object <- object
  # Capture all arguments including ... ========================================
  shared_args <- list(object = target_object, frequency = frequency, ...)
  # Resolve the requested target-strength model names ==========================
  model_names <- .resolve_target_strength_model_names(model)
  model <- model_names$model
  ts_model <- model_names$ts_model
  model_args <- .normalize_target_strength_model_args(model_args, model)
  # Initialize the requested model families ===================================
  target_object <- .initialize_target_strength_models(
    target_object = target_object,
    model = model,
    ts_model = ts_model,
    shared_args = shared_args,
    model_args = model_args,
    verbose = verbose
  )
  # Run the requested model families ===========================================
  .run_target_strength_models(
    target_object = target_object,
    model = model,
    ts_model = ts_model,
    verbose = verbose
  )
}

#' Normalize model names and exported solver names for target_strength()
#' @keywords internal
#' @noRd
.resolve_target_strength_model_names <- function(model) {
  # Normalize model labels to the naming used throughout the package ===========
  model <- tolower(model)
  ts_model <- gsub("(_.*)", "\\L\\1", paste0(toupper(model)), perl = TRUE)
  ts_model <- ifelse(
    ts_model %in% c("CALIBRATION", "HIGH_pass_stanton"),
    tolower(ts_model),
    ts_model
  )

  list(model = model, ts_model = ts_model)
}

#' Build a compact object label for target_strength() status messages
#' @keywords internal
#' @noRd
.target_strength_object_label <- function(object) {
  paste0(class(object), "-object: ", extract(object, "metadata")$ID)
}

#' Copy initialized model slots from a model-specific initialize() call
#' @keywords internal
#' @noRd
.store_target_strength_model_state <- function(target_object,
                                               object_copy,
                                               ts_model_name) {
  # Copy the initialized parameter recipe back onto the output object ==========
  methods::slot(target_object, "model_parameters")[ts_model_name] <-
    extract(object_copy, "model_parameters")[ts_model_name]
  methods::slot(target_object, "model")[ts_model_name] <-
    extract(object_copy, "model")[ts_model_name]

  target_object
}

#' Initialize all requested model families for target_strength()
#' @keywords internal
#' @noRd
.initialize_target_strength_models <- function(target_object,
                                               model,
                                               ts_model,
                                               shared_args,
                                               model_args,
                                               verbose = FALSE) {
  # Initialize one model family at a time =====================================
  for (idx in seq_along(model)) {
    model_name <- paste0(model[idx], "_initialize")

    if (!exists(model_name)) {
      stop(
        "Initialization function ",
        model_name,
        " not found for model ",
        model[idx]
      )
    }

    arg_pull <- .merge_target_strength_args(shared_args, model_args[[model[idx]]])
    true_args <- .filter_shape_args(model_name, arg_pull)
    object_copy <- do.call(model_name, true_args)
    target_object <- .store_target_strength_model_state(
      target_object = target_object,
      object_copy = object_copy,
      ts_model_name = ts_model[idx]
    )

    if (verbose) {
      cat(
        toupper(model[idx]), "model for",
        .target_strength_object_label(target_object),
        "initialized.\n\n"
      )
    }
  }

  target_object
}

#' Run all requested model families for target_strength()
#' @keywords internal
#' @noRd
.run_target_strength_models <- function(target_object,
                                        model,
                                        ts_model,
                                        verbose = FALSE) {
  # Evaluate the requested target-strength models sequentially ================
  for (idx in seq_along(model)) {
    if (verbose) {
      cat(
        "Beginning TS modeling via", toupper(model[idx]),
        "model for", .target_strength_object_label(target_object), "\n"
      )
    }

    if (!exists(ts_model[idx])) {
      stop("Model function ", ts_model[idx], " not found")
    }

    target_object <- do.call(ts_model[idx], list(object = target_object))

    if (verbose) {
      cat(
        toupper(model[idx]), "TS model predictions for",
        .target_strength_object_label(target_object),
        "complete.\n\n"
      )
    }
  }

  target_object
}

#' Normalize per-model argument bundles supplied to target_strength()
#' @keywords internal
#' @noRd
.normalize_target_strength_model_args <- function(model_args,
                                                  requested_models) {
  # Gather models to validate kwargs ===========================================
  requested_models <- unique(tolower(requested_models))
  # Revert to default values ===================================================
  if (is.null(model_args)) {
    out <- vector("list", length(requested_models))
    names(out) <- requested_models
    return(out)
  }
  # This internal helper only handles lists ====================================
  if (!is.list(model_args)) {
    stop("'model_args' must be a named list.", call. = FALSE)
  }
  # Validate specific kwargs ===================================================
  model_arg_names <- names(model_args)
  if (is.null(model_arg_names) || any(!nzchar(model_arg_names))) {
    stop("'model_args' must be a named list keyed by model name.",
      call. = FALSE
    )
  }
  # Check for argname conflicts ================================================
  normalized_names <- tolower(model_arg_names)
  if (anyDuplicated(normalized_names)) {
    stop(
      "'model_args' contains duplicate model entries after case normalization.",
      call. = FALSE
    )
  }
  # Superfluous kwargs check ===================================================
  unknown_models <- setdiff(normalized_names, requested_models)
  if (length(unknown_models) > 0) {
    stop(
      "'model_args' contains entries for model(s) not requested in 'model': ",
      paste(sprintf("'%s'", unknown_models), collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  # Push model-specific kwarg lists into vector and coerce into expected format
  out <- vector("list", length(requested_models))
  names(out) <- requested_models
  for (i in seq_along(model_args)) {
    out[[normalized_names[i]]] <- .coerce_target_strength_model_argset(
      model_args[[i]],
      model_name = normalized_names[i]
    )
  }

  out
}

#' Coerce one target_strength() per-model argument entry into a named list
#' @keywords internal
#' @noRd
.coerce_target_strength_model_argset <- function(x, model_name) {
  # Return empty set ===========================================================
  if (is.null(x)) {
    return(list())
  }
  # Compile and validate list-supplied kwargs ==================================
  if (is.list(x)) {
    if (!is.null(names(x)) && any(!nzchar(names(x)))) {
      stop(
        "The 'model_args' entry for model '", model_name,
        "' must use named arguments only.",
        call. = FALSE
      )
    }
    return(Filter(Negate(is.null), x))
  }
  # Non-list named kwarg converted to list =====================================
  if (is.atomic(x) && !is.null(names(x)) && all(nzchar(names(x)))) {
    return(as.list(x))
  }

  stop(
    "The 'model_args' entry for model '", model_name,
    "' must be either a named list or a named atomic vector.",
    call. = FALSE
  )
}

#' Merge shared and per-model argument bundles for target_strength()
#' @keywords internal
#' @noRd
.merge_target_strength_args <- function(shared_args,
                                        model_specific_args = NULL) {
  # Consolidate model-specific and list-grouped model kwargs ===================
  out <- shared_args
  if (is.null(model_specific_args) || length(model_specific_args) == 0) {
    return(out)
  }
  for (nm in names(model_specific_args)) {
    out[[nm]] <- model_specific_args[[nm]]
  }

  out
}
