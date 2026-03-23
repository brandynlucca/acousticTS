################################################################################
################################################################################
# INITIALIZATION AND MATERIAL-PROPERTY UTILITIES
################################################################################
################################################################################
#' Derive density and sound speed contrasts from absolute properties
#'
#' Prefer explicit contrasts `g`/`h` when present; otherwise compute them
#' relative to the provided medium properties when absolute density and/or
#' sound speed are available.
#'
#' @param body List-like object containing optional `g`, `h`, `density`,
#'   and `sound_speed` entries.
#' @param medium_sound_speed Surrounding-medium sound speed.
#' @param medium_density Surrounding-medium density.
#' @return Named list with elements `g` and `h` (may be NA if insufficient
#'   information was provided).
#' @keywords internal
#' @noRd
.derive_contrasts <- function(body, medium_sound_speed, medium_density) {
  h <- if (!is.null(body$h)) {
    body$h
  } else if (!is.null(body$sound_speed)) {
    body$sound_speed / medium_sound_speed
  } else {
    NA
  }

  g <- if (!is.null(body$g)) {
    body$g
  } else if (!is.null(body$density)) {
    body$density / medium_density
  } else {
    NA
  }

  list(g = g, h = h)
}

#' Extract common initialization components
#'
#' Internal helper to extract shape, body, and medium parameters from scatterer
#' objects during model initialization.
#'
#' @param object Scatterer object
#' @param sound_speed_sw Sound speed in seawater (m/s)
#' @param density_sw Density of seawater (kg/m^3)
#' @return List with shape, body, and medium components
#' @keywords internal
#' @noRd
.init_common <- function(object, sound_speed_sw = 1500, density_sw = 1026) {
  list(
    shape = extract(object, "shape_parameters"),
    body = extract(object, "body"),
    medium = .init_medium_params(sound_speed_sw, density_sw)
  )
}

#' Build surrounding-medium parameter frame
#' @param sound_speed_sw Surrounding-medium sound speed.
#' @param density_sw Surrounding-medium density.
#' @return Single-row data frame with medium properties.
#' @keywords internal
#' @noRd
.init_medium_params <- function(sound_speed_sw = 1500, density_sw = 1026) {
  data.frame(
    sound_speed = sound_speed_sw,
    density = density_sw
  )
}

#' Build acoustic recipe with named wavenumber columns
#'
#' @param frequency Frequency vector in Hz.
#' @param ... Named sound-speed vectors or scalars. Each name becomes a
#'   wavenumber column in the returned data frame.
#' @return Data frame containing `frequency` and named wavenumber columns.
#' @keywords internal
#' @noRd
.init_acoustics_df <- function(frequency, ...) {
  sound_speeds <- list(...)
  acoustics <- data.frame(frequency = frequency)

  if (length(sound_speeds) == 0) {
    return(acoustics)
  }

  for (nm in names(sound_speeds)) {
    acoustics[[nm]] <- wavenumber(frequency, sound_speeds[[nm]])
  }

  acoustics
}

#' Complete contrast and absolute material properties
#'
#' @param component List-like component containing optional `g`, `h`, `density`,
#'   and `sound_speed`.
#' @param medium_sound_speed Surrounding-medium sound speed.
#' @param medium_density Surrounding-medium density.
#' @return Updated component with both contrast and absolute properties filled
#'   when possible.
#' @keywords internal
#' @noRd
.complete_material_props <- function(component,
                                     medium_sound_speed,
                                     medium_density) {
  if (is.null(component$g) && !is.null(component$density)) {
    component$g <- component$density / medium_density
  }
  if (is.null(component$h) && !is.null(component$sound_speed)) {
    component$h <- component$sound_speed / medium_sound_speed
  }
  if (is.null(component$density) && !is.null(component$g)) {
    component$density <- component$g * medium_density
  }
  if (is.null(component$sound_speed) && !is.null(component$h)) {
    component$sound_speed <- component$h * medium_sound_speed
  }

  component
}

#' Fill body contrast properties from absolute values when needed
#'
#' @param body List-like body component.
#' @param medium_sound_speed Surrounding-medium sound speed.
#' @param medium_density Surrounding-medium density.
#' @return Updated body component with `g` and `h` populated when possible.
#' @keywords internal
#' @noRd
.hydrate_contrasts <- function(body, medium_sound_speed, medium_density) {
  contrasts <- .derive_contrasts(body, medium_sound_speed, medium_density)
  body$h <- contrasts$h
  body$g <- contrasts$g
  body
}

#' Store model parameters and initialize the matching results slot
#'
#' @param object Scatterer object.
#' @param model_name Name of the model slot to write.
#' @param frequency Frequency vector used for the model run.
#' @param model_parameters List to store under `object@model_parameters`.
#' @param result_cols Character vector of result-column names to initialize
#'   with `NA_real_`.
#' @return Updated scatterer object.
#' @keywords internal
#' @noRd
.init_model_slots <- function(object,
                              model_name,
                              frequency,
                              model_parameters,
                              result_cols = c("sigma_bs")) {
  methods::slot(object, "model_parameters")[[model_name]] <- model_parameters

  model_frame <- data.frame(frequency = frequency)
  for (col in result_cols) {
    model_frame[[col]] <- rep(NA_real_, length(frequency))
  }

  methods::slot(object, "model")[[model_name]] <- model_frame
  object
}

#' Calculate all relevant wavenumbers
#'
#' Internal helper to calculate wavenumbers for multiple sound speeds.
#'
#' @param frequency Frequency vector (Hz)
#' @param sound_speeds Named list of sound speeds (m/s)
#' @return Named list of wavenumber vectors
#' @keywords internal
#' @noRd
.calc_wavenumbers <- function(frequency, sound_speeds) {
  stats::setNames(
    lapply(sound_speeds, function(c) wavenumber(frequency, c)),
    names(sound_speeds)
  )
}

#' Extract material properties with fallback logic
#'
#' Internal helper to extract material properties from ESS components with
#' contrast-based fallback calculations.
#'
#' @param component Material component (shell, fluid, etc.)
#' @param ref_sound_speed Reference sound speed for contrast conversion
#' @param ref_density Reference density for contrast conversion
#' @return List of material properties
#' @keywords internal
#' @noRd
.extract_material_props <- function(component, ref_sound_speed = NULL,
                                    ref_density = NULL) {
  get_prop <- function(direct, contrast, ref_val, default = NA) {
    if (!is.null(component[[direct]])) {
      return(component[[direct]])
    }
    if (!is.null(component[[contrast]]) && !is.null(ref_val)) {
      return(component[[contrast]] * ref_val)
    }
    return(default)
  }

  list(
    sound_speed = get_prop("sound_speed", "h", ref_sound_speed),
    density = get_prop("density", "g", ref_density),
    nu = component$nu %||% NA,
    K = component$K %||% NA,
    E = component$E %||% NA,
    G = component$G %||% NA,
    lambda = component$lambda %||% NA
  )
}
