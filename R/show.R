################################################################################
# Show functions
################################################################################
################################################################################
# Methods for "show(...)" for each scattering class object
################################################################################
#' Safe display helper: returns "NA" string for NULL/all-NA vectors
#' @noRd
.show_mean <- function(x) {
  # Return a printable NA token for missing inputs =============================
  if (is.null(x) || all(is.na(x))) return("NA")
  # Otherwise return the rounded mean ==========================================
  round(mean(x, na.rm = TRUE), 4)
}
# Centralized named-property rendering for show() methods
.show_property_block <- function(properties,
                                 label_map = character(),
                                 unit_map = character(),
                                 none_text = "   None specified") {
  # Return the fallback text for empty property blocks =========================
  if (length(properties) == 0) {
    return(none_text)
  }
  # Format each named property with labels and units ===========================
  prop_strings <- mapply(function(name, value) {
    clean_name <- if (name %in% names(label_map)) {
      label_map[[name]]
    } else {
      gsub("_", " ", name)
    }
    units <- if (name %in% names(unit_map)) unit_map[[name]] else ""
    paste0(clean_name, ": ", round(value, 4), units)
  }, names(properties), properties, SIMPLIFY = FALSE)
  # Join the formatted property lines ==========================================
  paste("   ", prop_strings, collapse = "\n ")
}

#' Show emit helper
#' @noRd
.show_emit <- function(lines) {
  # Emit one preformatted show() block =========================================
  cat(paste(lines, collapse = "\n"))
}

.show_section_lines <- function(title, ...) {
  # Prefix a block of lines with the section title =============================
  c(paste0(title, ":"), ...)
}

.show_header_lines <- function(object, descriptor, meta) {
  # Build the shared object header lines =======================================
  c(
    paste0(methods::is(object)[[1]], "-object"),
    paste0(" ", descriptor),
    paste0(" ID:", meta$ID)
  )
}

#' Show summary appending helper
#' @noRd
.show_summary_lines <- function(object, descriptor, meta, ...) {
  # Combine the header and body sections =======================================
  c(
    .show_header_lines(object, descriptor, meta),
    ...
  )
}

#' Show dimensions helper
#' @noRd
.show_dimension_lines <- function(section_name,
                                  length_value,
                                  units,
                                  n_segments = NULL,
                                  mean_radius = NULL,
                                  max_radius = NULL,
                                  segment_label = "cylinders") {
  # Build the primary length line ==============================================
  length_line <- paste0(
    " Length:",
    round(length_value, 3),
    " ",
    units,
    if (!is.null(n_segments)) {
      paste0("(n = ", n_segments, " ", segment_label, ")")
    } else {
      ""
    }
  )
  # Append the optional radius summaries =======================================
  lines <- c(paste0(section_name, ":"), length_line)
  if (!is.null(mean_radius)) {
    lines <- c(lines,
      paste0(" Mean radius:", round(mean_radius, 4), " ", units)
    )
  }
  if (!is.null(max_radius)) {
    lines <- c(lines,
      paste0(" Max radius:", round(max_radius, 4), " ", units)
    )
  }

  lines
}

#' Show orientation line helper
#' @noRd
.show_orientation_line <- function(label, theta, theta_units) {
  # Format one orientation line ================================================
  paste0(label, round(theta, 3), " ", theta_units)
}
.show_density_speed_line <- function(prefix,
                                     density,
                                     sound_speed,
                                     density_units = "kg m^-3",
                                     sound_speed_units = "m s^-1") {
  # Format one fluid-material summary line =====================================
  paste0(
    prefix,
    " Density: ", .show_mean(density),
    " ", density_units, " | Sound speed: ", .show_mean(sound_speed),
    " ", sound_speed_units
  )
}

#' Show elastic material properties helper
#' @noRd
.show_elastic_speed_line <- function(prefix,
                                     density,
                                     sound_speed_longitudinal,
                                     sound_speed_transversal,
                                     density_units = "kg m^-3",
                                     speed_units = "m s^-1") {
  # Format one elastic-material summary line ===================================
  paste0(
    prefix,
    " Density: ", .show_mean(density),
    " ", density_units,
    " | cL: ", .show_mean(sound_speed_longitudinal),
    " ", speed_units,
    " | cT: ", .show_mean(sound_speed_transversal),
    " ", speed_units
  )
}

.show_shape_section <- function(section_name,
                                shape_meta,
                                units,
                                radius_values = NULL,
                                segment_label = "segments") {
  # Summarize the available radius statistics ==================================
  mean_radius <- if (!is.null(radius_values)) {
    mean(radius_values, na.rm = TRUE)
  } else {
    NULL
  }
  max_radius <- if (!is.null(radius_values)) {
    max(radius_values, na.rm = TRUE)
  } else {
    NULL
  }
  # Delegate the shared dimension formatting ===================================
  .show_dimension_lines(
    section_name = section_name,
    length_value = shape_meta$length,
    units = units,
    n_segments = shape_meta$n_segments,
    mean_radius = mean_radius,
    max_radius = max_radius,
    segment_label = segment_label
  )
}

#' Show component radius helper
#' @noRd
.show_component_radius_values <- function(component,
                                          fallback_context = "component") {
  # Prefer explicitly stored radius values =====================================
  if (!is.null(component$radius) && !all(is.na(component$radius))) {
    return(component$radius)
  }
  # Otherwise reconstruct the radius profile from geometry =====================
  .shape_radius_profile(
    body = component,
    row_major = TRUE,
    error_context = fallback_context
  )
}

.show_fluid_material_section <- function(section_name,
                                         density = NULL,
                                         sound_speed = NULL,
                                         g = NULL,
                                         h = NULL) {
  # Prefer absolute density and sound-speed values when present ================
  if ((!is.null(density) && !all(is.na(density))) ||
      (!is.null(sound_speed) && !all(is.na(sound_speed)))) {
    return(.show_section_lines(
      section_name,
      .show_density_speed_line(
        prefix = "",
        density = density,
        sound_speed = sound_speed
      )
    ))
  }
  # Otherwise print the stored contrast values =================================
  .show_section_lines(
    section_name,
    paste0(" g: ", .show_mean(g)),
    paste0(" h: ", .show_mean(h))
  )
}

#' Show summary elastic material properties helper
#' @noRd
.show_elastic_material_section <- function(section_name,
                                           density,
                                           sound_speed_longitudinal,
                                           sound_speed_transversal) {
  # Format one elastic-material section ========================================
  .show_section_lines(
    section_name,
    .show_elastic_speed_line(
      prefix = "",
      density = density,
      sound_speed_longitudinal = sound_speed_longitudinal,
      sound_speed_transversal = sound_speed_transversal
    )
  )
}
################################################################################
#' Generic function for show(...) for different scatterers.
#' @param object Scattering object.
#' @importFrom methods setMethod show
#' @keywords internal
#' @export
setMethod(
  f = "show",
  signature = "Scatterer",
  definition = function(object) {
    # Detect scatterer type ====================================================
    sc_type <- class(object)
    # Toggle through scatterer types ===========================================
    switch(sc_type,
      FLS = fls_show(object),
      SBF = sbf_show(object),
      BBF = bbf_show(object),
      CAL = cal_show(object),
      GAS = gas_show(object),
      ESS = ess_show(object)
    )
  }
)
################################################################################
#' show(...) for FLS-class objects.
#' @param object FLS-class object.
#' @keywords internal
#' @export
fls_show <- function(object) {
  meta <- acousticTS::extract(object, "metadata")
  shape <- acousticTS::extract(object, "shape_parameters")
  body <- acousticTS::extract(object, "body")
  radius_values <- .show_component_radius_values(body, "FLS body")

  .show_emit(.show_summary_lines(
    object, "Fluid-like scatterer ", meta,
    .show_shape_section(
      section_name = "Body dimensions",
      shape_meta = shape,
      units = shape$length_units,
      radius_values = radius_values,
      segment_label = "cylinders"
    ),
    "Shape parameters:",
    paste0(" Defined shape:", shape$shape),
    paste0(" L/a ratio:", round(shape$length / shape$radius, 1)),
    paste0(" Taper order:", shape$taper_order %||% "N/A"),
    .show_fluid_material_section(
      section_name = "Material properties",
      density = body$density,
      sound_speed = body$sound_speed,
      g = body$g,
      h = body$h
    ),
    .show_orientation_line(
      "Body orientation (relative to transducer face/axis):",
      body$theta,
      shape$theta_units
    )
  ))
}
#' show(...) for GAS_class objects
#' @param object GAS-class object
#' @keywords internal
#' @export
gas_show <- function(object) {
  meta <- acousticTS::extract(object, "metadata")
  shape <- acousticTS::extract(object, "shape_parameters")
  body <- acousticTS::extract(object, "body")

  .show_emit(.show_summary_lines(
    object, "Gas- and fluid-filled scatterer ", meta,
    .show_section_lines(
      "Body dimensions",
      paste0(" Diameter:", shape$radius * 2, " ", shape$radius_units),
      paste0(" Radius:", shape$radius, " ", shape$radius_units)
    ),
    .show_fluid_material_section(
      section_name = "Material properties",
      density = body$density,
      sound_speed = body$sound_speed,
      g = body$g,
      h = body$h
    )
  ))
}
#' show(...) for SBF-class objects.
#' @param object SBF_class object.
#' @keywords internal
#' @export
sbf_show <- function(object) {
  meta <- acousticTS::extract(object, "metadata")
  shape <- acousticTS::extract(object, "shape_parameters")
  body <- acousticTS::extract(object, "body")
  bladder <- acousticTS::extract(object, "bladder")

  .show_emit(.show_summary_lines(
    object, "Swimbladdered fish (SBF) ", meta,
    .show_shape_section(
      section_name = "Body dimensions",
      shape_meta = shape$body,
      units = shape$length_units,
      radius_values = .show_component_radius_values(body, "SBF body"),
      segment_label = "cylinders"
    ),
    .show_shape_section(
      section_name = "Bladder dimensions",
      shape_meta = shape$bladder,
      units = shape$length_units,
      radius_values = .show_component_radius_values(bladder, "SBF bladder"),
      segment_label = "cylinders"
    ),
    .show_fluid_material_section(
      section_name = "Body material properties",
      density = body$density,
      sound_speed = body$sound_speed,
      g = body$g,
      h = body$h
    ),
    .show_fluid_material_section(
      section_name = "Bladder fluid material properties",
      density = bladder$density,
      sound_speed = bladder$sound_speed,
      g = bladder$g,
      h = bladder$h
    ),
    .show_orientation_line(
      "Body orientation (relative to transducer face/axis):",
      body$theta,
      shape$theta_units
    )
  ))
}
#' show(...) for BBF-class objects.
#' @param object BBF-class object.
#' @keywords internal
#' @export
bbf_show <- function(object) {
  meta <- acousticTS::extract(object, "metadata")
  shape <- acousticTS::extract(object, "shape_parameters")
  body <- acousticTS::extract(object, "body")
  backbone <- acousticTS::extract(object, "backbone")

  .show_emit(.show_summary_lines(
    object, "Backboned fish (BBF) ", meta,
    .show_shape_section(
      section_name = "Body dimensions",
      shape_meta = shape$body,
      units = shape$length_units,
      radius_values = .show_component_radius_values(body, "BBF body"),
      segment_label = "segments"
    ),
    .show_shape_section(
      section_name = "Backbone dimensions",
      shape_meta = shape$backbone,
      units = shape$length_units,
      radius_values = .show_component_radius_values(backbone, "BBF backbone"),
      segment_label = "segments"
    ),
    .show_fluid_material_section(
      section_name = "Body material properties",
      density = body$density,
      sound_speed = body$sound_speed,
      g = body$g,
      h = body$h
    ),
    .show_elastic_material_section(
      section_name = "Backbone elastic properties",
      density = backbone$density,
      sound_speed_longitudinal = backbone$sound_speed_longitudinal,
      sound_speed_transversal = backbone$sound_speed_transversal
    ),
    paste0(
      "Body orientation:", round(body$theta, 3), " ", shape$theta_units,
      " | Backbone orientation:",
      round(backbone$theta, 3), " ", shape$theta_units
    )
  ))
}
#' show(...) for CAL-class objects.
#' @param object CAL-class object.
#' @keywords internal
#' @export
cal_show <- function(object) {
  meta <- acousticTS::extract(object, "metadata")
  shape <- acousticTS::extract(object, "shape_parameters")
  body <- acousticTS::extract(object, "body")

  .show_emit(.show_summary_lines(
    object, "Calibration sphere", meta,
    paste0("Material:", meta$Material),
    paste0(" Sphere longitudinal sound speed:",
           body$sound_speed_longitudinal,
           "m/s"),
    paste0(" Sphere transversal sound speed:",
           body$sound_speed_transversal,
           "m/s"),
    paste0(" Sphere density:", body$density, "kg/m^3"),
    paste0("Diameter:", shape$diameter, " ", shape$diameter_units),
    paste0(" Radius:", shape$radius, " ", shape$diameter_units),
    .show_orientation_line(
      "Propagation direction of the incident sound wave:",
      body$theta,
      shape$theta_units
    )
  ))
}
#' show(...) for ESS-class objects.
#' @param object ESS-class object.
#' @keywords internal
#' @export
ess_show <- function(object) {
  # Parse metadata =============================================================
  meta <- acousticTS::extract(
    object,
    "metadata"
  )
  # Parse shape ================================================================
  shape <- acousticTS::extract(
    object,
    "shape_parameters"
  )
  # Parse shell ================================================================
  shell <- acousticTS::extract(
    object,
    "shell"
  )
  # Parse fluid ================================================================
  fluid <- acousticTS::extract(
    object,
    "fluid"
  )
  # Create the shell material summary ==========================================
  shell_material_props <- shell[
    names(shell) %in% c(
      "sound_speed", "density", "g", "h", "K",
      "E", "G", "nu"
    )
  ]
  shell_material_text <- .show_property_block(
    shell_material_props,
    label_map = c(
      density = "Density",
      sound_speed = "Sound speed",
      K = "Bulk modulus (K)",
      E = "Young's modulus (E)",
      G = "Shear modulus (G)",
      nu = "Poisson's ratio"
    ),
    unit_map = c(
      density = " kg m^-3",
      sound_speed = " m s^-1",
      K = " Pa",
      E = " Pa",
      G = " Pa"
    )
  )
  # Create the internal-fluid material summary =================================
  fluid_material_props <- fluid[
    names(fluid) %in% c("sound_speed", "density", "g", "h")
  ]
  fluid_material_text <- .show_property_block(
    fluid_material_props,
    label_map = c(
      density = "Density",
      sound_speed = "Sound speed"
    ),
    unit_map = c(
      density = " kg m^-3",
      sound_speed = " m s^-1"
    )
  )
  # Print object summary information ===========================================
  cat(
    paste0(methods::is(object)[[1]], "-object"), "\n",
    "Elastic-shelled scatterer", "\n",
    " ID:",
    paste0(meta$ID), "\n",
    "Material:",
    paste0(meta$Material), "\n",
    "  Shell: \n",
    shell_material_text, " \n",
    "  Internal fluid-like body: \n",
    fluid_material_text, " \n",
    "Shape: \n",
    "  Shell: \n",
    "    Radius:", paste0(
      shape$shell$radius,
      " ",
      shape$shell$length_units
    ), " \n",
    "    Diameter:", paste0(
      shape$shell$diameter,
      " ",
      shape$shell$length_units
    ), " \n",
    "    Outer thickness:", paste0(
      shell$shell_thickness,
      " ",
      shape$shell$length_units
    ), "\n",
    "  Internal fluid-like body: \n",
    "    Radius:", paste0(
      shape$fluid$radius,
      " ",
      shape$fluid$length_units
    ), " \n",
    "    Diameter:", paste0(
      shape$fluid$diameter,
      " ",
      shape$fluid$length_units
    ), " \n",
    "Propagation direction of the incident sound wave:",
    paste0(round(shell$theta, 3), " radians")
  )
}
