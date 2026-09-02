################################################################################
# Show functions
################################################################################
################################################################################
# Methods for "show(...)" for each scattering class object
################################################################################
#' Format a numeric mean for display
#' @noRd
.show_mean <- function(x) {
  # Return a printable NA token for missing inputs =============================
  if (is.null(x) || all(is.na(x))) {
    return("NA")
  }
  # Otherwise return the rounded mean ==========================================
  round(mean(x, na.rm = TRUE), 4)
}
#' Format named properties for display
#' @noRd
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

#' Assemble the common scatterer summary header
#' @noRd
.show_summary_lines <- function(object, descriptor, meta, ...) {
  # Combine the header and body sections =======================================
  c(
    paste0(methods::is(object)[[1]], "-object"),
    paste0(" ", descriptor),
    paste0(" ID:", meta$ID),
    ...
  )
}

#' Format a shape section for display
#' @noRd
.show_shape_section <- function(section_name,
                                shape_meta,
                                units,
                                radius_values = NULL,
                                segment_label = "segments") {
  # Build the primary length line ==============================================
  length_line <- paste0(
    " Length:",
    round(shape_meta$length, 3),
    " ",
    units,
    if (!is.null(shape_meta$n_segments)) {
      paste0("(n = ", shape_meta$n_segments, " ", segment_label, ")")
    } else {
      ""
    }
  )
  # Append the optional radius summaries =======================================
  lines <- c(paste0(section_name, ":"), length_line)
  if (!is.null(radius_values)) {
    lines <- c(
      lines,
      paste0(" Mean radius:", round(mean(radius_values, na.rm = TRUE), 4),
        " ", units),
      paste0(" Max radius:", round(max(radius_values, na.rm = TRUE), 4),
        " ", units)
    )
  }

  lines
}

#' Resolve component radius values for display
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

#' Format fluid material properties for display
#' @noRd
.show_fluid_material_section <- function(section_name,
                                         density = NULL,
                                         sound_speed = NULL,
                                         g = NULL,
                                         h = NULL) {
  # Prefer absolute density and sound-speed values when present ================
  if ((!is.null(density) && !all(is.na(density))) ||
    (!is.null(sound_speed) && !all(is.na(sound_speed)))) {
    return(c(
      paste0(section_name, ":"),
      paste0(
        " Density: ", .show_mean(density),
        " kg m^-3 | Sound speed: ", .show_mean(sound_speed),
        " m s^-1"
      )
    ))
  }
  # Otherwise print the stored contrast values =================================
  c(
    paste0(section_name, ":"),
    paste0(" g: ", .show_mean(g)),
    paste0(" h: ", .show_mean(h))
  )
}

#' Format elastic material properties for display
#' @noRd
.show_elastic_material_section <- function(section_name,
                                           density,
                                           sound_speed_longitudinal,
                                           sound_speed_transversal) {
  # Format one elastic-material section ========================================
  c(
    paste0(section_name, ":"),
    paste0(
      " Density: ", .show_mean(density),
      " kg m^-3 | cL: ", .show_mean(sound_speed_longitudinal),
      " m s^-1 | cT: ", .show_mean(sound_speed_transversal),
      " m s^-1"
    )
  )
}
################################################################################
#' Display a scatterer object
#' @param object Scattering object.
#' @return Called for its side effect of printing a formatted summary;
#'   invisibly returns \code{NULL}.
#' @examples
#' show(cal_generate(material = "WC", n_segments = 20))
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
#' Display an FLS object
#' @noRd
fls_show <- function(object) {
  meta <- acousticTS::extract(object, "metadata")
  shape <- acousticTS::extract(object, "shape_parameters")
  body <- acousticTS::extract(object, "body")
  radius_values <- .show_component_radius_values(body, "FLS body")

  cat(paste(.show_summary_lines(
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
    paste0(
      "Body orientation (relative to transducer face/axis):",
      round(body$theta, 3),
      " ",
      shape$theta_units
    )
  ), collapse = "\n"))
}
#' Display a GAS object
#' @noRd
gas_show <- function(object) {
  meta <- acousticTS::extract(object, "metadata")
  shape <- acousticTS::extract(object, "shape_parameters")
  body <- acousticTS::extract(object, "body")

  cat(paste(.show_summary_lines(
    object, "Gas- and fluid-filled scatterer ", meta,
    c(
      "Body dimensions:",
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
  ), collapse = "\n"))
}
#' Display an SBF object
#' @noRd
sbf_show <- function(object) {
  meta <- acousticTS::extract(object, "metadata")
  shape <- acousticTS::extract(object, "shape_parameters")
  body <- acousticTS::extract(object, "body")
  bladder <- acousticTS::extract(object, "bladder")

  cat(paste(.show_summary_lines(
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
    paste0(
      "Body orientation (relative to transducer face/axis):",
      round(body$theta, 3),
      " ",
      shape$theta_units
    )
  ), collapse = "\n"))
}
#' Display a BBF object
#' @noRd
bbf_show <- function(object) {
  meta <- acousticTS::extract(object, "metadata")
  shape <- acousticTS::extract(object, "shape_parameters")
  body <- acousticTS::extract(object, "body")
  backbone <- acousticTS::extract(object, "backbone")

  cat(paste(.show_summary_lines(
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
  ), collapse = "\n"))
}
#' Display a CAL object
#' @noRd
cal_show <- function(object) {
  meta <- acousticTS::extract(object, "metadata")
  shape <- acousticTS::extract(object, "shape_parameters")
  body <- acousticTS::extract(object, "body")

  cat(paste(.show_summary_lines(
    object, "Calibration sphere", meta,
    paste0("Material:", meta$Material),
    paste0(
      " Sphere longitudinal sound speed:",
      body$sound_speed_longitudinal,
      "m/s"
    ),
    paste0(
      " Sphere transversal sound speed:",
      body$sound_speed_transversal,
      "m/s"
    ),
    paste0(" Sphere density:", body$density, "kg/m^3"),
    paste0("Diameter:", shape$diameter, " ", shape$diameter_units),
    paste0(" Radius:", shape$radius, " ", shape$diameter_units),
    paste0(
      "Propagation direction of the incident sound wave:",
      round(body$theta, 3),
      " ",
      shape$theta_units
    )
  ), collapse = "\n"))
}
#' Display an ESS object
#' @noRd
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
