################################################################################
################################################################################
# GENERIC PLOT FUNCTION FOR BODY SHAPES AND MODELS
################################################################################
################################################################################
#' Plot scatterer geometry, stored model results, or stored TMM scattering views
#'
#' @description
#' S3 plotting method for `Scatterer` objects. Depending on `type`, the method
#' draws the target geometry, one or more stored model curves, or a stored
#' single-target `TMM` scattering slice / map when the object was computed with
#' retained T-matrix state.
#'
#' @param x Scatterer-class object.
#' @param y Ignored (required for the base `plot()` method signature).
#' @param type Plot mode. Use `"shape"` for the stored geometry, `"model"` for
#'   the currently stored model output, or `"scattering"` for stored `TMM`
#'   angular products.
#' @param nudge_y Multiplicative vertical padding applied to automatically
#'   derived y-axis limits.
#' @param nudge_x Multiplicative horizontal padding applied to automatically
#'   derived x-axis limits.
#' @param aspect_ratio Aspect-ratio control for supported shape plots. Use
#'   `"manual"` to honor `nudge_x` / `nudge_y`, or `"auto"` to let the plotting
#'   helper derive the aspect ratio directly.
#' @param x_units Horizontal-axis convention for `type = "model"`. Supported
#'   values depend on the stored model family and typically include
#'   `"frequency"` and geometry-scaled wavenumber variants such as `"ka"`.
#' @param y_units Stored response quantity for `type = "model"` or
#'   `type = "scattering"`. Typical values include `"TS"` and `"sigma_bs"`.
#' @param ... Additional arguments passed through to the class-specific plotting
#'   helpers. For stored `TMM` scattering plots, this includes controls such as
#'   `frequency`, `vary`, `polar`, and `heatmap`.
#'
#' @return Called for its side effect of drawing a plot; returns the input
#'   invisibly.
#'
#' @details
#' This method dispatches to the relevant scatterer-class plotting helper:
#' `cal_plot()`, `ess_plot()`, `sbf_plot()`, `bbf_plot()`, `fls_plot()`, or
#' `gas_plot()`. The supported `type` values therefore depend on what has been
#' stored on the object. For example, `type = "model"` requires the object to
#' already contain model output, and `type = "scattering"` currently applies
#' only to stored `TMM` results.
#'
#' @seealso [extract()], [target_strength()], [tmm_scattering()],
#'   [tmm_scattering_grid()]
#'
#' @rdname plot.Scatterer
#' @aliases plot.Scatterer
#'
#' @keywords utility, plotting
#' @export
plot.Scatterer <- function(x,
                           y = NULL,
                           type = "shape",
                           nudge_y = 1.1,
                           nudge_x = 1.05,
                           aspect_ratio = "manual",
                           x_units = "frequency",
                           y_units = "TS", ...) {
  # Detect scatterer type ======================================================
  sc_type <- class(x)
  # Switch to sub-class-specific plotting method ===============================
  switch(sc_type,
         CAL = cal_plot(
           object = x, type = type, nudge_y = nudge_y, nudge_x = nudge_x,
           x_units = x_units, ...
         ),
         ESS = ess_plot(
           object = x, type = type, nudge_y = nudge_y, nudge_x = nudge_x,
           x_units = x_units, y_units = y_units, ...
         ),
         SBF = sbf_plot(
           object = x, type = type, nudge_y = nudge_y, nudge_x = nudge_x,
           aspect_ratio = aspect_ratio, x_units = x_units, y_units = y_units, ...
         ),
         BBF = bbf_plot(
           object = x, type = type, nudge_y = nudge_y, nudge_x = nudge_x,
           aspect_ratio = aspect_ratio, x_units = x_units, y_units = y_units, ...
         ),
         FLS = fls_plot(
           object = x, type = type, nudge_y = nudge_y, nudge_x = nudge_x,
           aspect_ratio = aspect_ratio, x_units = x_units, y_units = y_units, ...
         ),
         GAS = gas_plot(
           object = x, type = type, nudge_y = nudge_y, nudge_x = nudge_x,
           x_units = x_units, y_units = y_units, ...
         )
  )
}
################################################################################
#' Base plotting color palette
#' @description
#' Color palette vector referenced for plotting multiple models simultaneously.
#' @rdname model_palette
#' @noRd
model_palette <- c(
  "black", "azure4", "royalblue4", "firebrick3", "steelblue3",
  "orangered3", "darkolivegreen", "sienna3", "cadetblue"
)
################################################################################
################################################################################
# ELASTIC SHELLED SCATTERERS
################################################################################
################################################################################
#' Plotting for CAL-class objects
#' @param object CAL-class object.
#' @param type Toggle between body shape (`"shape"`), modeling results
#'   (`"model"`), or stored TMM angular slices (`"scattering"`).
#' @param x_units If "model" is selected, then toggle between frequency
#'    ("frequency", kHz) or ka ("ka").
#' @param nudge_y y-axis nudge.
#' @param nudge_x x-axis nudge.
#' @param ... Additional plot inputs
#' @keywords internal
#' @noRd
cal_plot <- function(object,
                     type = "shape",
                     nudge_y = 1.01,
                     nudge_x = 1.01,
                     x_units = "frequency", ...) {
  # Respect caller-managed panel layouts by avoiding a full par() reset =======
  if (type == "shape") {
    body <- acousticTS::extract(
      object,
      "body"
    )
    .plot_axisymmetric_outline(
      position_matrix = body$rpos,
      shape_parameters = acousticTS::extract(object, "shape_parameters"),
      xlab = "Semi-major diameter (m)",
      ylab = "Semi-minor diameter (m)",
      lwd = 4,
      nudge_y = nudge_y
    )
  } else if (type == "model") {
    if (length(extract(object, "model")) == 0) {
      stop("ERROR: no model results detected in object.")
    } else {
      # Extract body shape information ============================
      shape <- extract(object, "body")
      # Extract model results ====================================
      TS <- extract(object, "model")$calibration$TS
      x_axis_domain <- extract(
        object,
        "model_parameters"
      )$calibration$parameters
      if (x_units == "frequency") {
        x_axis <- x_axis_domain$acoustics$frequency * 1e-3
        x_lab <- "Frequency (kHz)"
      } else if (x_units == "k_sw") {
        x_axis <- x_axis_domain$acoustics$k_sw * shape$radius
        x_lab <- expression(italic(k[sw] * a))
      } else if (x_units == "k_l") {
        x_axis <- x_axis_domain$acoustics$k_l * shape$radius
        x_lab <- expression(italic(k[l] * a))
      } else {
        x_axis <- x_axis_domain$acoustics$k_t * shape$radius
        x_lab <- expression(italic(k[t] * a))
      }
      # Plot results ===============================================
      graphics::par(
        ask = FALSE,
        mar = c(4, 4.5, 1, 1)
      )
      plot(
        x = x_axis,
        y = TS,
        type = "l",
        xlab = x_lab,
        ylab = expression(Target ~ strength ~ (dB ~ re. ~ 1 ~ m^2)),
        lwd = 2.5,
        cex.lab = 1.2,
        cex.axis = 1.2,
        xlim = c(
          abs(min(x_axis) * (1 - (1 + nudge_x))),
          max(x_axis) * (nudge_x)
        ),
        ylim = c(
          min(TS) * (1 - (1 - nudge_y)),
          max(TS) * (1 + (1 - nudge_y))
        ),
        xaxs = "i",
        yaxs = "i"
      )
    }
  }
  invisible()
}
#' Plotting for ESS-class objects
#' @param object ESS-class object.
#' @param type Toggle between body shape ("shape") or modeling results ("model")
#' @param x_units If "model" is selected, then toggle between frequency
#'    ("frequency", kHz) or ka ("ka").
#' @param nudge_y y-axis nudge.
#' @param nudge_x x-axis nudge.
#' @param y_units y-axis data selection (e.g. TS, sigma_bs -- defaults to TS).
#' @param ... Additional plot inputs
#' @noRd
ess_plot <- function(object,
                     type = "shape",
                     nudge_y = 1.02,
                     nudge_x = 1.01,
                     x_units = "frequency",
                     y_units = "TS", ...) {
  # Respect caller-managed panel layouts by avoiding a full par() reset =======
  if (type == "shape") {
    shell <- extract(object, "shell")
    shell_shape <- acousticTS::extract(object, "shape_parameters")$shell
    .plot_axisymmetric_outline(
      position_matrix = shell$rpos,
      shape_parameters = shell_shape,
      xlab = "Semi-major diameter (m)",
      ylab = "Semi-minor diameter (m)",
      nudge_y = nudge_y,
      center_x = TRUE,
      mar = c(4, 4.5, 1, 2)
    )
    fluid <- extract(object, "fluid")
    if ("rpos" %in% names(fluid)) {
      fluid_shape <- acousticTS::extract(object, "shape_parameters")$fluid
      .plot_axisymmetric_outline(
        position_matrix = fluid$rpos,
        shape_parameters = fluid_shape,
        col = "red",
        segment_col = "red",
        center_x = TRUE,
        init = FALSE
      )
      graphics::legend(
        x = "bottomright",
        legend = c("Shell", "Fluid-like body"),
        lty = c(1, 1),
        lwd = c(4, 3.5),
        col = c("black", "red"),
        cex = 0.95
      )
    }
  } else if (type == "model") {
    .plot_model_results(
      .collect_model_plot_data(object, x_units = x_units, y_units = y_units),
      nudge_y = nudge_y,
      nudge_x = nudge_x
    )
  }
  invisible()
}
################################################################################
################################################################################
# FLUID-LIKE SCATTERERS
################################################################################
################################################################################
#' Plotting for FLS-class objects
#' @param object FLS-class object.
#' @param type Toggle between body shape (`"shape"`), modeling results
#'   (`"model"`), or stored TMM angular slices (`"scattering"`).
#' @param x_units If "model" is selected, then toggle between frequency
#'    ("frequency", kHz) or ka ("ka").
#' @param nudge_y y-axis nudge.
#' @param nudge_x x-axis nudge.
#' @param aspect_ratio Aspect ratio setting ( defaults to "manual" for nudge_y
#' and nudge_x to apply; otherwise, input "auto").
#' @param y_units y-axis data selection (e.g. TS, sigma_bs -- defaults to TS)
#' @param ... Additional plot inputs
#' @noRd
fls_plot <- function(object,
                     type = "shape",
                     nudge_y = 1.05,
                     nudge_x = 1.01,
                     aspect_ratio = "manual",
                     x_units = "frequency",
                     y_units = "TS", ...) {
  # Respect caller-managed panel layouts by avoiding a full par() reset =======
  if (type == "shape") {
    body <- acousticTS::extract(object, "body")
    shape_params <- acousticTS::extract(object, "shape_parameters")
    .plot_row_major_segmented_body(
      rpos = body$rpos,
      shape_parameters = shape_params,
      nudge_y = nudge_y,
      aspect_ratio = aspect_ratio
    )
  } else if (type == "model") {
    .plot_model_results(
      .collect_model_plot_data(object, x_units = x_units, y_units = y_units),
      nudge_y = nudge_y,
      nudge_x = nudge_x
    )
  } else if (type == "scattering") {
    .plot_tmm_scattering_slice(
      object = object,
      nudge_y = nudge_y,
      nudge_x = nudge_x,
      ...
    )
  }
  invisible()
}
################################################################################
################################################################################
# GAS-BEARING SCATTERERS
################################################################################
################################################################################
#' Plotting for GAS-class objects
#' @param object GAS-class object.
#' @param type Toggle between body shape ("shape") or modeling results ("model")
#' @param x_units If "model" is selected, then toggle between frequency
#'    ("frequency", kHz) or ka ("ka").
#' @param nudge_y y-axis nudge.
#' @param nudge_x x-axis nudge.
#' @param y_units y-axis data selection (e.g. TS, sigma_bs -- defaults to TS).
#' @param ... Additional plot inputs
#' @noRd
gas_plot <- function(object,
                     type = "shape",
                     nudge_y = 1.01,
                     nudge_x = 1.01,
                     x_units = "frequency",
                     y_units = "TS", ...) {
  # Respect caller-managed panel layouts by avoiding a full par() reset =======
  if (type == "shape") {
    body <- extract(object, "body")
    shape <- acousticTS::extract(object, "shape_parameters")
    .plot_axisymmetric_outline(
      position_matrix = body$rpos,
      shape_parameters = shape,
      xlab = "Semi-major axis (m)",
      ylab = "Semi-minor radius (m)",
      nudge_y = nudge_y
    )
  } else if (type == "model") {
    .plot_model_results(
      .collect_model_plot_data(
        object,
        x_units = x_units,
        y_units = y_units
      ),
      nudge_y = nudge_y,
      nudge_x = nudge_x
    )
  } else if (type == "scattering") {
    .plot_tmm_scattering_slice(
      object = object,
      nudge_y = nudge_y,
      nudge_x = nudge_x,
      ...
    )
  }
  invisible()
}
################################################################################
#' Plotting for SBF-class objects
#' @param object SBF-class object.
#' @param type Toggle between body shape ("shape") or modeling results ("model")
#' @param x_units If "model" is selected, then toggle between frequency
#'    ("frequency", kHz) or ka ("ka").
#' @param nudge_y y-axis nudge.
#' @param nudge_x x-axis nudge.
#' @noRd
sbf_plot <- function(object,
                     type = "shape",
                     nudge_y = 1.05,
                     nudge_x = 1.01,
                     aspect_ratio = "manual",
                     x_units = "frequency",
                     y_units = "TS", ...) {
  # Respect caller-managed panel layouts by avoiding a full par() reset =======
  if (type == "shape") {
    .plot_profile_component_panels(
      primary_rpos = extract(object, "body")$rpos,
      secondary_rpos = extract(object, "bladder")$rpos,
      primary_label = "Body",
      secondary_label = "Resonant feature",
      secondary_col = "red",
      nudge_y = nudge_y,
      nudge_x = nudge_x
    )
  } else if (type == "model") {
    .plot_model_results(
      .collect_model_plot_data(
        object,
        x_units = x_units,
        y_units = y_units
      ),
      nudge_y = nudge_y,
      nudge_x = nudge_x
    )
  }
  invisible()
}
################################################################################
#' Plotting for BBF-class objects
#' @param object BBF-class object.
#' @param type Toggle between body shape (`"shape"`), modeling results
#'   (`"model"`), or stored TMM angular slices (`"scattering"`).
#' @param x_units If "model" is selected, then toggle between frequency
#'    ("frequency", kHz) or ka ("ka").
#' @param nudge_y y-axis nudge.
#' @param nudge_x x-axis nudge.
#' @param aspect_ratio Aspect ratio setting.
#' @param y_units y-axis data selection (e.g. TS, sigma_bs -- defaults to TS).
#' @param ... Additional plot inputs
#' @noRd
bbf_plot <- function(object,
                     type = "shape",
                     nudge_y = 1.05,
                     nudge_x = 1.01,
                     aspect_ratio = "manual",
                     x_units = "frequency",
                     y_units = "TS", ...) {
  # Respect caller-managed panel layouts by avoiding a full par() reset =======

  if (type == "shape") {
    body_component <- extract(object, "body")
    backbone_component <- extract(object, "backbone")
    .plot_profile_component_panels(
      primary_rpos = body_component$rpos,
      primary_radius = body_component$radius,
      secondary_rpos = backbone_component$rpos,
      secondary_radius = backbone_component$radius,
      primary_label = "Body",
      secondary_label = "Backbone",
      secondary_col = "dodgerblue3",
      nudge_y = nudge_y,
      nudge_x = nudge_x
    )
  } else if (type == "model") {
    .plot_model_results(
      .collect_model_plot_data(
        object,
        x_units = x_units,
        y_units = y_units
      ),
      nudge_y = nudge_y,
      nudge_x = nudge_x
    )
  }

  invisible()
}
