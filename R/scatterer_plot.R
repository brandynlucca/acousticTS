################################################################################
################################################################################
# GENERIC PLOT FUNCTION FOR BODY SHAPES AND MODELS
################################################################################
################################################################################
#' Method for what is printed for objects.
#' @param x Scatterer-class object.
#' @param y Ignored (required for plot method signature).
#' @param ... Additional arguments passed to plotting functions
#' @param type Toggle between body shape ("shape") or modeling results ("model")
#' @param x_units If "model" is selected, then toggle between frequency
#'    ("frequency", kHz) or ka ("ka").
#' @param nudge_y y-axis nudge.
#' @param nudge_x x-axis nudge.
#' @param aspect_ratio Aspect ratio setting ( defaults to "manual" for nudge_y
#' and nudge_x to apply; otherwise, input "auto").
#' @param y_units y-axis data selection (e.g. TS, sigma_bs -- defaults to TS)
#' @param ... Additional plot inputs
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
#' @param type Toggle between body shape ("shape") or modeling results ("model")
#' @param x_units If "model" is selected, then toggle between frequency
#'    ("frequency", kHz) or ka ("ka").
#' @param nudge_y y-axis nudge.
#' @param nudge_x x-axis nudge.
#' @param ... Additional plot inputs
#' @export
cal_plot <- function(object,
                     type = "shape",
                     nudge_y = 1.01,
                     nudge_x = 1.01,
                     x_units = "frequency", ...) {
  # Retrieve default plot window parameters ====================================
  opar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(opar))
  if (type == "shape") {
    body <- acousticTS::extract(
      object,
      "body"
    )
    .plot_axisymmetric_outline(
      position_matrix = body$rpos,
      shape_parameters = list(radius = body$radius),
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
  # Retrieve default plot window parameters ====================================
  opar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(opar))
  if (type == "shape") {
    shell <- extract(object, "shell")
    .plot_axisymmetric_outline(
      position_matrix = shell$rpos,
      shape_parameters = list(radius = shell$radius),
      xlab = "Semi-major diameter (m)",
      ylab = "Semi-minor diameter (m)",
      nudge_y = nudge_y,
      center_x = TRUE,
      mar = c(4, 4.5, 1, 2)
    )
    fluid <- extract(object, "fluid")
    if ("rpos" %in% names(fluid)) {
      .plot_axisymmetric_outline(
        position_matrix = fluid$rpos,
        shape_parameters = list(radius = fluid$radius),
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
#' @param type Toggle between body shape ("shape") or modeling results ("model")
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
  # Retrieve default plot window parameters ====================================
  opar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(opar))
  if (type == "shape") {
    # Extract body shape information ===========================================
    body <- acousticTS::extract(
      object,
      "body"
    )
    # Extract generic shape information ========================================
    shape_params <- acousticTS::extract(object, "shape_parameters")
    # Define plot margins ======================================================
    graphics::par(
      ask = FALSE,
      oma = c(1, 1, 1, 0),
      mar = c(5.0, 4.5, 1.5, 2)
    )
    # Center shape =============================================================
    body$rpos[3, ] <- body$rpos[3, ] - stats::median(body$rpos[3, ])
    # Get radius ===============================================================
    if ("radius_shape" %in% names(shape_params)) {
      radius <- shape_params$radius_shape
    } else if (!is.null(body$radius) && !all(is.na(body$radius))) {
      radius <- body$radius
    } else if ("zU" %in% rownames(body$rpos)) {
      # Arbitrary shape with dorsal/ventral rows: use zU as half-thickness
      radius <- body$rpos["zU", ]
    } else {
      radius <- rep(0, ncol(body$rpos))
    }
    # Sort index ===============================================================
    if (body$rpos[1, 1] > body$rpos[1, ncol(body$rpos)]) {
      body$rpos <- body$rpos[, rev(seq_len(ncol(body$rpos))),
        drop = FALSE
      ]
      if (!is.null(radius)) radius <- rev(radius)
      # sort_idx <- order(body$rpos[1, ])
      # body$rpos <- body$rpos[, sort_idx, drop = FALSE]
      # if (!is.null(body$radius)) body$radius <- body$radius[sort_idx]
    }
    # Adjust axes ==============================================================
    if (aspect_ratio == "manual") {
      vert_lims <- c(
        min(body$rpos[3, ] - radius) * (1 - (1 - nudge_y)),
        max(body$rpos[3, ] + radius) * nudge_y
      )
    } else {
      vert_lims <- c(
        -max(body$rpos[1, ]) * 0.10,
        max(body$rpos[1, ]) * 0.10
      )
    }
    # Begin plotting ===========================================================
    # graphics::plot( x = body$rpos[ 1 , ] ,
    #                 y = body$rpos[ 3 , ] ,
    #                 type = 'l' ,
    #                 lwd = 4 ,
    #                 cex.lab = 1.2 ,
    #                 cex.axis = 1.2 ,
    #                 xlab = "Length (m)" ,
    #                 ylab = "Thickness (m)" ,
    #                 ylim = vert_lims )
    plot(body$rpos[1, ], body$rpos[3, ],
      type = "n",
      xlab = "Length (m)", ylab = "Thickness (m)",
      ylim = c(
        min(body$rpos[3, ] - radius) * 1.1,
        max(body$rpos[3, ] + radius) * 1.1
      )
    )
    # Add lower perimeter of shape =============================================
    # graphics::lines( x = body$rpos[ 1 , ] ,
    #                  y = body$rpos[ 3 , ] - body$radius ,
    #                  lty = 1 ,
    #                  lwd = 4 )
    # Add upper perimeter of shape =============================================
    # graphics::lines( x = body$rpos[ 1 , ] ,
    #                  y = body$rpos[ 3 , ] + body$radius ,
    #                  lty = 1 ,
    #                  lwd = 4 )
    # Add body segments ========================================================
    # graphics::segments( x0 = body$rpos[ 1 , ] ,
    #                     x1 = body$rpos[ 1 , ] ,
    #                     y0 = body$rpos[ 3 , ] - body$radius ,
    #                     y1 = body$rpos[ 3 , ] + body$radius ,
    #                     lty = 3 ,
    #                     lwd = 1.25 )
    # Draw angled "cylinders" for each segment =================================
    n_segments <- shape_params$n_segments
    # count <- 0
    # rc <- c()
    # ---- Iterate =============================================================
    for (i in 1:(n_segments)) {
      # ---- Leading coordinates
      x0 <- body$rpos[1, i]
      y0 <- body$rpos[3, i]
      # ---- Trailing coordinates
      x1 <- body$rpos[1, i + 1]
      y1 <- body$rpos[3, i + 1]
      # ---- Thickness
      r <- radius[i]
      # rc <- c(rc, r)
      if (r == 0) next
      # cat(sprintf("Segment %d [%d, %d]: (%.6f, %.6f) -> (%.6f, %.6f),
      # r=%.6f\n",
      #             i, order(body$rpos[1,])[i], order(body$rpos[1,])[i+1],
      # x0, y0, x1, y1, r))
      # ---- Direction vector
      dx <- x1 - x0
      dy <- y1 - y0
      len <- sqrt(dx^2 + dy^2)
      # ---- Perpendicular vector (unit)
      px <- -dy / len
      py <- dx / len
      # ---- Polygon nodes
      xA <- x0 + r * px
      yA <- y0 + r * py
      xB <- x1 + r * px
      yB <- y1 + r * py
      xC <- x1 - r * px
      yC <- y1 - r * py
      xD <- x0 - r * px
      yD <- y0 - r * py
      # ---- Draw polygon
      graphics::polygon(
        x = c(xA, xB, xC, xD),
        y = c(yA, yB, yC, yD),
        col = grDevices::adjustcolor("gray50", alpha.f = 0.6),
        border = "black", lwd = 1
      )
      # count <- count + 1
    }
    # cat("Polygons drawn:", count, "\n")
    # Draw centerline and points
    graphics::lines(body$rpos[1, ], body$rpos[3, ], lwd = 3, col = "gray90")
    graphics::points(body$rpos[1, ], body$rpos[3, ],
      pch = 1, col = "black",
      cex = 0.8
    )
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
  # Retrieve default plot window parameters ====================================
  opar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(opar))
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
  # Retrieve default plot window parameters ===================
  opar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(opar))
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
#' @param type Toggle between body shape ("shape") or modeling results ("model")
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
  opar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(opar))

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
