################################################################################
# METHOD FUNCTIONS
################################################################################
################################################################################
# Methods for "show(...)" for each scattering class object
################################################################################
#' show(...) for FLS-class objects.
#' @param object FLS-class object.
#' @export
fls_show <- function(object) {
  # Print out informational text ===================
  if("SDWBA" %in% names(extract(object, "model_parameters"))) {
    SDWBA_text <-  paste0("\n SDWBA-initialized: L0 = ",
                          extract(object, "model_parameters")$SDWBA$parameters$lparameters$ength_init,
                          " ",
                          extract(object, "shape_parameters")$body$length_units,
                          "; N0 = ",
                          extract(object, "model_parameters")$SDWBA$parameters$parameters$ncyl_init,
                          " cylinders; phi_sd0 = ",
                          round(extract(object, "model_parameters")$SDWBA$parameters$parameters$phase_sd_init, 3),
                          "; iterations = ",
                          extract(object, "model_parameters")$SDWBA$parameters$parameters$n_iterations)
  } else {
    SDWBA_text <- ""
  }
  cat(paste0(is(object)[[1]], " object"), "\n",
      "Fluid-like scatterer \n",
      "ID:",
      extract(object, "metadata")$ID, "\n",
      "Body length:",
      round(extract(object, "shape_parameters")$body$length, 3),
      extract(object, "shape_parameters")$body$length_units,
      "(n =",
      extract(object, "shape_parameters")$body$ncyl,
      " cylinders)", "\n",
      "Maximum radius: ",
      round(max(extract(object, "body")$radius), 3),
      paste0(extract(object, "shape_parameters")$body$length_units, "; "),
      "Mean radius: ",
      round(mean(extract(object, "body")$radius), 3),
      extract(object, "shape_parameters")$body$length_units,
      "\n",
      "Body orientation (relative to transducer axis): ",
      round(extract(object, "body")$theta, 3),
      extract(object, "shape_parameters")$body$theta_units, "\n",
      "Material properties (body): g =",
      paste0(extract(object, "body")$g,";"),
      "h =",
      extract(object, "body")$h, SDWBA_text)
  }
#' show(...) for SBF-class objects.
#' @param object SBF_class object.
#' @export
sbf_show <- function(object) {
  cat(paste0(is(object)[[1]], " object"), "\n",
      "Gas-filled swimbladdered scatterer \n",
      "ID:",
      extract(object, "metadata")$ID, "\n",
      "Body length:",
      round(extract(object, "shape_parameters")$body$length, 3),
      extract(object, "shape_parameters")$length_units,
      "(n =",
      extract(object, "shape_parameters")$body$ncyl,
      "cylinders) \n",
      "Bladder length:",
      round(extract(object, "shape_parameters")$bladder$length, 3),
      extract(object, "shape_parameters")$length_units,
      "(n =",
      extract(object, "shape_parameters")$bladder$ncyl,
      "cylinders) \n",
      "Body orientation (relative to transducer axis):",
      round(extract(object, "body")$theta, 3),
      extract(object, "shape_parameters")$theta_units, "\n",
      "Bladder orientation (relative to transducer axis):",
      round(extract(object, "bladder")$theta, 3),
      extract(object, "shape_parameters")$theta_units, "\n",
      "Material properties (body): density =",
      extract(object, "body")$density, "(kg/m^3);",
      "sound speed =",
      extract(object, "body")$sound_speed, "(m/s) \n",
      "Material properties (bladder): density =",
      extract(object, "bladder")$density, "(kg/m^3);",
      "sound speed =",
      extract(object, "bladder")$sound_speed, "(m/s)")
  }
#' Generic function for show(...) for different scatterers.
#' @param object Scattering object.
#' @export
setMethod("show",
          signature("scatterer"),
          function(object) {
            # Detect scatterer type ============================
            sc_type <- class(object)

            switch(sc_type,
                   FLS = fls_show(object),
                   SBF = sbf_show(object))
          })
################################################################################
# Methods for "plot(...)" for each scattering class object
################################################################################
#' Plotting for SBF-class objects
#' @param object SBF-class object.
#' @param type Toggle between body shape ("shape") or modeling results ("model")
#' @param x_units If "model" is selected, then toggle between frequency
#'    ("frequency", kHz) or ka ("ka").
#' @param nudge_y y-axis nudge.
#' @param nudge_x x-axis nudge.
#' @export
sbf_plot <- function(object,
                     type = "shape",
                     nudge_y = 1.05,
                     nudge_x = 1.01,
                     x_units = "frequency") {
  # Retrieve default plot window parameters ===================
  opar <- par(no.readonly = TRUE)
  on.exit(par(opar))
  if(type == "shape") {
    # Extract body shape information ============================
    body <- extract(object, "body")$rpos
    # Extract bladder information ===============================
    bladder <- extract(object, "bladder")$rpos
    # Set up plot limits for each shape component ===============
    ylow_lim <- min(body[4, ]) * (1 - (1 - nudge_y))
    yhi_lim <- max(body[3, ]) * nudge_y
    xlow_lim <- min(body[1, ]) * (1 - (1 - nudge_x))
    xhi_lim <- max(body[1, ]) * nudge_x
    # Plot dorsoventral view =====================================
    par(ask = F,
        mfrow = c(2, 1),
        mar = c(2, 4.5, 1, 2),
        oma = c(1.5, 0.5, 0.5, 0))
    # Dorsal =====================================================
    plot(x = body[1, ],
         y = body[3, ],
         type = 'l',
         ylim = c(ylow_lim, yhi_lim),
         xlim = c(xlow_lim, xhi_lim),
         xlab = "",
         ylab = "Height (m)",
         lwd = 4,
         xaxt = "n",
         cex.lab = 1.2,
         cex.axis = 1.2)
    # Ventral ====================================================
    lines(x = body[1, ],
          y = body[4, ],
          lty = 1,
          lwd = 4)
    # Connect anterior/posterior gaps =============================
    segments(x0 = c(min(body[1, ]), max(body[1, ])),
             x1 = c(min(body[1, ]), max(body[1, ])),
             y0 = c(body[4, ][body[1, ] == min(body[1, ])],
                    body[4, ][body[1, ] == max(body[1, ])]),
             y1 = c(body[3, ][body[1, ] == min(body[1, ])],
                    body[3, ][body[1, ] == max(body[1, ])]),
             lty = 1,
             lwd = 4,
             col = 'black')
    # Plot body segments  =========================================
    segments(x0 = body[1, ], x1 = body[1, ], y0 = body[4, ], y1 = body[3, ],
             lty = 3,
             lwd = 1,
             col = 'gray30')
    # Ventral bladder ==============================================
    lines(x = bladder[1, ],
          y = bladder[4, ],
          lty = 1,
          lwd = 3.5,
          col = 'red')
    # Dorsal bladder ================================================
    lines(x = bladder[1, ],
          y = bladder[3, ],
          lty = 1,
          lwd = 3.5,
          col = 'red')
    # Connect anterior/posterior gaps =============================
    segments(x0 = c(min(bladder[1, ]), max(bladder[1, ])),
             x1 = c(min(bladder[1, ]), max(bladder[1, ])),
             y0 = c(bladder[4, ][bladder[1, ] == min(bladder[1, ])],
                    bladder[4, ][bladder[1, ] == max(bladder[1, ])]),
             y1 = c(bladder[3, ][bladder[1, ] == min(bladder[1, ])],
                    bladder[3, ][bladder[1, ] == max(bladder[1, ])]),
             lty = 1,
             lwd = 3.5,
             col = 'red')
    # Plot bladder segments ==========================================
    segments(x0 = bladder[1, ], x1 = bladder[1, ],
             y0 = bladder[4, ], y1 = bladder[3, ],
             lty = 3,
             lwd = 1,
             col = 'red')
    # Create legend ==================================================
    legend(x = "bottomright",
           legend = c("Body","Resonant feature"),
           lty = c(1, 1),
           lwd = c(4, 3.5),
           col = c('black', 'red'),
           cex = 0.95)
    # End dorsoventral view =====================================
    left_limit <- -max(body[2, ] / 2) * (1 - (1 - nudge_y))
    right_limit <- max(body[2, ] / 2) * (1 - (1 - nudge_y))
    # Plot ventral view =========================================
    plot(x = body[1, ],
         y = body[2, ] / 2,
         type = 'l',
         ylim = c(left_limit, right_limit),
         ylab = "Width (m)",
         lwd = 4,
         cex.lab = 1.2,
         cex.axis = 1.2)
    # Opposite side ====================================================
    lines(x = body[1, ],
          y = -body[2, ] / 2,
          lty = 1,
          lwd = 4)
    # Close ends ====================================================
    segments(x0 = c(min(body[1, ]), max(body[1, ])),
             x1 = c(min(body[1, ]), max(body[1, ])),
             y0 = c(-body[2, ][body[1, ] == min(body[1, ])] / 2,
                    -body[2, ][body[1, ] == max(body[1, ])] / 2),
             y1 = c(body[2, ][body[1, ] == min(body[1, ])] / 2,
                    body[2, ][body[1, ] == max(body[1, ])] / 2),
             lty = 1,
             lwd = 4,
             col ='black')
    # Plot body segments ==========================================
    segments(x0 = body[1, ], x1 = body[1, ],
             y0 = -body[2, ] / 2, y1 = body[2, ] / 2,
             lty = 3,
             lwd = 1,
             col = 'gray30')
    # Left bladder ==================================================
    lines(x = bladder[1, ],
          y = bladder[2, ] / 2,
          lty = 1,
          lwd = 3.5,
          col = 'red')
    # Right bladder =================================================
    lines(x = bladder[1, ],
          y = -bladder[2, ]/2,
          lty = 1,
          lwd = 3.5,
          col = 'red')
    # Close end of bladders =========================================
    segments(x0 = c(min(bladder[1, ]), max(bladder[1, ])),
             x1 = c(min(bladder[1, ]), max(bladder[1, ])),
             y0 = c(bladder[2, ][bladder[1, ] == min(bladder[1, ])]/ 2,
                    bladder[2, ][bladder[1, ] == max(bladder[1, ])] / 2),
             y1 = c(-bladder[2, ][bladder[1, ] == min(bladder[1, ])] / 2,
                    -bladder[2, ][bladder[1, ] == max(bladder[1, ])] / 2),
             lty = 1,
             lwd = 3.5,
             col = 'red')
    # Plot bladder segments ==========================================
    segments(x0 = bladder[1 ,], x1 = bladder[1, ],
             y0 = -bladder[2, ] / 2, y1 = bladder[2, ] / 2,
             lty = 3,
             lwd = 1,
             col ='red')
    # Re-add x-axis text ==========================================
    mtext('Along-body axis (m)',
          side = 1,
          outer = TRUE,
          line = 0,
          cex = 1.2)
  } else if(type == "model") {
    if(length(extract(object, "model")) == 0) {
      stop("ERROR: no model results detected in object.")
    } else {
      # Extract body shape information ============================
      shape <- extract(object, "body")
      # Extract model results ====================================
      TS <- extract(object, "model")$KRM$TS
      x_axis_domain <- extract(object, "model_parameters")$KRM$parameters
      if(x_units == "frequency") {
        x_axis <- x_axis_domain$acoustics$frequency * 1e-3
        x_lab <- "Frequency (kHz)"
      } else if(x_units == "k_b") {
        x_axis <- x_axis_domain$acoustics$k_sw * max(shape$rpos[1, ])
        x_lab <- expression(italic(k[b]*L))
      }
      # Plot results ===============================================
      par(ask = F,
          mar = c(4, 4.5, 1, 1))
      plot(x = x_axis,
           y = TS,
           type = 'l',
           xlab = x_lab,
           ylab = expression(Target~strength~(dB~re.~1~m^2)),
           lwd = 2.5,
           cex.lab = 1.2,
           cex.axis = 1.2,
           xlim = c(min(x_axis) * (1 - nudge_x),
                    max(x_axis) * (nudge_x)),
           ylim = c(min(TS) * (1 - (1 - nudge_y)),
                    max(TS) * (1 + (1 - nudge_y))),
           xaxs = "i",
           yaxs = "i")
    }
  }
  invisible()
}
################################################################################
# Methods for "plot(...)" for each scattering class object
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
#' @import graphics
#' @export
cal_plot <- function(object,
                     type = "shape",
                     nudge_y = 1.00,
                     nudge_x = 1.00,
                     x_units = "frequency", ...) {
  # Retrieve default plot window parameters ===================
  opar <- par(no.readonly = TRUE)
  on.exit(par(opar))
  if(type == "shape") {
    # Extract body shape information ============================
    shape <- extract(object, "body")
    # Center sphere =============================================
    shape$rpos[1, ] <- shape$rpos[1, ]
    # Plot sphere ===============================================
    par(ask = F,
        oma = c(1, 1, 1, 0),
        mar = c(3, 4.5, 1, 2))
    plot(x = shape$rpos[1, ],
         y = shape$rpos[3, ],
         type = 'l',
         ylim = c(min(shape$rpos[3, ]) * (1 - (1 - nudge_y)),
                  max(-shape$rpos[3, ] * nudge_y)),
         xlab = "Diameter (m)",
         ylab = "Diameter (m)",
         lwd = 4,
         cex.lab = 1.2,
         cex.axis = 1.2)
    lines(x = shape$rpos[1, ],
          y = -rev(shape$rpos[3, ]),
          lty = 1,
          lwd = 4)
    segments(x0 = shape$rpos[1, ],
             x1 = shape$rpos[1, ],
             y0 = shape$rpos[3, ],
             y1 = -rev(shape$rpos[3, ]),
             lty = 3,
             lwd = 1)
  } else if (type == "model") {
    if(length(extract(object, "model")) == 0) {
      stop("ERROR: no model results detected in object.")
    } else {
      # Extract body shape information ============================
      shape <- extract(object, "body")
      # Extract model results ====================================
      TS <- extract(object, "model")$TS
      x_axis_domain <- extract(object, "model_parameters")$calibration$parameters
      if(x_units == "frequency") {
        x_axis <- x_axis_domain$acoustics$frequency * 1e-3
        x_lab <- "Frequency (kHz)"
      } else if (x_units == "k_sw") {
        x_axis <- x_axis_domain$acoustics$k_sw * shape$radius
        x_lab <- expression(italic(k[sw]*a))
      } else if (x_units == "k_l") {
        x_axis <- x_axis_domain$acoustics$k_l * shape$radius
        x_lab <- expression(italic(k[l]*a))
      } else{
        x_axis <- x_axis_domain$acoustics$k_t * shape$radius
        x_lab <- expression(italic(k[t]*a))
      }
      # Plot results ===============================================
      par(ask = F,
          mar = c(4, 4.5, 1, 1))
      plot(x = x_axis,
           y = TS,
           type = 'l',
           xlab = x_lab,
           ylab = expression(Target~strength~(dB~re.~1~m^2)),
           lwd = 2.5,
           cex.lab = 1.2,
           cex.axis = 1.2,
           xlim = c(min(x_axis) * (1 - nudge_x),
                    max(x_axis) * (nudge_x)),
           ylim = c(min(TS) * (1 - (1 - nudge_y)),
                    max(TS) * (1 + (1 - nudge_y))),
           xaxs = "i",
           yaxs = "i")
    }
  }
  invisible()
}
################################################################################
# Methods for "plot(...)" for each scattering class object
################################################################################
#' Plotting for FLS-class objects
#' @param object FLS-class object.
#' @param type Toggle between body shape ("shape") or modeling results ("model")
#' @param x_units If "model" is selected, then toggle between frequency
#'    ("frequency", kHz) or ka ("ka").
#' @param nudge_y y-axis nudge.
#' @param nudge_x x-axis nudge.
#' @param ... Additional plot inputs
#' @import graphics
#' @import stats
#' @export
fls_plot <- function(object,
                     type = "shape",
                     nudge_y = 1.05,
                     nudge_x = 1.01,
                     x_units = "frequency", ...) {
  # Retrieve default plot window parameters ===================
  opar <- par(no.readonly = TRUE)
  on.exit(par(opar))
  if(type == "shape") {
    # Extract body shape information ============================
    shape <- extract(object, "body")
    # Center shape=============================================
    shape$rpos[3, ] <- shape$rpos[3, ] - median(shape$rpos[3, ])
    # Upper part of body
    z_up <- shape$rpos[3, ] + shape$radius
    # Lower part of body
    z_lo <- shape$rpos[3, ] - shape$radius
    # Plot sphere ===============================================
    par(ask = F,
        oma = c(1, 1, 1, 0),
        mar = c(3, 4.5, 1, 2))
    plot(x = shape$rpos[1, ],
         y = shape$rpos[3, ],
         type = 'l',
         ylim = c(min(shape$rpos[3, ] - shape$radius) * (1 - (1 - nudge_y)),
                  max(shape$rpos[3, ] + shape$radius) * nudge_y),
         xlab = "Diameter (m)",
         ylab = "Diameter (m)",
         lwd = 4,
         cex.lab = 1.2,
         cex.axis = 1.2)
    lines(x = shape$rpos[1, ],
          y = shape$rpos[3, ] - shape$radius,
          lty = 1,
          lwd = 4)
    lines(x = shape$rpos[1, ],
          y = shape$rpos[3, ] + shape$radius,
          lty = 1,
          lwd = 4)
    # Close ends ====================================================
    segments(y0 = c(z_lo[shape$rpos[1, ] == min(shape$rpos[1, ], na.rm = T)],
                    z_lo[shape$rpos[1, ] == max(shape$rpos[1, ], na.rm = T)]),
             y1 = c(z_up[shape$rpos[1, ] == min(shape$rpos[1, ], na.rm = T)],
                    z_up[shape$rpos[1, ] == max(shape$rpos[1, ], na.rm = T)]),
             x0 = c(0, max(shape$rpos[1, ], na.rm = T)),
             x1 = c(0, max(shape$rpos[1, ], na.rm = T)),
             lty = 1,
             lwd = 4,
             col = 'black')
    # Segments ====================================================
    segments(y0 = z_lo,
             y1 = z_up,
             x0 = shape$rpos[1, ],
             x1 = shape$rpos[1, ],
             lty = 3,
             lwd = 1,
             col = "gray30")
  } else if (type == "model") {
    if(length(extract(object, "model")) == 0) {
      stop("ERROR: no model results detected in object.")
    } else {
      # Extract body shape information ============================
      shape <- extract(object, "body")
      # Extract model results ====================================
      mods <- extract(object, "model")
      nms <- names(mods)
      TS <- mods[[nms]]$TS
      x_axis_domain <- extract(object, "model_parameters")[[nms]]$parameters
      if(x_units == "frequency") {
        x_axis <- x_axis_domain$acoustics$frequency * 1e-3
        x_lab <- "Frequency (kHz)"
      } else if (x_units == "k_sw") {
        x_axis <- x_axis_domain$acoustics$k_sw * shape$radius
        x_lab <- expression(italic(k[sw]*a))
      } else if (x_units == "k_b") {
        x_axis <- x_axis_domain$acoustics$k_b * shape$radius
        x_lab <- expression(italic(k*a))
      }
      # Plot results ===============================================
      par(ask = F,
          mar = c(4, 4.5, 1, 1))
      plot(x = x_axis,
           y = TS,
           type = 'l',
           xlab = x_lab,
           ylab = expression(Target~strength~(dB~re.~1~m^2)),
           lwd = 2.5,
           cex.lab = 1.2,
           cex.axis = 1.2,
           xlim = c(min(x_axis) * (1 - nudge_x),
                    max(x_axis) * (nudge_x)),
           # ylim = c(min(TS) * (1 - (1 - nudge_y)),
           #          max(TS) * (1 + (1 - nudge_y))),
           xaxs = "i",
           yaxs = "i", ...)
    }
  }
  invisible()
}
################################################################################
# Generic for "plot(...)" for each scattering class object
################################################################################
#' Method for what is printed for objects.
#'
#' @param x Scatterer-class object.
#' @param type Toggle between body shape ("shape") or modeling results ("model")
#' @param x_units If "model" is selected, then toggle between frequency
#'    ("frequency", kHz) or ka ("ka").
#' @param nudge_y y-axis nudge.
#' @param nudge_x x-axis nudge.
#' @param ... Additional plot inputs
#' @export
setMethod("plot",
          signature(x = "scatterer", y = "missing"),
          definition = function(x,
                                type = "shape",
                                nudge_y = 1.1,
                                nudge_x = 1.05,
                                x_units = "frequency", ...) {
            # Detect scatterer type ============================
            sc_type <- class(x)

            switch(sc_type,
                   CAL = cal_plot(x, type, nudge_y, nudge_x, x_units, ...),
                   SBF = sbf_plot(x, type, nudge_y, nudge_x, x_units, ...),
                   FLS = fls_plot(x, type, nudge_y, nudge_x, x_units, ...))
          })

#' Initialize object for modeling using the DCM.
#'
#' @param object FLS-class object.
#' @param frequency Transmit frequency (kHz)
#' @param radius_cylinder Optional input to override current shape radius.
#' @param radius_curvature_ratio Ratio between body length and the radius of
#'    curvature. Defaults to 3.0.
#' @param radius_cylinder_fun Defines which radius value will be used from the
#'    radius vector. Defaults to "max", but also accepts "mean" and "median".
#' @param length Body length (m).
#' @param g Density contrast.
#' @param h Sound speed contrast.
#' @param theta Body orientation relative to the indicent sound wave.
#' @param sound_speed_sw Seawater sound speed
#' (\ifelse{html}{\out{c<sub>body</sub>}}{\eqn{c_{body}}},
#' m \ifelse{html}{\out{s<sup>-1</sup>}}{\eqn{s^{-1}}}).
#' @param density_sw Seawater density
#' (\ifelse{html}{\out{&rho;<sub>body</sub>}}{\eqn{\rho_{body}}},
#' kg \ifelse{html}{\out{m<sup>3</sup>}}{\eqn{m^3}})
#' @param alpha_B Numerical coefficient (\ifelse{html}{\out{&alpha;<sub>B</sub>}}{\eqn{\alpha_{B}}}).
#' @export
dcm_initialize <- function(object,
                           # Required inputs
                           frequency,
                           # Optional inputs
                           # Shape (radius)
                           radius_cylinder = NULL,
                           radius_curvature_ratio = 3.0,
                           radius_cylinder_fun = "max",
                           # Body parameters
                           length = NULL,
                           g = NULL,
                           h = NULL,
                           theta = NULL,
                           # Medium
                           sound_speed_sw = 1500,
                           density_sw = 1026,
                           # Numerical coefficient
                           alpha_B = 0.8) {
  # Extract model body parameters =============================================
  body <- extract(object, "body")
  length <- ifelse(!is.null(length),
                   length,
                   extract(object, "shape_parameters")$body$length)
  # Determine radius to be used for DCM =======================================
  radius_uniform <- ifelse(!is.null(radius_cylinder),
                           radius_cylinder,
                           switch(radius_cylinder_fun,
                                  max = max(body$radius, na.rm = T),
                                  mean = mean(body$radius, na.rm = T),
                                  median = median(body$radius, na.rm = T)))
  # Calculate radius of curvature either based on user input or ratio ========
  radius_curvature <- ifelse(!is.null(radius_curvature),
                             radius_curvature,
                             length * radius_curvature_ratio)
  # Fill in remaining model inputs required by DCM ===========================
  # Orientation ==============================================================
  theta <- ifelse(!is.null(theta), theta, body$theta)
  # Density contrast, g ======================================================
  g <- ifelse(!is.null(g), g, body$g)
  # Sound speed contrast, h ==================================================
  h <- ifelse(!is.null(h), h, body$h)
  # Define model dataframe for body parameterization =========================
  body_params <- data.frame(length = length,
                            radius = radius_uniform,
                            radius_function = ifelse(!is.null(radius_cylinder),
                                                     "none",
                                                     radius_cylinder_fun),
                            radius_curvature = radius_curvature,
                            theta = theta,
                            g = g,
                            h = h,
                            alpha_B = alpha_B)
  # Define model dataframe for medium parameterization =======================
  medium_params <- data.frame(sound_speed = sound_speed_sw,
                              density = density_sw)
  # Define model dataframe for acoustic parameterization =====================
  acoustics <- data.frame(frequency = frequency,
                          k_sw = kcalc(frequency, sound_speed_sw),
                          k_b = kcalc(frequency, sound_speed_sw * h))
  # Tidy up model parameters to insert into object ===========================
  model_params <- list(body = body_params,
                       medium = medium_params,
                       acoustics = acoustics)
  slot(object, "model_parameters")$DCM <- model_params
  return(object)
}

#' Resizing function for targets.
#' @param object Scatterer-class object.
#' @param ... Additional inputs.
#' @rdname reforge
#' @export
setGeneric("reforge", function(object, ...)
  standardGeneric("reforge"))

#' Resizing function for swimbladdered targets
#' @param object SBF-class object.
#' @param full_size new factored resize of body.
#' @param bladder_size new factored resize of bladder.
#' @export
setMethod("reforge",
          signature(object="SBF"),
          function(object,
                   full_size=1.0,
                   bladder_size=1.0) {

            if(full_size != 1.0){
              lscale <- full_size
              body_rpos <- extract(object, "body")$rpos
              bladder_rpos <- extract(object, "bladder")$rpos
              bladder_origin <- (bladder_rpos[1,1] - min(body_rpos[1, ])) /
                (max(body_rpos[1, ]) - min(body_rpos[1, ]))
              bladder_rpos[1, ] <- bladder_rpos[1, ] - min(bladder_rpos[1, ])
              mscale <- cbind(c(1,0,0,0),
                              c(0,1,0,0),
                              c(0,0,1,0),
                              c(0,0,0,1)) * lscale
              body_rpos <- t(t(body_rpos) %*% mscale)
              bladder_rpos <- t(t(bladder_rpos) %*% mscale)
              bladder_rpos[1, ] <- bladder_rpos[1, ] + bladder_origin*max(body_rpos[1, ])
              slot(object, "body")$rpos <- body_rpos
              slot(object, "bladder")$rpos <- bladder_rpos
            }

            if(bladder_size != 1.0){
              sbscale <- 1 - bladder_size
              body_rpos <- extract(object, "body")$rpos
              bladder_rpos <- extract(object, "bladder")$rpos
              bladder_origin <- (bladder_rpos[1,1] - min(body_rpos[1, ])) /
                (max(body_rpos[1, ]) - min(body_rpos[1, ]))

              vert_dist <- (bladder_rpos[3, ] - bladder_rpos[4, ]) * sbscale
              bladder_rpos[3, ] <- bladder_rpos[3, ] - 0.5*vert_dist
              bladder_rpos[4, ] <- bladder_rpos[4, ] + 0.5*vert_dist
              bladder_rpos[2, ] <- bladder_rpos[2, ] * sbscale

              #diagnostics


              # if(sum(body_rpos[4, body_rpos[1, ] >= min(bladder_rpos[1, ]) &
              #                  body_rpos[1, ] <= max(bladder_rpos[1, ])] -
              #        bladder_rpos[4, ] > 0) > 0 |
              #    sum(bladder_rpos[3, ] -
              #        body_rpos[3, body_rpos[1, ] >= min(bladder_rpos[1, ]) &
              #                  body_rpos[1, ] <= min(bladder_rpos[1, ])] > 0) > 0){
              #   stop("Swimbladder boundary exceeds boundary of scatterer body shape.
              #      Consider adjusting the swimbladder rescaling factor to amend.")
              # }else{
              #   return(object)
              # }

              slot(object, "bladder")$rpos <- bladder_rpos
            }

            return(object)
          })
#' Reforge FLS-class object.
#' @param object FLS-class object.
#' @param full_size New body length resize.
#' @export
setMethod("reforge",
          signature(object="FLS"),
          function(object,
                   full_size=1.0) {
            lscale <- full_size / extract(object, "shape_parameters")$body$length
            body_rpos <- extract(object, "body")$rpos
            radius <- extract(object, "body")$radius
            mscale <- cbind(c(1,0,0),
                            c(0,1,0),
                            c(0,0,1)) * lscale
            body_rpos <- t(t(body_rpos) %*% mscale)
            radius <- radius * lscale
            slot(object, "body")$rpos <- body_rpos
            slot(object, "body")$radius <- radius
            return(object)
          })

#' Initialize SBF-class object for KRM calculations.
#'
#' @param object SBF-class object
#' @param frequency Frequency (Hz).
#' @param sound_speed_sw Seawater sound speed.
#' @param density_sw Seawater density.
#' @param density_body Optional flesh density input.
#' @param density_swimbladder Optional gas density input.
#' @param sound_speed_body Optional flesh sound speed input.
#' @param sound_speed_swimbladder Optional gas sound speed input.
#' @param theta Optional orientation input (relative to incident sound wave).
#' @param body_resize Factor resize for body.
#' @param swimbladder_resize Factor resize for swimbladder.
#' @export
krm_initialize <- function(object,
                           frequency,
                           sound_speed_sw=1500,
                           density_sw=1026,
                           density_body=NULL,
                           density_swimbladder=NULL,
                           sound_speed_body=NULL,
                           sound_speed_swimbladder=NULL,
                           theta=NULL,
                           body_resize=1.0,
                           swimbladder_resize=1.0){
  medium_params <- data.frame(sound_speed=sound_speed_sw,
                              density=density_sw)

  model_params <- list(acoustics=data.frame(frequency=frequency,
                                            k_sw=kcalc(frequency, sound_speed_sw),
                                            k_b=kcalc(frequency,
                                                      ifelse(is.null(sound_speed_body),
                                                             extract(object, "body")$sound_speed,
                                                             sound_speed_body))),
                       parameters=data.frame(ncyl_b=
                                               extract(object, "shape_parameters")$body$ncyl,
                                             ncyl_sb=
                                               extract(object, "shape_parameters")$bladder$ncyl))

  length_body <- extract(object,
                            "shape_parameters")$body$length * body_resize

  body_params <- data.frame(length=length_body,
                            theta=ifelse(is.null(theta),
                                         extract(object, "body")$theta,
                                         theta),
                            density=ifelse(is.null(density_body),
                                           extract(object, "body")$density,
                                           density_body),
                            sound_speed=ifelse(is.null(sound_speed_body),
                                               extract(object, "body")$sound_speed,
                                               sound_speed_body))

  swimbladder_params <- data.frame(inflation_factor=swimbladder_resize,
                                   theta=ifelse(is.null(theta),
                                                extract(object, "bladder")$theta,
                                                theta),
                                   density=ifelse(is.null(density_swimbladder),
                                                  extract(object, "bladder")$density,
                                                  density_swimbladder),
                                   sound_speed=ifelse(is.null(sound_speed_swimbladder),
                                                      extract(object, "bladder")$sound_speed,
                                                      sound_speed_swimbladder))

  slot(object, "model_parameters")$KRM <- list(parameters=model_params,
                                                   medium=medium_params,
                                                   scatterer=list(body=body_params,
                                                                  bladder=swimbladder_params))
  slot(object, "model")$KRM$sigma_bs <- data.frame(frequency=frequency,
                                                       # ka=rep(NA, length(frequency)),
                                                       sigma_bs=rep(NA, length(frequency)))

  object <- reforge(object,
                        full_size=body_resize,
                        bladder_size=swimbladder_resize)

  return(object)
}

#' Initialize CAL-class object for modeling.
#' @param object CAL-class object.
#' @inheritParams krm_initialize
#' @export
calibration_initialize <- function(object,
                                   frequency,
                                   sound_speed_sw=1500,
                                   density_sw=1026) {
  medium_params <- data.frame(sound_speed=sound_speed_sw,
                              density=density_sw)

  model_params <- list(acoustics = data.frame(frequency = frequency,
                                              k_sw = kcalc(frequency, sound_speed_sw),
                                              k_l = kcalc(frequency,
                                                          extract(object, "body")$sound_speed_longitudinal),
                                              k_t = kcalc(frequency,
                                                          extract(object, "body")$sound_speed_transversal)),
                       parameters = data.frame(ncyl_b = extract(object, "shape_parameters")$body$ncyl))

  body_params <- data.frame(diameter = extract(object, "body")$diameter,
                            radius = extract(object, "body")$radius,
                            sound_speed_longitudinal = extract(object, "body")$sound_speed_longitudinal,
                            sound_speed_transversal = extract(object, "body")$sound_speed_transversal,
                            density = extract(object, "body")$density)

  slot(object, "model_parameters")$calibration  <- list(parameters = model_params,
                                                            medium = medium_params,
                                                            scatterer = list(body = body_params))
  slot(object, "model")$calibration$sigma_bs <- data.frame(frequency=frequency,
                                                               sigma_bs=rep(NA, length(frequency)))

  return(object)
}



#' Model backscatter from a given target.
#' @param object Scatterer-class object.
#' @param frequency Frequency (Hz).
#' @param model Model name.
#' @param ... Additional optional model inputs/parameters.
#' @export
target_strength <- function(object,
                            frequency,
                            model,
                            ...) {
  object <- switch(model,
                   anderson = anderson_initialize(object,
                                                  frequency, ...),
                   calibration = calibration_initialize(object,
                                                        frequency, ...),
                   DCM = dcm_initialize(object,
                                        frequency, ...),
                   DWBA = dwba_initialize(object,
                                          frequency, ...),
                   KRM = krm_initialize(object,
                                        frequency, ...),
                   stanton_high_pass = stanton_high_pass_initialize(object,
                                                                    frequency, ...))

  object <- switch(model,
                   anderson = anderson_model(object),
                   calibration = calibration(object),
                   DCM = DCM(object),
                   DWBA = DWBA(object),
                   KRM = KRM(object),
                   stanton_high_pass = stanton_high_pass(object))
  return(object)
}

#' Initialize FLS-class object for TS modeling.
#' @param object FLS-class object.
#' @inheritParams krm_initialize
#' @export
dwba_initialize <- function(object,
                            frequency,
                            sound_speed_sw = 1500,
                            density_sw = 1026,
                            theta = pi / 2) {
  # Define medium parameters =========================================
  medium_params <- data.frame(sound_speed = sound_speed_sw,
                              density = density_sw)
  # Define acoustic parameters =========================================
  acoustics <- data.frame(frequency = frequency,
                          k_sw = kcalc(frequency, sound_speed_sw))
  acoustics$k_b <- acoustics$k_sw * extract(object, "body")$h
  # Define object and model parameters ==================================
  length_body <- extract(object, "shape_parameters")$body$length
  model_params <- list(acoustics = acoustics,
                       parameters = data.frame(sound_speed_sw = sound_speed_sw,
                                               density_sw = density_sw))
  # Append to object ==================================
  slot(object, "model_parameters")$DWBA <- list(parameters = model_params,
                                                medium = medium_params,
                                                scatterer = extract(object, "body"))
  return(object)
}

#' Initialize GAS-object for modal series solution.
#' @param object GAS-class object.
#' @param radius_body Radius of sphere (m).
#' @param g_body Density contrast for gas.
#' @param h_body Sound speed contrast for gas.
#' @param ka_limit Modal series limit (i.e. max "m"). The default is the maximum
#'    ka + 10.
#' @inheritParams krm_initialize
#' @export
anderson_initialize <- function(object,
                                # Required input
                                frequency,
                                # Optional body parameters
                                radius_body = NULL,
                                g_body = NULL,
                                h_body = NULL,
                                # Optional medium inputs
                                sound_speed_sw = 1500,
                                density_sw = 1026,
                                # Optional modal series limit
                                ka_limit = NULL) {
  # Extract model body parameters =============================================
  body <- extract(object, "body")
  # Define medium parameters =========================================
  medium_params <- data.frame(sound_speed = sound_speed_sw,
                              density = density_sw)
  # Define acoustic parameters =========================================
  acoustics <- data.frame(frequency = frequency,
                          k_sw = kcalc(frequency, sound_speed_sw))
  acoustics$k_b <- acoustics$k_sw * extract(object, "body")$h
  ka_limit <- ifelse(!is.null(ka_limit),
                     ka_limit,
                     round(max(acoustics$k_sw) * body$radius) + 10)
  # Fill in remaining model inputs required by DCM ===========================
  radius <- ifelse(!is.null(radius_body), radius_body, body$radius)
  # Density contrast, g ======================================================
  g <- ifelse(!is.null(g), g, body$g)
  # Sound speed contrast, h ==================================================
  h <- ifelse(!is.null(h), h, body$h)
  # Define model dataframe for body parameterization =========================
  body_params <- data.frame(radius = radius,
                            diameter = radius * 2,
                            g = g,
                            h = h)
  # Tidy up model parameters to insert into object ===========================
  model_params <- list(body = body_params,
                       medium = medium_params,
                       acoustics = acoustics,
                       modal = data.frame(ka_limit = ka_limit))
  # Define object and model parameters =======================================
  slot(object, "model_parameters")$anderson <- model_params
  return(object)
}
#' Initialize object for Stanton high-pass approximation
#' @param object ESS-class object.
#' @param radius_shell Radius of shell.
#' @param g_shell Optional shell density contrast.
#' @param h_shell Optional shell sound speed contrast.
#' @inheritParams krm_initialize
#' @export
stanton_high_pass_initialize <- function(object,
                                         frequency,
                                         radius_shell = NULL,
                                         g_shell = NULL,
                                         h_shell = NULL,
                                         sound_speed_sw = 1500,
                                         density_sw = 1026) {
  # Extract model body parameters =============================================
  shell <- extract(object, "shell")
  # Define medium parameters ==================================================
  medium_params <- data.frame(sound_speed = sound_speed_sw,
                              density = density_sw)
  # Define acoustic parameters ================================================
  acoustics <- data.frame(frequency = frequency,
                          k_sw = kcalc(frequency, sound_speed_sw))
  acoustics$k_b <- acoustics$k_sw * extract(object, "shell")$h
  # Fill in remaining model inputs required by DCM ===========================
  radius <- ifelse(!is.null(radius_shell), radius_shell, shell$radius)
  # Density contrast, g ======================================================
  g <- ifelse(!is.null(g_shell), g_shell, shell$g)
  # Sound speed contrast, h ==================================================
  h <- ifelse(!is.null(h_shell), h_shell, shell$h)
  # Define model dataframe for body parameterization =========================
  shell_params <- data.frame(radius = radius,
                             diameter = radius * 2,
                             g = g,
                             h = h)
  # Tidy up model parameters to insert into object ===========================
  model_params <- list(shell = shell_params,
                       medium = medium_params,
                       acoustics = acoustics)
  # Define object and model parameters =======================================
  slot(object, "model_parameters")$stanton_high_pass <- model_params
  return(object)
}
