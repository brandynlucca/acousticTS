################################################################################
################################################################################
# PLOTTING AND DISPLAY UTILITIES
################################################################################
################################################################################
#' Collect model data for generic scatterer result plotting
#' @param object Scatterer object.
#' @param x_units Column name used for the x-axis.
#' @param y_units Column name used for the y-axis.
#' @return List containing split plot vectors and axis labels.
#' @keywords internal
#' @noRd
.collect_model_plot_data <- function(object,
                                     x_units = "frequency",
                                     y_units = "TS") {
  # Recover all model outputs currently stored on the scatterer ================
  models <- extract(object, "model")
  model_names <- names(models)

  if (length(model_names) == 0) {
    stop("ERROR: no model results detected in object.", call. = FALSE)
  }

  # Validate the requested plotting columns for each model output ==============
  required_cols <- unique(c("frequency", x_units, y_units))
  model_frames <- lapply(seq_along(model_names), function(i) {
    model_df <- models[[i]]
    missing_cols <- setdiff(required_cols, colnames(model_df))

    if (length(missing_cols) > 0) {
      stop(
        sprintf(
          "Model '%s' does not contain the requested plotting column(s): %s",
          model_names[i],
          paste(missing_cols, collapse = ", ")
        ),
        call. = FALSE
      )
    }

    model_df <- model_df[, required_cols, drop = FALSE]
    model_df$model <- model_names[i]
    model_df
  })

  # Combine the model outputs and build the plotting labels ====================
  models_df <- do.call("rbind", model_frames)
  models_df$model <- factor(models_df$model, levels = model_names)

  x_axis <- models_df[[x_units]]
  y_axis <- models_df[[y_units]]

  # Return the split vectors and axis labels used by the generic plotter =======
  list(
    x_axis = x_axis,
    y_axis = y_axis,
    x_mat = split(x_axis, models_df$model),
    y_mat = split(y_axis, models_df$model),
    model_names = model_names,
    x_lab = switch(x_units,
      frequency = "Frequency (Hz)",
      k_sw = expression(italic(k[sw] * a)),
      x_units
    ),
    y_lab = switch(y_units,
      TS = expression("Target" ~ "strength" ~ ("dB" ~ "re." ~ 1 ~ "m"^2)),
      sigma_bs = expression(sigma[bs]),
      y_units
    )
  )
}

#' Detect whether the current graphics device is using a multi-panel layout
#' @keywords internal
#' @noRd
.has_active_multi_panel_layout <- function() {
  layout_dims <- graphics::par("mfg")[3:4]
  prod(layout_dims) > 1
}

#' Plot model-result curves for scatterers with multiple model outputs
#' @param plot_data Output from `.collect_model_plot_data()`.
#' @param nudge_y y-axis nudge factor.
#' @param nudge_x x-axis nudge factor.
#' @param mar Plot margins.
#' @keywords internal
#' @noRd
.plot_model_results <- function(plot_data,
                                nudge_y = 1.05,
                                nudge_x = 1.01,
                                mar = c(5.0, 5.5, 2.0, 3.5)) {
  # Initialize the plotting device and empty plotting window ===================
  if (.has_active_multi_panel_layout()) {
    graphics::par(ask = FALSE)
  } else {
    graphics::par(ask = FALSE, mar = mar)
  }

  graphics::plot(
    x = seq(
      from = min(plot_data$x_axis),
      to = max(plot_data$x_axis),
      length.out = 2
    ),
    y = seq(
      from = min(plot_data$y_axis),
      to = max(plot_data$y_axis),
      length.out = 2
    ),
    xlim = c(
      min(plot_data$x_axis) * (1 - nudge_x),
      max(plot_data$x_axis) * nudge_x
    ),
    ylim = c(
      min(plot_data$y_axis) * (1 - (1 - nudge_y)),
      max(plot_data$y_axis) * (1 + (1 - nudge_y))
    ),
    xlab = plot_data$x_lab,
    ylab = plot_data$y_lab,
    yaxs = "i",
    xaxs = "i",
    xaxt = "n",
    type = "n",
    cex.axis = 1.3,
    cex.lab = 1.5
  )

  # Draw the x-axis ticks explicitly to avoid scientific notation noise ========
  atx <- seq(
    graphics::par("xaxp")[1], graphics::par("xaxp")[2],
    (graphics::par("xaxp")[2] - graphics::par("xaxp")[1]) /
      graphics::par("xaxp")[3]
  )
  graphics::axis(
    1,
    at = atx,
    labels = format(atx, scientific = FALSE),
    cex.axis = 1.3
  )

  # Draw one result curve per stored model =====================================
  invisible(mapply(
    graphics::lines,
    plot_data$x_mat,
    plot_data$y_mat,
    col = model_palette[seq_len(length(plot_data$model_names))],
    lwd = 4
  ))

  # Add the model legend after drawing all curves ==============================
  graphics::legend(
    "bottomright",
    title = expression(bold("TS" ~ "model")),
    title.adj = 0.05,
    legend = plot_data$model_names,
    lty = rep(1, length(plot_data$model_names)),
    col = model_palette[seq_len(length(plot_data$model_names))],
    cex = 1.05,
    lwd = 4
  )
}

#' Format named material properties as indented display lines
#'
#' @param properties Named list or vector of properties.
#' @param label_map Optional named character vector mapping property names to
#'   display labels.
#' @param unit_map Optional named character vector mapping property names to
#'   display units.
#' @param digits Number of digits passed to `round()`.
#' @param indent Leading indentation string.
#' @return Character vector of formatted lines.
#' @keywords internal
#' @noRd
.format_property_lines <- function(properties,
                                   label_map = character(),
                                   unit_map = character(),
                                   digits = 4L,
                                   indent = "    ") {
  # Return early when there are no properties to format ========================
  if (length(properties) == 0) {
    return(character())
  }

  # Format each named property with the requested labels and units =============
  vapply(names(properties), function(name) {
    value <- properties[[name]]
    label <- if (name %in% names(label_map)) {
      label_map[[name]]
    } else {
      gsub("_", " ", name)
    }
    units <- if (name %in% names(unit_map)) unit_map[[name]] else ""
    paste0(indent, label, ": ", round(value, digits), units)
  }, character(1))
}

#' Build centered x/radius data for axisymmetric outline plotting
#' @param position_matrix Column-major position matrix.
#' @param shape_parameters Optional shape-parameter list.
#' @param center_x Logical; whether to center x around its midpoint.
#' @param error_context Character context for radius recovery.
#' @keywords internal
#' @noRd
.axisymmetric_outline_data <- function(position_matrix,
                                       shape_parameters = NULL,
                                       center_x = FALSE,
                                       error_context = "shape") {
  # Recover the centered x coordinates and nodewise radius profile =============
  x <- .shape_x(position_matrix, row_major = FALSE)
  radius <- .shape_radius_profile(
    position_matrix = position_matrix,
    shape_parameters = shape_parameters,
    row_major = FALSE,
    error_context = error_context
  )

  # Center the x axis around its midpoint when requested =======================
  if (center_x) {
    x <- x - stats::median(x, na.rm = TRUE)
  }

  # Return the outline coordinates used by downstream plotting helpers =========
  list(x = x, radius = radius)
}

#' Plot an axisymmetric outline from a column-major position matrix
#' @param position_matrix Column-major position matrix.
#' @param shape_parameters Optional shape-parameter list.
#' @param xlab x-axis label.
#' @param ylab y-axis label.
#' @param nudge_y y-axis nudge factor.
#' @param lwd Line width.
#' @param col Line color.
#' @param segment_col Segment color.
#' @param segment_lwd Segment line width.
#' @param center_x Logical; whether to center the x-axis.
#' @param init Logical; when `TRUE`, initialize the plot; otherwise add to the
#'   current device.
#' @param mar Plot margins.
#' @param cex_lab Axis-label expansion.
#' @param cex_axis Axis tick-label expansion.
#' @param ... Additional arguments passed to `graphics::plot()`.
#' @keywords internal
#' @noRd
.plot_axisymmetric_outline <- function(position_matrix,
                                       shape_parameters = NULL,
                                       xlab = "Semi-major axis (m)",
                                       ylab = "Semi-minor radius (m)",
                                       nudge_y = 1.01,
                                       lwd = 4,
                                       col = "black",
                                       segment_col = col,
                                       segment_lwd = 1.25,
                                       center_x = FALSE,
                                       init = TRUE,
                                       mar = c(3, 4.5, 1, 2),
                                       cex_lab = 1.2,
                                       cex_axis = 1.2,
                                       ...) {
  # Build the outline coordinates for the supplied axisymmetric geometry =======
  outline <- .axisymmetric_outline_data(
    position_matrix = position_matrix,
    shape_parameters = shape_parameters,
    center_x = center_x
  )

  # Initialize the plot when requested =========================================
  if (init) {
    if (.has_active_multi_panel_layout()) {
      graphics::par(ask = FALSE)
    } else {
      graphics::par(
        ask = FALSE,
        oma = c(1, 1, 1, 0),
        mar = mar
      )
    }
    graphics::plot(
      x = outline$x,
      y = outline$radius,
      type = "l",
      ylab = ylab,
      xlab = xlab,
      lwd = lwd,
      cex.lab = cex_lab,
      cex.axis = cex_axis,
      ylim = c(
        min(-outline$radius, na.rm = TRUE) * (1 - (1 - nudge_y)),
        max(outline$radius, na.rm = TRUE) * nudge_y
      ),
      col = col,
      ...
    )
  } else {
    # Add the upper outline onto an existing plotting device ===================
    graphics::lines(
      x = outline$x,
      y = outline$radius,
      lty = 1,
      lwd = lwd,
      col = col
    )
  }

  # Draw the mirrored lower outline and the segment guide lines ================
  graphics::lines(
    x = outline$x,
    y = -outline$radius,
    lty = 1,
    lwd = lwd,
    col = col
  )
  graphics::segments(
    x0 = outline$x,
    x1 = outline$x,
    y0 = -outline$radius,
    y1 = outline$radius,
    lty = 3,
    lwd = segment_lwd,
    col = segment_col
  )

  # Return the outline coordinates invisibly ===================================
  invisible(outline)
}

#' Plot a segmented row-major body outline using the canonical geometry helpers
#' @param rpos Row-major position matrix.
#' @param shape_parameters Optional shape-parameter list.
#' @param nudge_y y-axis nudge factor.
#' @param aspect_ratio Plot aspect-ratio mode.
#' @param xlab x-axis label.
#' @param ylab y-axis label.
#' @keywords internal
#' @noRd
.plot_row_major_segmented_body <- function(rpos,
                                           shape_parameters = NULL,
                                           nudge_y = 1.05,
                                           aspect_ratio = "manual",
                                           xlab = "Length (m)",
                                           ylab = "Thickness (m)") {
  # Canonicalize the row-major geometry and recover its key profiles ===========
  rpos <- .canonicalize_position_matrix(rpos, row_major = TRUE)
  x <- .shape_x(rpos, row_major = TRUE)
  radius <- .shape_radius_profile(
    position_matrix = rpos,
    shape_parameters = shape_parameters,
    row_major = TRUE,
    error_context = "segmented body plot"
  )
  zU <- .geometry_axis_values(
    rpos,
    axis = "zU",
    row_major = TRUE,
    default = NULL,
    context = "segmented body plot"
  )
  zL <- .geometry_axis_values(
    rpos,
    axis = "zL",
    row_major = TRUE,
    default = NULL,
    context = "segmented body plot"
  )
  z_center <- if (!is.null(rownames(rpos)) &&
      any(c("z", "z_body", "z_bladder") %in% rownames(rpos))) {
    .extract_shape_component_row(
      rpos,
      c("z", "z_body", "z_bladder"),
      default = rep(0, ncol(rpos))
    )
  } else if (!is.null(zU) && !is.null(zL)) {
    (zU + zL) / 2
  } else {
    rep(0, length(x))
  }

  # Resolve the plotting limits for the requested aspect-ratio mode ============
  if (aspect_ratio == "manual") {
    vert_lims <- c(
      min(z_center - radius, na.rm = TRUE) * (1 - (1 - nudge_y)),
      max(z_center + radius, na.rm = TRUE) * nudge_y
    )
  } else {
    max_x <- max(abs(x), na.rm = TRUE)
    vert_lims <- c(-max_x * 0.10, max_x * 0.10)
  }

  # Initialize the plot before drawing each body segment polygon ===============
  if (.has_active_multi_panel_layout()) {
    graphics::par(ask = FALSE)
  } else {
    graphics::par(
      ask = FALSE,
      oma = c(1, 1, 1, 0),
      mar = c(5.0, 4.5, 1.5, 2)
    )
  }
  graphics::plot(
    x,
    z_center,
    type = "n",
    xlab = xlab,
    ylab = ylab,
    ylim = vert_lims
  )

  # Draw one quadrilateral per body segment using the local normal direction ===
  for (i in seq_len(length(x) - 1L)) {
    x0 <- x[i]
    y0 <- z_center[i]
    x1 <- x[i + 1L]
    y1 <- z_center[i + 1L]
    r0 <- radius[i]
    r1 <- radius[i + 1L]

    if ((!is.finite(r0) || r0 < 0) && (!is.finite(r1) || r1 < 0)) {
      next
    }

    if ((!is.finite(r0) || r0 <= 0) && (!is.finite(r1) || r1 <= 0)) {
      next
    }

    dx <- x1 - x0
    dy <- y1 - y0
    len <- sqrt(dx^2 + dy^2)
    if (!is.finite(len) || len <= 0) {
      next
    }

    px <- -dy / len
    py <- dx / len
    r0 <- if (is.finite(r0) && r0 > 0) r0 else 0
    r1 <- if (is.finite(r1) && r1 > 0) r1 else 0
    xA <- x0 + r0 * px
    yA <- y0 + r0 * py
    xB <- x1 + r1 * px
    yB <- y1 + r1 * py
    xC <- x1 - r1 * px
    yC <- y1 - r1 * py
    xD <- x0 - r0 * px
    yD <- y0 - r0 * py

    graphics::polygon(
      x = c(xA, xB, xC, xD),
      y = c(yA, yB, yC, yD),
      col = grDevices::adjustcolor("gray50", alpha.f = 0.6),
      border = "black",
      lwd = 1
    )
  }

  # Overlay the segmented centerline and node markers ==========================
  graphics::lines(x, z_center, lwd = 3, col = "gray90")
  graphics::points(x, z_center, pch = 1, col = "black", cex = 0.8)
}

#' Resolve plotting fields for a row-major profile component
#' @param rpos Row-major profile matrix.
#' @param radius Optional explicit radius profile.
#' @keywords internal
#' @noRd
.component_profile_outline <- function(rpos, radius = NULL) {
  # Recover the canonical x/z/w fields for the requested component =============
  fields <- .profile_fields(rpos)

  # Fall back to the radius-derived width when no width profile is present =====
  if (is.null(radius)) {
    radius <- .shape_radius_profile(
      position_matrix = rpos,
      row_major = TRUE,
      error_context = "profile component"
    )
  }

  width <- fields$w
  if (all(is.na(width)) || max(abs(width), na.rm = TRUE) <=
      sqrt(.Machine$double.eps)) {
    width <- 2 * radius
  }

  # Return the outline fields used by the two-panel component plot =============
  list(
    x = fields$x,
    zU = fields$zU,
    zL = fields$zL,
    half_width = width / 2
  )
}

#' Plot the shared two-panel shape layout used by profile-based composite targets
#' @param primary_rpos Primary row-major profile matrix.
#' @param primary_radius Optional explicit primary radius profile.
#' @param secondary_rpos Optional secondary row-major profile matrix.
#' @param secondary_radius Optional explicit secondary radius profile.
#' @param primary_label Primary legend label.
#' @param secondary_label Secondary legend label.
#' @param secondary_col Secondary plotting color.
#' @param nudge_y y-axis nudge factor.
#' @param nudge_x x-axis nudge factor.
#' @keywords internal
#' @noRd
.plot_profile_component_panels <- function(
    primary_rpos,
    primary_radius = NULL,
    secondary_rpos = NULL,
    secondary_radius = NULL,
    primary_label = "Body",
    secondary_label = "Internal component",
    secondary_col = "red",
    nudge_y = 1.05,
    nudge_x = 1.01
  ) {
  # Recover the outline data for the primary and optional secondary component ==
  primary <- .component_profile_outline(primary_rpos, radius = primary_radius)
  secondary <- if (!is.null(secondary_rpos)) {
    .component_profile_outline(secondary_rpos, radius = secondary_radius)
  } else {
    NULL
  }

  # Resolve the plotting limits shared by both component panels ================
  ylow_lim <- min(primary$zL, na.rm = TRUE) * (1 - (1 - nudge_y))
  yhi_lim <- max(primary$zU, na.rm = TRUE) * nudge_y
  xlow_lim <- min(primary$x, na.rm = TRUE) * (1 - (1 - nudge_x))
  xhi_lim <- max(primary$x, na.rm = TRUE) * nudge_x

  # Initialize the shared two-panel layout =====================================
  graphics::par(
    ask = FALSE,
    mfrow = c(2, 1),
    mar = c(2, 4.5, 1, 2),
    oma = c(1.5, 0.5, 0.5, 0)
  )

  # Draw the height panel and optional secondary component =====================
  graphics::plot(
    x = primary$x,
    y = primary$zU,
    type = "l",
    ylim = c(ylow_lim, yhi_lim),
    xlim = c(xlow_lim, xhi_lim),
    xlab = "",
    ylab = "Height (m)",
    lwd = 4,
    xaxt = "n",
    cex.lab = 1.2,
    cex.axis = 1.2
  )
  graphics::lines(primary$x, primary$zL, lty = 1, lwd = 4)
  graphics::segments(
    x0 = primary$x, x1 = primary$x,
    y0 = primary$zL, y1 = primary$zU,
    lty = 3, lwd = 1, col = "gray30"
  )

  if (!is.null(secondary)) {
    graphics::lines(
      secondary$x, secondary$zL,
      lty = 1, lwd = 3.5, col = secondary_col
    )
    graphics::lines(
      secondary$x, secondary$zU,
      lty = 1, lwd = 3.5, col = secondary_col
    )
    graphics::segments(
      x0 = secondary$x, x1 = secondary$x,
      y0 = secondary$zL, y1 = secondary$zU,
      lty = 3, lwd = 1, col = secondary_col
    )
    graphics::legend(
      x = "bottomright",
      legend = c(primary_label, secondary_label),
      lty = c(1, 1),
      lwd = c(4, 3.5),
      col = c("black", secondary_col),
      cex = 0.95
    )
  }

  # Draw the width panel and optional secondary component ======================
  left_limit <- -max(primary$half_width, na.rm = TRUE) * (1 - (1 - nudge_y))
  right_limit <- max(primary$half_width, na.rm = TRUE) * (1 - (1 - nudge_y))

  graphics::plot(
    x = primary$x,
    y = primary$half_width,
    type = "l",
    ylim = c(left_limit, right_limit),
    ylab = "Width (m)",
    lwd = 4,
    cex.lab = 1.2,
    cex.axis = 1.2
  )
  graphics::lines(primary$x, -primary$half_width, lty = 1, lwd = 4)

  if (!is.null(secondary)) {
    graphics::lines(
      secondary$x, secondary$half_width,
      lty = 1, lwd = 3.5, col = secondary_col
    )
    graphics::lines(
      secondary$x, -secondary$half_width,
      lty = 1, lwd = 3.5, col = secondary_col
    )
  }

  graphics::mtext("Along-body axis (m)",
    side = 1,
    outer = TRUE,
    line = 0,
    cex = 1.2
  )
}
