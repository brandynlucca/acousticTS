################################################################################
# Transition matrix method (TMM) plotting helpers
################################################################################

# Normalize the user-facing TMM post-processing quantity selection while
# keeping the older `level_dB` name as a compatibility alias.
#' @noRd
.tmm_plot_quantity <- function(quantity) {
  # Fall back to the default scattering-level display quantity =================
  if (missing(quantity) || is.null(quantity)) {
    return("sigma_scat_dB")
  }
  if (identical(quantity, "level_dB")) {
    return("sigma_scat_dB")
  }
  # Validate the requested plotting quantity ===================================
  match.arg(quantity, c("sigma_scat_dB", "sigma_scat", "mod_f", "phase"))
}
# Build one angular scattering slice from the stored TMM blocks.
#' @noRd
.tmm_scattering_slice_data <- function(object,
                                       frequency = NULL,
                                       vary = "theta_scatter",
                                       quantity = "sigma_scat_dB",
                                       n_points = 181,
                                       theta_body = NULL,
                                       phi_body = NULL,
                                       theta_scatter = NULL,
                                       phi_scatter = NULL) {
  # Recover the stored TMM state and the requested frequency slice =============
  model_params <- .tmm_require_stored_blocks(object)
  parameters <- model_params$parameters
  defaults <- model_params$body
  acoustics <- parameters$acoustics
  shape_parameters <- acousticTS::extract(object, "shape_parameters")
  idx <- .tmm_plot_frequency_index(frequency, acoustics$frequency)

  if (!is.numeric(n_points) || length(n_points) != 1 || !is.finite(n_points) ||
      n_points < 2 || n_points %% 1 != 0) {
    stop("'n_points' must be a single integer >= 2.", call. = FALSE)
  }

  # Resolve the varying angle and normalize the plotting quantity ==============
  vary <- match.arg(vary, c("theta_scatter", "phi_scatter", "theta_body",
                            "phi_body"))
  quantity <- .tmm_plot_quantity(quantity)

  theta_body_default <- defaults$theta_body
  phi_body_default <- defaults$phi_body %||% pi

  grid <- if (vary %in% c("theta_body", "theta_scatter")) {
    seq(0, pi, length.out = n_points)
  } else {
    seq(0, 2 * pi, length.out = n_points)
  }

  # Expand the incident and receive geometry over the slice grid ===============
  theta_body_vec <- rep(.tmm_scalar_angle(theta_body, theta_body_default,
                                          "theta_body"), n_points)
  phi_body_vec <- rep(.tmm_scalar_angle(phi_body, phi_body_default,
                                        "phi_body"), n_points)
  theta_scatter_vec <- if (is.null(theta_scatter)) {
    pi - theta_body_vec
  } else {
    rep(.tmm_scalar_angle(theta_scatter, pi - theta_body_vec[1],
                          "theta_scatter"), n_points)
  }
  phi_scatter_vec <- if (is.null(phi_scatter)) {
    phi_body_vec + pi
  } else {
    rep(.tmm_scalar_angle(phi_scatter, phi_body_vec[1] + pi,
                          "phi_scatter"), n_points)
  }

  if (vary == "theta_body") {
    theta_body_vec <- grid
    if (is.null(theta_scatter)) {
      theta_scatter_vec <- pi - grid
    }
  } else if (vary == "phi_body") {
    phi_body_vec <- grid
    if (is.null(phi_scatter)) {
      phi_scatter_vec <- grid + pi
    }
  } else if (vary == "theta_scatter") {
    theta_scatter_vec <- grid
  } else if (vary == "phi_scatter") {
    phi_scatter_vec <- grid
  }

  # Evaluate the retained operator and map it onto the plotted quantity ========
  f_vals <- .tmm_scattering_points(
    model_params = model_params,
    frequency_idx = idx,
    shape_parameters = shape_parameters,
    theta_body = theta_body_vec,
    phi_body = phi_body_vec,
    theta_scatter = theta_scatter_vec,
    phi_scatter = phi_scatter_vec
  )

  y_vals <- switch(
    quantity,
    sigma_scat_dB = db(.sigma_bs(f_vals)),
    sigma_scat = .sigma_bs(f_vals),
    mod_f = Mod(f_vals),
    phase = Arg(f_vals)
  )

  # Resolve axis labels for the selected angular slice =========================
  x_lab <- switch(
    vary,
    theta_body = expression(theta[body] ~ "(rad)"),
    phi_body = expression(phi[body] ~ "(rad)"),
    theta_scatter = expression(theta[scatter] ~ "(rad)"),
    phi_scatter = expression(phi[scatter] ~ "(rad)")
  )
  y_lab <- switch(
    quantity,
    sigma_scat_dB = expression(10 * log[10](sigma[scat] / (1 ~ m^2)) ~ "(dB)"),
    sigma_scat = expression(sigma[scat]),
    mod_f = expression("|" * f[scat] * "|"),
    phase = "Scattering phase (rad)"
  )

  # Return the fully prepared slice data =======================================
  list(
    frequency = acoustics$frequency[idx],
    x = grid,
    y = y_vals,
    x_lab = x_lab,
    y_lab = y_lab,
    vary = vary,
    quantity = quantity
  )
}

# Resolve finite plotting data and labels for a stored TMM scattering grid.
#' @noRd
.tmm_grid_plot_data <- function(grid, quantity = "sigma_scat_dB") {
  # Normalize the requested plotted quantity ===================================
  quantity <- .tmm_plot_quantity(quantity)
  z_vals <- switch(
    quantity,
    sigma_scat_dB = grid$sigma_scat_dB,
    sigma_scat = grid$sigma_scat,
    mod_f = Mod(grid$f_scat),
    phase = Arg(grid$f_scat)
  )
  z_lab <- switch(
    quantity,
    sigma_scat_dB = expression(10 * log[10](sigma[scat] / (1 ~ m^2)) ~ "(dB)"),
    sigma_scat = expression(sigma[scat]),
    mod_f = expression("|" * f[scat] * "|"),
    phase = "Scattering phase (rad)"
  )

  # Return the grid values and legend label used by the plotting helpers =======
  list(z = z_vals, label = z_lab, quantity = quantity)
}

# Build a symmetric grid-edge vector for cell-based plotting.
#' @noRd
.plot_tmm_scattering_heatmap <- function(grid, quantity = "sigma_scat_dB") {
  # Resolve the plotted grid values and legend label ===========================
  plot_data <- .tmm_grid_plot_data(grid, quantity = quantity)

  # Draw the filled contour heatmap with TMM-specific axis labels ==============
  graphics::filled.contour(
    x = grid$phi_scatter,
    y = grid$theta_scatter,
    z = t(plot_data$z),
    xlab = "",
    ylab = "",
    key.title = graphics::title(main = plot_data$label, cex.main = 0.82),
    plot.axes = {
      graphics::axis(1)
      graphics::axis(2)
    },
    plot.title = graphics::title(
      main = sprintf("TMM scattering grid at %.1f kHz", grid$frequency * 1e-3),
      xlab = expression(phi[scatter] ~ "(rad)"),
      ylab = expression(theta[scatter] ~ "(rad)")
    ),
    color.palette = function(n) grDevices::hcl.colors(n, "Spectral", rev = TRUE)
  )

  # Return the input grid invisibly for chaining ===============================
  invisible(grid)
}

# Plot a polar-style scattering map from a stored TMM grid.
#' @noRd
.plot_tmm_scattering_polar <- function(grid, quantity = "sigma_scat_dB") {
  # Resolve the plotted quantity and finite dynamic range ======================
  plot_data <- .tmm_grid_plot_data(grid, quantity = quantity)
  z_vals <- plot_data$z
  z_finite <- z_vals[is.finite(z_vals)]

  if (!length(z_finite)) {
    stop(
      "The requested TMM scattering grid does not contain finite plotted ",
      "values.",
      call. = FALSE
    )
  }

  # Build the polar cell geometry, palette, and legend scale ===================
  theta_edges <- .tmm_grid_edges(grid$theta_scatter, lower = 0, upper = pi)
  phi_edges <- .tmm_grid_edges(grid$phi_scatter, lower = 0, upper = 2 * pi)
  palette <- grDevices::hcl.colors(128, "Spectral", rev = TRUE)
  z_breaks <- seq(min(z_finite), max(z_finite),
                  length.out = length(palette) + 1)
  opar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(opar))
  graphics::par(
    ask = FALSE,
    fig = c(0, 0.82, 0, 1),
    mar = c(2.0, 2.0, 2.4, 1.1)
  )
  graphics::plot(
    NA,
    NA,
    xlim = c(-pi, pi),
    ylim = c(-pi, pi),
    asp = 1,
    axes = FALSE,
    xlab = "",
    ylab = "",
    xaxs = "i",
    yaxs = "i"
  )
  graphics::mtext(
    sprintf("TMM polar scattering map at %.1f kHz", grid$frequency * 1e-3),
    side = 3,
    line = 0.35,
    font = 2
  )

  # Fill each polar cell with the corresponding scattering value ===============
  for (j in seq_len(length(grid$phi_scatter))) {
    for (i in seq_len(length(grid$theta_scatter))) {
      z_ij <- z_vals[i, j]
      if (!is.finite(z_ij)) {
        next
      }

      col_idx <- findInterval(z_ij, z_breaks, all.inside = TRUE)
      r_inner <- theta_edges[i]
      r_outer <- theta_edges[i + 1]
      phi_left <- phi_edges[j]
      phi_right <- phi_edges[j + 1]

      x_poly <- c(
        r_inner * cos(phi_left),
        r_outer * cos(phi_left),
        r_outer * cos(phi_right),
        r_inner * cos(phi_right)
      )
      y_poly <- c(
        r_inner * sin(phi_left),
        r_outer * sin(phi_left),
        r_outer * sin(phi_right),
        r_inner * sin(phi_right)
      )

      graphics::polygon(
        x_poly,
        y_poly,
        col = palette[col_idx],
        border = NA
      )
    }
  }

  # Overlay the polar guide circles, spokes, and angular labels ================
  for (r in seq(pi / 4, pi, by = pi / 4)) {
    graphics::symbols(
      0,
      0,
      circles = r,
      inches = FALSE,
      add = TRUE,
      fg = "grey75",
      bg = NA
    )
  }
  for (phi_ref in c(0, pi / 2, pi, 3 * pi / 2)) {
    graphics::segments(
      0,
      0,
      pi * cos(phi_ref),
      pi * sin(phi_ref),
      col = "grey75",
      lty = 3
    )
  }
  graphics::symbols(0, 0, circles = pi, inches = FALSE, add = TRUE,
                    fg = "black", bg = NA)
  graphics::text(0, 0, labels = "0", cex = 0.85, font = 2)
  graphics::text(
    x = c(pi + 0.18, 0, -pi - 0.18, 0),
    y = c(0, pi - 0.12, 0, -pi + 0.12),
    labels = c(expression(0), expression(pi / 2), expression(pi),
               expression(3 * pi / 2)),
    xpd = NA,
    cex = 0.9
  )
  graphics::text(
    x = c(pi / 4, pi / 2, 3 * pi / 4, pi),
    y = 0,
    labels = c(expression(pi / 4), expression(pi / 2),
               expression(3 * pi / 4), expression(pi)),
    pos = 3,
    col = "grey35",
    cex = 0.8
  )

  # Draw the companion colorbar for the plotted scattering quantity ============
  graphics::par(fig = c(0.845, 0.965, 0.15, 0.88),
                mar = c(1.8, 0.2, 2.2, 3.0),
                new = TRUE)
  graphics::plot.new()
  graphics::plot.window(xlim = c(0, 1), ylim = range(z_breaks))
  y_edges <- seq(min(z_breaks), max(z_breaks), length.out = length(palette) + 1)
  for (i in seq_along(palette)) {
    graphics::rect(
      xleft = 0,
      ybottom = y_edges[i],
      xright = 1,
      ytop = y_edges[i + 1],
      col = palette[i],
      border = NA
    )
  }
  graphics::axis(4, las = 1)
  graphics::mtext(plot_data$label, side = 4, line = 2.05, cex = 0.8)

  # Return the input grid invisibly for chaining ===============================
  invisible(grid)
}

# Plot a stored TMM angular scattering slice.
#' @noRd
.plot_tmm_scattering_slice <- function(object,
                                       nudge_x = 1.01,
                                       nudge_y = 1.05,
                                       polar = FALSE,
                                       heatmap = FALSE,
                                       ...) {
  # Collect plotting options passed through from the public plot method ========
  extra_args <- list(...)

  if (isTRUE(polar) && isTRUE(heatmap)) {
    stop(
      "Only one of 'polar' or 'heatmap' can be TRUE for TMM scattering plots.",
      call. = FALSE
    )
  }

  # Route map-style requests through the 2D grid plotting helpers ==============
  if (isTRUE(polar) || isTRUE(heatmap)) {
    grid <- do.call(tmm_scattering_grid, c(list(object = object), extra_args))
    if (isTRUE(polar)) {
      return(invisible(.plot_tmm_scattering_polar(
        grid,
        quantity = extra_args$quantity %||% "sigma_scat_dB"
      )))
    }
    return(invisible(.plot_tmm_scattering_heatmap(
      grid,
      quantity = extra_args$quantity %||% "sigma_scat_dB"
    )))
  }

  # Build the one-dimensional slice and its plotting limits ====================
  slice <- do.call(.tmm_scattering_slice_data, c(list(object = object),
                                                 extra_args))
  y_finite <- slice$y[is.finite(slice$y)]

  if (!length(y_finite)) {
    stop(
      "The requested scattering slice does not contain any finite plotted ",
      "values.",
      call. = FALSE
    )
  }

  y_lim <- range(y_finite)
  if (diff(y_lim) == 0) {
    pad <- max(1e-12, abs(y_lim[1]) * 0.05)
    y_lim <- y_lim + c(-pad, pad)
  } else {
    y_lim <- c(
      y_lim[1] * (1 - (1 - nudge_y)),
      y_lim[2] * (1 + (1 - nudge_y))
    )
  }

  # Draw the angular scattering slice with the resolved axis labels ============
  graphics::par(ask = FALSE, mar = c(5.0, 5.5, 2.0, 2.0))
  graphics::plot(
    x = slice$x,
    y = slice$y,
    type = "l",
    lwd = 3,
    xlab = slice$x_lab,
    ylab = slice$y_lab,
    xlim = c(min(slice$x) * (1 - (1 - nudge_x)), max(slice$x) * nudge_x),
    ylim = y_lim,
    xaxs = "i",
    yaxs = "i"
  )
  graphics::title(main = sprintf("TMM scattering slice at %.1f kHz",
                                 slice$frequency * 1e-3))

  # Return the slice data invisibly for downstream inspection ==================
  invisible(slice)
}
