
################################################################################
# Transition matrix method (TMM) diagnostics helpers
################################################################################

# Resolve zero, one, or many stored frequencies for diagnostics-style helpers.
#' @noRd
.tmm_frequency_indices <- function(frequency, available) {
  # Default to all stored frequencies ==========================================
  if (is.null(frequency)) {
    return(seq_along(available))
  }
  # Validate the requested frequencies =========================================
  if (!is.numeric(frequency) || !length(frequency) || any(!is.finite(frequency))) {
    stop("'frequency' must be NULL or a numeric vector of finite values in Hz.", call. = FALSE)
  }

  idx <- vapply(
    frequency,
    function(f) .tmm_plot_frequency_index(f, available),
    integer(1)
  )
  # Return the resolved stored-frequency indices ===============================
  as.integer(idx)
}

# Deterministic angle pairs used by the reciprocity diagnostic when the user
# does not supply a custom set.
#' @noRd
.tmm_default_reciprocity_pairs <- function() {
  data.frame(
    theta_body = c(pi / 6, pi / 3, pi / 2, 2 * pi / 3),
    phi_body = c(pi / 8, pi / 3, pi / 2, 5 * pi / 4),
    theta_scatter = c(pi / 4, pi / 2, 2 * pi / 3, pi / 3),
    phi_scatter = c(pi / 2, 7 * pi / 6, 0.8 * pi, 11 * pi / 6)
  )
}

# Validate an optional user-supplied reciprocity-angle table.
#' @noRd
.tmm_validate_reciprocity_pairs <- function(reciprocity_pairs) {
  # Fall back to the package-default reciprocity pairs =========================
  if (is.null(reciprocity_pairs)) {
    return(.tmm_default_reciprocity_pairs())
  }
  # Validate the required columns and angle values =============================
  required_cols <- c("theta_body", "phi_body", "theta_scatter", "phi_scatter")
  missing_cols <- setdiff(required_cols, names(reciprocity_pairs))
  if (length(missing_cols)) {
    stop(
      "'reciprocity_pairs' must contain the columns: ",
      paste(required_cols, collapse = ", "),
      ". Missing: ",
      paste(missing_cols, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  for (nm in required_cols) {
    if (!is.numeric(reciprocity_pairs[[nm]]) || any(!is.finite(reciprocity_pairs[[nm]]))) {
      stop("'", nm, "' in 'reciprocity_pairs' must be finite numeric angles in radians.", call. = FALSE)
    }
  }

  reciprocity_pairs
}

# Equivalent-volume sphere radius for the supported axisymmetric TMM shapes.
#' @noRd
.tmm_equivalent_volume_radius <- function(shape_parameters) {
  # Return the exact sphere radius when already spherical ======================
  shape <- shape_parameters$shape
  if (identical(shape, "Sphere")) {
    return(shape_parameters$radius)
  }

  if (identical(shape, "ProlateSpheroid")) {
    semi_major <- shape_parameters$length / 2
    semi_minor <- shape_parameters$radius
    return((semi_major * semi_minor^2)^(1 / 3))
  }

  if (identical(shape, "OblateSpheroid")) {
    semi_equatorial <- shape_parameters$radius
    semi_polar <- shape_parameters$length / 2
    return((semi_equatorial^2 * semi_polar)^(1 / 3))
  }

  NA_real_
}

# Target aspect ratio for the supported spheroidal TMM branches.
#' @noRd
.tmm_spheroidal_aspect_ratio <- function(shape_parameters) {
  # Return the prolate aspect ratio ============================================
  if (identical(shape_parameters$shape, "ProlateSpheroid")) {
    return((shape_parameters$length / 2) / shape_parameters$radius)
  }
  if (identical(shape_parameters$shape, "OblateSpheroid")) {
    return(shape_parameters$radius / (shape_parameters$length / 2))
  }
  # Non-spheroidal shapes have no spheroidal aspect ratio ======================
  NA_real_
}

# Convert an equal-volume sphere radius and target aspect ratio into the shape
# parameters needed by the public spheroid generators.
#' @noRd
.tmm_spheroid_axes_from_equal_volume <- function(shape_type,
                                                 radius_eq,
                                                 aspect_ratio) {
  # Convert the equal-volume radius into prolate axes ==========================
  if (identical(shape_type, "ProlateSpheroid")) {
    semi_minor <- radius_eq / aspect_ratio^(1 / 3)
    semi_major <- radius_eq * aspect_ratio^(2 / 3)
    return(list(length_body = 2 * semi_major, radius_body = semi_minor))
  }

  if (identical(shape_type, "OblateSpheroid")) {
    semi_equatorial <- radius_eq * aspect_ratio^(1 / 3)
    semi_polar <- radius_eq / aspect_ratio^(2 / 3)
    return(list(length_body = 2 * semi_polar, radius_body = semi_equatorial))
  }
  # Reject unsupported target shapes ===========================================
  stop("Unsupported shape for spheroidal continuation.", call. = FALSE)
}

# Rebuild one fluid- or gas-like scatterer around a new canonical shape while
# preserving the original material specification.
#' @noRd
.tmm_rebuild_shape_like <- function(object, shape, body) {
  # Rebuild a gas-filled surrogate when needed =================================
  if (methods::is(object, "GAS")) {
    args <- list(shape = shape, theta_body = body$theta %||% (pi / 2))
    if (!is.null(body$density)) {
      args$density_fluid <- body$density
    } else if (!is.null(body$g)) {
      args$g_fluid <- body$g
    }
    if (!is.null(body$sound_speed)) {
      args$sound_speed_fluid <- body$sound_speed
    } else if (!is.null(body$h)) {
      args$h_fluid <- body$h
    }
    return(do.call(gas_generate, args))
  }
  # Otherwise rebuild a fluid-like surrogate ===================================
  args <- list(shape = shape, theta_body = body$theta %||% (pi / 2))
  if (!is.null(body$density)) {
    args$density_body <- body$density
  } else if (!is.null(body$g)) {
    args$g_body <- body$g
  }
  if (!is.null(body$sound_speed)) {
    args$sound_speed_body <- body$sound_speed
  } else if (!is.null(body$h)) {
    args$h_body <- body$h
  }
  do.call(fls_generate, args)
}

# Build the equal-volume sphere-to-spheroid continuation path for the stored
# target so aspect-ratio drift can be checked as an internal TMM sanity test.
#' @noRd
.tmm_sphere_to_spheroid_path <- function(object,
                                         model_params,
                                         continuation_steps = 6L) {
  # Skip non-spheroidal continuation targets ===================================
  shape_parameters <- acousticTS::extract(object, "shape_parameters")
  if (!(shape_parameters$shape %in% c("ProlateSpheroid", "OblateSpheroid"))) {
    return(NULL)
  }

  target_aspect_ratio <- .tmm_spheroidal_aspect_ratio(shape_parameters)
  radius_eq <- .tmm_equivalent_volume_radius(shape_parameters)
  if (!is.finite(target_aspect_ratio) || !is.finite(radius_eq)) {
    return(NULL)
  }
  # Validate the continuation path length ======================================
  if (!is.numeric(continuation_steps) || length(continuation_steps) != 1 ||
      !is.finite(continuation_steps) || continuation_steps < 0 ||
      continuation_steps %% 1 != 0) {
    stop("'continuation_steps' must be a single integer >= 0.", call. = FALSE)
  }
  if (continuation_steps < 2) {
    return(NULL)
  }
  # Recover the baseline object and medium state ===============================
  body <- acousticTS::extract(object, "body")
  medium <- model_params$medium
  acoustics <- model_params$parameters$acoustics
  boundary <- model_params$parameters$boundary
  frequency <- acoustics$frequency
  n_segments <- shape_parameters$n_segments %||% 80L

  aspect_path <- unique(c(
    1,
    seq(
      from = 1 + (target_aspect_ratio - 1) / max(continuation_steps - 1, 1),
      to = target_aspect_ratio,
      length.out = max(continuation_steps - 1, 1)
    )
  ))
  if (target_aspect_ratio <= 1 + sqrt(.Machine$double.eps)) {
    aspect_path <- 1
  }
  # Rebuild and solve each continuation step ===================================
  path_rows <- lapply(
    seq_along(aspect_path),
    function(i) {
      aspect_i <- aspect_path[i]
      if (aspect_i <= 1 + sqrt(.Machine$double.eps)) {
        shape_i <- sphere(radius_body = radius_eq, n_segments = n_segments)
        shape_label <- "Sphere"
      } else {
        axes_i <- .tmm_spheroid_axes_from_equal_volume(
          shape_type = shape_parameters$shape,
          radius_eq = radius_eq,
          aspect_ratio = aspect_i
        )
        shape_i <- if (identical(shape_parameters$shape, "ProlateSpheroid")) {
          prolate_spheroid(
            length_body = axes_i$length_body,
            radius_body = axes_i$radius_body,
            n_segments = n_segments
          )
        } else {
          oblate_spheroid(
            length_body = axes_i$length_body,
            radius_body = axes_i$radius_body,
            n_segments = n_segments
          )
        }
        shape_label <- shape_parameters$shape
      }

      object_i <- .tmm_rebuild_shape_like(object, shape = shape_i, body = body)
      object_i <- target_strength(
        object = object_i,
        frequency = frequency,
        model = "tmm",
        boundary = boundary,
        density_sw = medium$density,
        sound_speed_sw = medium$sound_speed
      )

      data.frame(
        step = i,
        shape = shape_label,
        aspect_ratio = rep(aspect_i, length(frequency)),
        target_aspect_ratio = rep(target_aspect_ratio, length(frequency)),
        frequency = frequency,
        sigma_bs = object_i@model$TMM$sigma_bs,
        TS = object_i@model$TMM$TS
      )
    }
  )

  do.call(rbind, path_rows)
}

# Summarize the smoothness of one sphere-to-spheroid continuation path.
#' @noRd
.tmm_sphere_to_spheroid_summary <- function(path_df) {
  # Skip missing continuation paths ============================================
  if (is.null(path_df) || !nrow(path_df)) {
    return(NULL)
  }
  # Summarize first- and second-difference smoothness ==========================
  split_path <- split(path_df, path_df$frequency)
  out <- lapply(
    split_path,
    function(df) {
      df <- df[order(df$aspect_ratio), ]
      ts_diff <- diff(df$TS)
      ts_second_diff <- diff(ts_diff)
      data.frame(
        frequency = df$frequency[1],
        continuation_target_aspect_ratio = df$target_aspect_ratio[1],
        continuation_max_abs_step_TS = if (length(ts_diff)) max(abs(ts_diff)) else 0,
        continuation_max_abs_second_diff_TS = if (length(ts_second_diff)) max(abs(ts_second_diff)) else 0,
        continuation_any_nonfinite = any(!is.finite(df$TS))
      )
    }
  )

  do.call(rbind, out)
}

# Integrate one scattering grid over solid angle for the optical-theorem check.
#' @noRd
.tmm_total_scattering_cross_section <- function(theta_scatter,
                                                phi_scatter,
                                                sigma_scat) {
  # Build the quadrature-like interval weights =================================
  w_theta <- .tmm_interval_weights(
    theta_scatter,
    lower = min(theta_scatter),
    upper = max(theta_scatter),
    name = "theta_scatter"
  )
  w_phi <- .tmm_interval_weights(
    phi_scatter,
    lower = min(phi_scatter),
    upper = max(phi_scatter),
    name = "phi_scatter"
  )
  # Integrate the differential cross section over solid angle ==================
  sum(outer(w_theta * sin(theta_scatter), w_phi) * sigma_scat)
}

# Summarize one stored frequency's block-level conditioning indicators.
#' @noRd
.tmm_block_metrics <- function(t_blocks) {
  # Initialize the empty output template =======================================
  empty <- data.frame(
    m = numeric(0),
    n_terms = numeric(0),
    rcond = numeric(0),
    transpose_residual = numeric(0)
  )

  if (!is.list(t_blocks) || !length(t_blocks)) {
    return(empty)
  }

  valid_blocks <- Filter(
    function(block) is.list(block) && !is.null(block$T),
    t_blocks
  )
  if (!length(valid_blocks)) {
    return(empty)
  }
  # Summarize each retained T-matrix block =====================================
  out <- lapply(
    valid_blocks,
    function(block) {
      T_block <- block$T
      den <- sqrt(sum(Mod(T_block)^2))
      transpose_residual <- if (!is.finite(den) || den == 0) {
        0
      } else {
        sqrt(sum(Mod(T_block - t(T_block))^2)) / den
      }

      rcond_value <- block$rcond_lhs %||% tryCatch(
        rcond(T_block),
        error = function(...) NA_real_
      )

      data.frame(
        m = block$m,
        n_terms = length(block$n_seq),
        rcond = rcond_value,
        transpose_residual = transpose_residual
      )
    }
  )

  do.call(rbind, out)
}

# Compute reciprocity residuals by swapping the incident and receive angles at
# fixed frequency.
#' @noRd
.tmm_reciprocity_residual <- function(model_params,
                                      frequency_idx,
                                      shape_parameters,
                                      reciprocity_pairs) {
  # Evaluate scattering in the forward ordering ================================
  f_forward <- .tmm_scattering_points(
    model_params = model_params,
    frequency_idx = frequency_idx,
    shape_parameters = shape_parameters,
    theta_body = reciprocity_pairs$theta_body,
    phi_body = reciprocity_pairs$phi_body,
    theta_scatter = reciprocity_pairs$theta_scatter,
    phi_scatter = reciprocity_pairs$phi_scatter
  )
  f_reverse <- .tmm_scattering_points(
    model_params = model_params,
    frequency_idx = frequency_idx,
    shape_parameters = shape_parameters,
    theta_body = reciprocity_pairs$theta_scatter,
    phi_body = reciprocity_pairs$phi_scatter,
    theta_scatter = reciprocity_pairs$theta_body,
    phi_scatter = reciprocity_pairs$phi_body
  )
  # Compare the forward and reversed angle orderings ===========================
  denom <- pmax(Mod(f_forward), Mod(f_reverse), .Machine$double.eps)
  max(Mod(f_forward - f_reverse) / denom)
}

#' Diagnostics for stored single-target TMM solutions
#'
#' @description
#' Reuses the stored T-matrix blocks to compute a compact set of numerical and
#' physics-based diagnostics for one or more stored frequencies. The summary
#' combines:
#'
#' \itemize{
#'   \item monostatic reconstruction residuals,
#'   \item reciprocity residuals under incident/receive-angle exchange,
#'   \item an optical-theorem residual based on forward scattering and the
#'   integrated differential cross section,
#'   \item block-level conditioning indicators from the stored T-matrix blocks,
#'   \item and, for prolate/oblate targets, an equal-volume sphere-to-spheroid
#'   continuation path that checks whether the monostatic response deforms
#'   smoothly away from the exact sphere limit.
#' }
#'
#' These checks are meant to help distinguish "the post-processing is
#' self-consistent" from "the retained solve was also numerically comfortable,"
#' which is especially helpful for the newer nonspherical TMM branches.
#' For stored cylinders, the current diagnostics are intentionally limited to
#' exact monostatic reconstruction because a validated retained cylinder angular
#' operator is not yet available.
#'
#' @param object Scatterer-object previously evaluated with
#'   `target_strength(..., model = "TMM", store_t_matrix = TRUE)`.
#' @param frequency Optional stored frequency or vector of stored frequencies
#'   (Hz). Defaults to all stored frequencies.
#' @param theta_body Incident polar angle (radians) used for the monostatic and
#'   optical-theorem checks. Defaults to the stored TMM incident angle.
#' @param phi_body Incident azimuth angle (radians) used for the monostatic and
#'   optical-theorem checks. Defaults to the stored TMM incident angle.
#' @param reciprocity_pairs Optional data frame giving explicit reciprocity test
#'   angles. Must contain `theta_body`, `phi_body`, `theta_scatter`, and
#'   `phi_scatter` columns in radians.
#' @param n_theta Number of polar-angle grid points used by the
#'   optical-theorem integration check.
#' @param n_phi Number of azimuth-angle grid points used by the
#'   optical-theorem integration check.
#' @param continuation_steps Number of equal-volume aspect-ratio steps used for
#'   the sphere-to-spheroid continuation check on prolate and oblate targets.
#'   Set to `0` or `1` to skip the continuation path. Ignored for non-spheroidal
#'   targets.
#'
#' @return A list with components:
#' \describe{
#'   \item{\code{summary}}{Per-frequency diagnostic summary.}
#'   \item{\code{block_metrics}}{Per-frequency block-level conditioning and
#'   transpose-residual summaries.}
#'   \item{\code{continuation}}{Equal-volume sphere-to-spheroid continuation
#'   path for spheroidal targets, or `NULL` for other shapes.}
#' }
#'
#' @seealso \code{\link{tmm_scattering}}, \code{\link{tmm_scattering_grid}},
#'   \code{\link{tmm_products}}
#' @export
tmm_diagnostics <- function(object,
                            frequency = NULL,
                            theta_body = NULL,
                            phi_body = NULL,
                            reciprocity_pairs = NULL,
                            n_theta = 61,
                            n_phi = 121,
                            continuation_steps = 6L) {
  # Recover the stored TMM state ===============================================
  model_params <- .tmm_require_stored_blocks(object)
  .tmm_warn_exploratory_cylinder_blocks(object, model_params)
  parameters <- model_params$parameters
  acoustics <- parameters$acoustics
  defaults <- model_params$body
  shape_parameters <- acousticTS::extract(object, "shape_parameters")
  shape_type <- shape_parameters$shape
  # Validate the grid controls =================================================
  if (!is.numeric(n_theta) || length(n_theta) != 1 || !is.finite(n_theta) ||
      n_theta < 5 || n_theta %% 1 != 0) {
    stop("'n_theta' must be a single integer >= 5.", call. = FALSE)
  }
  if (!is.numeric(n_phi) || length(n_phi) != 1 || !is.finite(n_phi) ||
      n_phi < 5 || n_phi %% 1 != 0) {
    stop("'n_phi' must be a single integer >= 5.", call. = FALSE)
  }
  # Resolve the active angles, frequencies, and continuation path ==============
  theta_body <- .tmm_scalar_angle(theta_body, defaults$theta_body, "theta_body")
  phi_body <- .tmm_scalar_angle(phi_body, defaults$phi_body %||% pi, "phi_body")
  reciprocity_pairs <- .tmm_validate_reciprocity_pairs(reciprocity_pairs)
  freq_idx <- .tmm_frequency_indices(frequency, acoustics$frequency)
  continuation <- .tmm_sphere_to_spheroid_path(
    object = object,
    model_params = model_params,
    continuation_steps = continuation_steps
  )
  continuation_summary <- .tmm_sphere_to_spheroid_summary(continuation)
  # Initialize the per-frequency outputs =======================================
  summary_rows <- vector("list", length(freq_idx))
  block_metrics <- vector("list", length(freq_idx))
  # Evaluate the diagnostics for each stored frequency =========================
  for (i in seq_along(freq_idx)) {
    idx <- freq_idx[i]
    freq_val <- acoustics$frequency[idx]
    t_blocks <- parameters$t_matrix[[idx]]
    block_df <- .tmm_block_metrics(t_blocks)
    block_metrics[[i]] <- block_df

    f_mono <- .tmm_scattering_points(
      model_params = model_params,
      frequency_idx = idx,
      shape_parameters = shape_parameters,
      theta_body = theta_body,
      phi_body = phi_body,
      theta_scatter = pi - theta_body,
      phi_scatter = (phi_body + pi) %% (2 * pi)
    )[1]
    stored_f_bs <- acousticTS::extract(object, "model")$TMM$f_bs[idx]
    mono_abs <- Mod(f_mono - stored_f_bs)
    mono_rel <- mono_abs / max(Mod(stored_f_bs), .Machine$double.eps)

    reciprocity_rel <- if (identical(parameters$coordinate_system, "cylindrical")) {
      NA_real_
    } else {
      .tmm_reciprocity_residual(
        model_params = model_params,
        frequency_idx = idx,
        shape_parameters = shape_parameters,
        reciprocity_pairs = reciprocity_pairs
      )
    }

    sigma_total <- NA_real_
    sigma_ext <- NA_real_
    optical_sign <- NA_real_
    optical_rel <- NA_real_
    if (identical(parameters$coordinate_system, "spherical")) {
      theta_grid <- seq(0, pi, length.out = n_theta)
      phi_grid <- seq(0, 2 * pi, length.out = n_phi)
      grid <- .tmm_scattering_grid_matrix(
        model_params = model_params,
        frequency_idx = idx,
        shape_parameters = shape_parameters,
        theta_body = theta_body,
        phi_body = phi_body,
        theta_scatter = theta_grid,
        phi_scatter = phi_grid
      )
      sigma_grid <- .sigma_bs(grid)
      sigma_total <- .tmm_total_scattering_cross_section(
        theta_scatter = theta_grid,
        phi_scatter = phi_grid,
        sigma_scat = sigma_grid
      )
      f_forward <- .tmm_scattering_points(
        model_params = model_params,
        frequency_idx = idx,
        shape_parameters = shape_parameters,
        theta_body = theta_body,
        phi_body = phi_body,
        theta_scatter = theta_body,
        phi_scatter = phi_body
      )[1]
      sigma_ext_raw <- (4 * pi / acoustics$k_sw[idx]) * Im(f_forward)
      optical_sign <- if (
        abs(sigma_total - sigma_ext_raw) <= abs(sigma_total + sigma_ext_raw)
      ) 1 else -1
      sigma_ext <- optical_sign * sigma_ext_raw
      optical_rel <- abs(sigma_total - sigma_ext) /
        max(abs(sigma_total), abs(sigma_ext), .Machine$double.eps)
    }

    min_rcond <- if (nrow(block_df)) min(block_df$rcond, na.rm = TRUE) else NA_real_
    if (!is.finite(min_rcond)) {
      min_rcond <- NA_real_
    }
    max_cond <- if (is.na(min_rcond) || min_rcond <= 0) Inf else 1 / min_rcond
    max_transpose <- if (nrow(block_df)) max(block_df$transpose_residual, na.rm = TRUE) else NA_real_
    continuation_row <- if (!is.null(continuation_summary)) {
      continuation_summary[continuation_summary$frequency == freq_val, , drop = FALSE]
    } else {
      NULL
    }
    # Store the scalar diagnostics for this frequency ==========================
    summary_rows[[i]] <- data.frame(
      shape = shape_type,
      coordinate_system = parameters$coordinate_system,
      boundary = parameters$boundary,
      frequency = freq_val,
      monostatic_abs_residual = mono_abs,
      monostatic_rel_residual = mono_rel,
      reciprocity_rel_residual = reciprocity_rel,
      sigma_total = sigma_total,
      sigma_ext = sigma_ext,
      optical_theorem_sign = optical_sign,
      optical_theorem_rel_residual = optical_rel,
      min_block_rcond = min_rcond,
      max_block_cond_est = max_cond,
      max_block_transpose_residual = max_transpose,
      continuation_target_aspect_ratio = if (is.null(continuation_row) || !nrow(continuation_row)) NA_real_ else continuation_row$continuation_target_aspect_ratio[1],
      continuation_max_abs_step_TS = if (is.null(continuation_row) || !nrow(continuation_row)) NA_real_ else continuation_row$continuation_max_abs_step_TS[1],
      continuation_max_abs_second_diff_TS = if (is.null(continuation_row) || !nrow(continuation_row)) NA_real_ else continuation_row$continuation_max_abs_second_diff_TS[1],
      continuation_any_nonfinite = if (is.null(continuation_row) || !nrow(continuation_row)) NA else continuation_row$continuation_any_nonfinite[1]
    )
  }
  # Name the block metrics and return the full diagnostic bundle ===============
  names(block_metrics) <- as.character(acoustics$frequency[freq_idx])

  list(
    summary = do.call(rbind, summary_rows),
    block_metrics = block_metrics,
    continuation = continuation
  )
}
