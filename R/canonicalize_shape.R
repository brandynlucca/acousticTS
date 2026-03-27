#' Canonicalize one shape into a canonical surrogate
#'
#' @description
#' Build an explicit canonical surrogate shape from an existing `Shape` object.
#' This is intended for workflows where a measured or segmented body should be
#' approximated deliberately by one of the package's canonical geometries before
#' being passed to a shape-specific model family such as `SPHMS`, `PSMS`, or
#' `FCMS`.
#'
#' @param shape A `Shape` object to approximate.
#' @param to Canonical target shape. One of `"Sphere"`, `"Cylinder"`,
#'   `"ProlateSpheroid"`, or `"OblateSpheroid"`.
#' @param method Canonicalization rule. `"auto"` selects a shape-specific
#'   default: `"volume"` for spheres, `"length_volume"` for cylinders, and
#'   `"profile_l2"` for spheroids. `"volume"` preserves enclosed volume only.
#'   `"length_volume"` preserves body length and enclosed volume.
#'   `"profile_l2"` fits the target canonical profile directly to the source
#'   equivalent-radius profile by least squares.
#' @param n_segments Optional number of segments for the returned canonical
#'   shape. Defaults to the source shape segment count.
#' @param diagnostics Logical. If `FALSE` (default), return only the canonical
#'   `Shape`. If `TRUE`, return a list containing the canonical shape and fit
#'   diagnostics.
#'
#' @details
#' `canonicalize_shape()` is intentionally explicit. It does **not** run
#' automatically inside model calls such as `target_strength(..., model =
#' "psms")`, because there is no single defensible canonicalization rule for
#' every biological or segmented target.
#'
#' When the input shape contains both an explicit width profile and upper/lower
#' height profiles, the canonicalization step first reduces those local
#' cross-sections to an equal-area circular radius before fitting the canonical
#' surrogate. That keeps the reduction axisymmetric while still honoring the
#' original local cross-sectional area more closely than a height-only radius.
#'
#' The returned diagnostics report how the source and target compare in length,
#' enclosed volume, maximum equivalent radius, and equivalent-radius profile
#' root-mean-square error on the source x-grid. Those diagnostics are the
#' package's recommended way to judge whether a canonical surrogate is
#' defensible for the intended model comparison.
#'
#' @return
#' If `diagnostics = FALSE`, a canonical `Shape` object. If
#' `diagnostics = TRUE`, a list with elements:
#' \itemize{
#'   \item `shape`: the fitted canonical `Shape`
#'   \item `diagnostics`: a named list of source, target, and fit metrics
#' }
#'
#' @examples
#' bladder_like <- arbitrary(
#'   x_body = c(0, 0.01, 0.02, 0.03, 0.04),
#'   radius_body = c(0, 0.004, 0.006, 0.004, 0)
#' )
#'
#' psms_shape <- canonicalize_shape(
#'   bladder_like,
#'   to = "ProlateSpheroid"
#' )
#'
#' cyl_fit <- canonicalize_shape(
#'   bladder_like,
#'   to = "Cylinder",
#'   diagnostics = TRUE
#' )
#' cyl_fit$diagnostics$fit$radius_nrmse
#'
#' @seealso [create_shape()], [arbitrary()], [sphere()], [cylinder()],
#'   [prolate_spheroid()], [oblate_spheroid()]
#'
#' @keywords shape_generation
#' @export
canonicalize_shape <- function(shape,
                               to = c(
                                 "ProlateSpheroid",
                                 "Cylinder",
                                 "Sphere",
                                 "OblateSpheroid"
                               ),
                               method = c(
                                 "auto",
                                 "volume",
                                 "length_volume",
                                 "profile_l2"
                               ),
                               n_segments = NULL,
                               diagnostics = FALSE) {
  # Validate source shape ======================================================
  if (!methods::is(shape, "Shape")) {
    stop("'shape' must be a Shape object.", call. = FALSE)
  }
  # Match the canonical target and fitting rule ================================
  to <- match.arg(to)
  method <- match.arg(method)
  # Extract source shape fields ================================================
  position_matrix <- acousticTS::extract(shape, "position_matrix")
  shape_parameters <- acousticTS::extract(shape, "shape_parameters")
  # Resolve the output segmentation ============================================
  if (is.null(n_segments)) {
    n_segments <- shape_parameters$n_segments %||%
      .shape_segment_count(position_matrix = position_matrix)
  }
  n_segments <- as.integer(n_segments)
  if (!is.numeric(n_segments) || length(n_segments) != 1L ||
      is.na(n_segments) || n_segments < 2L) {
    stop(
      "'n_segments' must be one integer greater than or equal to 2.",
      call. = FALSE
    )
  }
  # Resolve the effective fitting method =======================================
  resolved_method <- .canonicalize_shape_method(method, to)
  source_metrics <- .canonicalize_shape_metrics(shape)
  # Validate the derived source metrics ========================================
  if (!is.finite(source_metrics$length) || source_metrics$length <= 0) {
    stop("Source shape must have positive length.", call. = FALSE)
  }
  if (!is.finite(source_metrics$volume) || source_metrics$volume <= 0) {
    stop("Source shape must have positive enclosed volume.", call. = FALSE)
  }
  if (!is.finite(source_metrics$max_radius) || source_metrics$max_radius <= 0) {
    stop(
      "Source shape must have positive radius somewhere along the profile.",
      call. = FALSE
    )
  }
  # Fit the requested canonical surrogate ======================================
  fit <- switch(
    resolved_method,
    volume = .canonicalize_shape_fit_volume(source_metrics, to),
    length_volume = .canonicalize_shape_fit_length_volume(source_metrics, to),
    profile_l2 = .canonicalize_shape_fit_profile_l2(source_metrics, to)
  )
  # Build the output shape and diagnostics =====================================
  out_shape <- .canonicalize_shape_build(
    to = to,
    fit = fit,
    n_segments = n_segments,
    source_parameters = shape_parameters
  )

  diag <- .canonicalize_shape_diagnostics(
    source_metrics = source_metrics,
    target_shape = out_shape,
    to = to,
    method = resolved_method,
    fit = fit
  )
  # Return the requested output format =========================================
  if (isTRUE(diagnostics)) {
    return(list(shape = out_shape, diagnostics = diag))
  }

  out_shape
}

#' Resolve the effective canonicalization rule for one target family
#' @keywords internal
#' @noRd
.canonicalize_shape_method <- function(method, to) {
  # Respect explicit method choices ============================================
  if (!identical(method, "auto")) {
    return(method)
  }
  # Assign the default rule for the target family ==============================
  switch(
    to,
    Sphere = "volume",
    Cylinder = "length_volume",
    ProlateSpheroid = "profile_l2",
    OblateSpheroid = "profile_l2"
  )
}

#' Estimate an equivalent axisymmetric radius profile from a shape matrix
#' @keywords internal
#' @noRd
.shape_equivalent_radius_profile <- function(position_matrix) {
  # Validate the input geometry matrix =========================================
  if (!is.matrix(position_matrix)) {
    stop("'position_matrix' must be a matrix.", call. = FALSE)
  }
  # Resolve the available axis names ===========================================
  axis_names <- colnames(position_matrix)
  if (is.null(axis_names)) {
    axis_names <- character()
  }
  # Build a helper for the first matching axis =================================
  get_first <- function(candidates) {
    idx <- match(candidates, axis_names, nomatch = 0L)
    idx <- idx[idx > 0L]
    if (length(idx) == 0L) {
      return(NULL)
    }
    as.numeric(position_matrix[, idx[1]])
  }
  # Pull candidate radius, width, and height fields ============================
  explicit_radius <- get_first(c(
    "a", "radius", "radius_body", "radius_bladder", "radius_shell",
    "radius_fluid", "radius_backbone"
  ))
  explicit_width <- get_first(c(
    "w", "w_body", "w_bladder", "w_shell", "w_fluid", "w_backbone"
  ))
  zU <- get_first(c(
    "zU", "zU_body", "zU_bladder", "zU_shell", "zU_fluid", "zU_backbone"
  ))
  zL <- get_first(c(
    "zL", "zL_body", "zL_bladder", "zL_shell", "zL_fluid", "zL_backbone"
  ))
  # Prefer the most explicit available radius representation ===================
  if (!is.null(explicit_radius) && !all(is.na(explicit_radius))) {
    return(abs(explicit_radius))
  }

  if (!is.null(explicit_width) && !is.null(zU) && !is.null(zL)) {
    semi_width <- abs(explicit_width) / 2
    semi_height <- abs(zU - zL) / 2
    return(sqrt(semi_width * semi_height))
  }

  if (!is.null(zU) && !is.null(zL)) {
    return(abs(zU - zL) / 2)
  }

  if (!is.null(explicit_width)) {
    return(abs(explicit_width) / 2)
  }
  # Fall back to the package-standard radius extraction ========================
  .shape_radius_profile(position_matrix = position_matrix)
}

#' Integrate an equivalent-volume profile from x nodes and radius nodes
#' @keywords internal
#' @noRd
.shape_equivalent_volume <- function(x, radius_eq) {
  # Short-circuit degenerate profiles ==========================================
  if (length(x) < 2L || length(radius_eq) < 2L) {
    return(0)
  }
  # Integrate the equivalent circular area profile =============================
  area <- pi * radius_eq^2
  dx <- diff(x)
  sum(dx * (utils::head(area, -1L) + utils::tail(area, -1L)) / 2)
}

#' Collect shape metrics used by the canonicalization helpers
#' @keywords internal
#' @noRd
.canonicalize_shape_metrics <- function(shape) {
  # Extract and sort the source profile ========================================
  position_matrix <- acousticTS::extract(shape, "position_matrix")
  x <- .shape_x(position_matrix = position_matrix)
  ord <- order(x)
  x <- x[ord]
  radius_eq <- .shape_equivalent_radius_profile(position_matrix)[ord]
  x <- x - min(x, na.rm = TRUE)
  # Summarize the source geometry on the canonical axis ========================
  list(
    class = class(shape)[1],
    x = as.numeric(x),
    radius_eq = as.numeric(radius_eq),
    length = .shape_length(position_matrix = position_matrix),
    volume = .shape_equivalent_volume(as.numeric(x), as.numeric(radius_eq)),
    max_radius = max(radius_eq, na.rm = TRUE),
    n_segments = .shape_segment_count(position_matrix = position_matrix)
  )
}

#' Canonical radius profile evaluated on one source x-grid
#' @keywords internal
#' @noRd
.canonical_radius_profile <- function(x, to, length, radius) {
  # Normalize the evaluation grid ==============================================
  x <- as.numeric(x)
  # Return zeros for invalid canonical dimensions ==============================
  if (!is.finite(length) || length <= 0 || !is.finite(radius) || radius <= 0) {
    return(rep(0, base::length(x)))
  }
  # Evaluate the target canonical profile ======================================
  switch(
    to,
    Sphere = sqrt(pmax(radius^2 - (x - radius)^2, 0)),
    Cylinder = ifelse(x >= 0 & x <= length, radius, 0),
    ProlateSpheroid = {
      half_length <- length / 2
      radius * sqrt(pmax(1 - ((x - half_length) / half_length)^2, 0))
    },
    OblateSpheroid = {
      half_length <- length / 2
      radius * sqrt(pmax(1 - ((x - half_length) / half_length)^2, 0))
    }
  )
}

#' Fit a pure equal-volume surrogate
#' @keywords internal
#' @noRd
.canonicalize_shape_fit_volume <- function(source_metrics, to) {
  # Restrict pure volume fitting to the sphere case ============================
  if (!identical(to, "Sphere")) {
    stop(
      "'method = \"volume\"' is currently only supported for 'to = ",
      "\"Sphere\"'.",
      call. = FALSE
    )
  }
  # Convert the source volume into an equal-volume sphere ======================
  radius <- (3 * source_metrics$volume / (4 * pi))^(1 / 3)
  list(length = 2 * radius, radius = radius, objective = 0)
}

#' Fit a length-plus-volume surrogate
#' @keywords internal
#' @noRd
.canonicalize_shape_fit_length_volume <- function(source_metrics, to) {
  # Preserve the source length and volume ======================================
  length_out <- source_metrics$length
  volume_out <- source_metrics$volume
  # Reuse the sphere-specific equal-volume rule ================================
  if (identical(to, "Sphere")) {
    return(.canonicalize_shape_fit_volume(source_metrics, to))
  }
  # Solve the equal-length equal-volume cylinder ===============================
  if (identical(to, "Cylinder")) {
    radius_out <- sqrt(volume_out / (pi * length_out))
    return(list(length = length_out, radius = radius_out, objective = 0))
  }
  # Solve the equal-length equal-volume prolate spheroid =======================
  if (identical(to, "ProlateSpheroid")) {
    semimajor <- length_out / 2
    semiminor <- sqrt((3 * volume_out) / (4 * pi * semimajor))
    if (semiminor > semimajor) {
      stop(
        "The source shape cannot preserve both length and volume while ",
        "remaining a prolate spheroid. Try 'method = \"profile_l2\"' or ",
        "canonicalize to 'OblateSpheroid' instead.",
        call. = FALSE
      )
    }
    return(list(length = length_out, radius = semiminor, objective = 0))
  }
  # Solve the equal-length equal-volume oblate spheroid ========================
  if (identical(to, "OblateSpheroid")) {
    semiminor <- length_out / 2
    semimajor <- sqrt((3 * volume_out) / (4 * pi * semiminor))
    if (semimajor < semiminor) {
      stop(
        "The source shape cannot preserve both length and volume while ",
        "remaining an oblate spheroid. Try 'method = \"profile_l2\"' or ",
        "canonicalize to 'ProlateSpheroid' instead.",
        call. = FALSE
      )
    }
    return(list(length = length_out, radius = semimajor, objective = 0))
  }

  stop("Unsupported canonical target.", call. = FALSE)
}

#' Fit a least-squares canonical surrogate to the source radius profile
#' @keywords internal
#' @noRd
.canonicalize_shape_fit_profile_l2 <- function(source_metrics, to) {
  # Unpack the source profile and initial scales ===============================
  x <- source_metrics$x
  r <- source_metrics$radius_eq
  L0 <- source_metrics$length
  a0 <- source_metrics$max_radius
  eps <- max(1e-9, L0 * 1e-6, a0 * 1e-6)
  # Optimize the equal-profile sphere radius ===================================
  if (identical(to, "Sphere")) {
    lower <- eps
    upper <- max(L0 * 2, a0 * 5, eps * 10)
    sphere_objective <- function(radius) {
      model <- .canonical_radius_profile(x, "Sphere", 2 * radius, radius)
      mean((r - model)^2)
    }
    opt <- stats::optimize(sphere_objective, interval = c(lower, upper))
    return(list(length = 2 * opt$minimum,
                radius = opt$minimum,
                objective = opt$objective))
  }
  # Shift the starting guess toward the admissible aspect-ratio side ===========
  if (identical(to, "ProlateSpheroid")) {
    a0 <- min(a0, max(L0 / 2 * 0.95, eps))
    if (L0 <= 2 * a0) {
      L0 <- 2 * a0 * 1.05
    }
  }

  if (identical(to, "OblateSpheroid")) {
    a0 <- max(a0, L0 / 2 * 1.05)
  }
  # Define the penalized least-squares objective ===============================
  penalty_scale <- 1e6 / max(a0^2, eps)
  profile_objective <- function(par_log) {
    length_fit <- exp(par_log[1])
    radius_fit <- exp(par_log[2])

    penalty <- 0
    if (identical(to, "ProlateSpheroid") && length_fit < 2 * radius_fit) {
      penalty <- penalty_scale * (2 * radius_fit - length_fit)^2
    }
    if (identical(to, "OblateSpheroid") && length_fit > 2 * radius_fit) {
      penalty <- penalty_scale * (length_fit - 2 * radius_fit)^2
    }

    model <- .canonical_radius_profile(x, to, length_fit, radius_fit)
    mean((r - model)^2) + penalty
  }
  # Set bounds and optimize on the log scale ===================================
  lower <- log(c(
    max(L0 / 10, eps),
    max(a0 / 10, eps)
  ))
  upper <- log(c(
    max(L0 * 10, a0 * 10, eps * 10),
    max(a0 * 10, L0 * 2, eps * 10)
  ))

  opt <- stats::optim(
    par = log(c(L0, a0)),
    fn = profile_objective,
    method = "L-BFGS-B",
    lower = lower,
    upper = upper
  )
  # Return the fitted canonical dimensions =====================================
  list(
    length = exp(opt$par[1]),
    radius = exp(opt$par[2]),
    objective = opt$value
  )
}

#' Build one canonical Shape from fitted length/radius parameters
#' @keywords internal
#' @noRd
.canonicalize_shape_build <- function(to,
                                      fit,
                                      n_segments,
                                      source_parameters) {
  # Recover the source length-unit metadata ====================================
  length_units <- source_parameters$length_units %||%
    source_parameters$diameter_units %||%
    "m"
  diameter_units <- source_parameters$diameter_units %||%
    source_parameters$length_units %||%
    "m"
  # Build the requested canonical surrogate ====================================
  switch(
    to,
    Sphere = sphere(
      radius_body = fit$radius,
      n_segments = n_segments,
      diameter_units = diameter_units
    ),
    Cylinder = cylinder(
      length_body = fit$length,
      radius_body = fit$radius,
      n_segments = n_segments,
      length_units = length_units
    ),
    ProlateSpheroid = prolate_spheroid(
      length_body = fit$length,
      radius_body = fit$radius,
      n_segments = n_segments,
      length_units = length_units
    ),
    OblateSpheroid = oblate_spheroid(
      length_body = fit$length,
      radius_body = fit$radius,
      n_segments = n_segments,
      length_units = length_units
    )
  )
}

#' Build canonicalization diagnostics from source/target metrics
#' @keywords internal
#' @noRd
.canonicalize_shape_diagnostics <- function(source_metrics,
                                            target_shape,
                                            to,
                                            method,
                                            fit) {
  # Re-evaluate the fitted target on the source x-grid =========================
  target_metrics <- .canonicalize_shape_metrics(target_shape)
  model_on_source_grid <- .canonical_radius_profile(
    source_metrics$x,
    to = to,
    length = target_metrics$length,
    radius = target_metrics$max_radius
  )
  rmse <- sqrt(mean((source_metrics$radius_eq - model_on_source_grid)^2))
  # Return source, target, and fit summaries ===================================
  list(
    source = list(
      shape = source_metrics$class,
      length = source_metrics$length,
      volume = source_metrics$volume,
      max_radius = source_metrics$max_radius,
      length_radius_ratio = source_metrics$length / source_metrics$max_radius
    ),
    target = list(
      shape = to,
      length = target_metrics$length,
      volume = target_metrics$volume,
      max_radius = target_metrics$max_radius,
      length_radius_ratio = target_metrics$length / target_metrics$max_radius
    ),
    fit = list(
      method = method,
      objective = fit$objective,
      radius_rmse = rmse,
      radius_nrmse = rmse / source_metrics$max_radius,
      length_ratio = target_metrics$length / source_metrics$length,
      volume_ratio = target_metrics$volume / source_metrics$volume,
      max_radius_ratio = target_metrics$max_radius / source_metrics$max_radius
    )
  )
}
