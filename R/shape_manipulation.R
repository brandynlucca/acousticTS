################################################################################
################################################################################
# SHAPE-MANIPULATION HELPERS
################################################################################
################################################################################
#' Detect the internal geometry storage format for a position matrix
#' @param position_matrix Numeric position matrix.
#' @keywords internal
#' @noRd
.geometry_storage <- function(position_matrix) {
  # Validate the candidate storage contracts ==================================
  row_major_ok <- tryCatch(
    {
      .validate_geometry_contract(
        position_matrix,
        storage = "profile_row_major",
        context = "Geometry"
      )
      TRUE
    },
    error = function(e) FALSE
  )
  col_major_ok <- tryCatch(
    {
      .validate_geometry_contract(
        position_matrix,
        storage = "shape_column_major",
        context = "Geometry"
      )
      TRUE
    },
    error = function(e) FALSE
  )

  # Return the resolved geometry-storage label ================================
  if (row_major_ok && !col_major_ok) {
    return("profile_row_major")
  }
  if (col_major_ok && !row_major_ok) {
    return("shape_column_major")
  }
  if (row_major_ok && col_major_ok) {
    if (!is.null(rownames(position_matrix)) && nrow(position_matrix) <= 6) {
      return("profile_row_major")
    }
    return("shape_column_major")
  }

  stop(
    "Unable to determine whether the geometry is stored as a shape matrix ",
    "or a row-major component profile.",
    call. = FALSE
  )
}

#' Resolve the default manipulator component for a scatterer object
#' @param object Scatterer-class object.
#' @keywords internal
#' @noRd
.default_shape_component <- function(object) {
  # Preserve shell as the primary ESS geometry ================================
  if (methods::is(object, "ESS")) {
    return("shell")
  }

  # Return body for all other supported scatterers ============================
  "body"
}

#' Resolve matrix metadata for shape-manipulation helpers
#' @param object Shape or Scatterer object.
#' @param component Optional scatterer component name.
#' @keywords internal
#' @noRd
.shape_manipulation_info <- function(object, component = NULL) {
  # Resolve the direct Shape pathway ==========================================
  if (methods::is(object, "Shape")) {
    position_matrix <- acousticTS::extract(object, "position_matrix")
    return(list(
      object_type = "Shape",
      component = NULL,
      position_matrix = position_matrix,
      storage = .geometry_storage(position_matrix)
    ))
  }

  # Validate the Scatterer pathway ============================================
  if (!methods::is(object, "Scatterer")) {
    stop("Expected either a Shape or Scatterer object.", call. = FALSE)
  }

  component <- component %||% .default_shape_component(object)
  if (!component %in% methods::slotNames(object)) {
    stop(
      "Scatterer does not contain a '", component, "' slot.",
      call. = FALSE
    )
  }

  component_value <- methods::slot(object, component)
  if (!is.list(component_value) ||
    is.null(component_value$rpos) ||
    !is.matrix(component_value$rpos)) {
    stop(
      "Scatterer component '", component, "' does not contain a valid ",
      "position matrix in 'rpos'.",
      call. = FALSE
    )
  }

  # Return scatterer-component manipulation metadata ==========================
  list(
    object_type = "Scatterer",
    component = component,
    position_matrix = component_value$rpos,
    storage = .geometry_storage(component_value$rpos)
  )
}

#' Resolve summary information for a manipulated geometry matrix
#' @param position_matrix Numeric position matrix.
#' @param storage Geometry storage label.
#' @keywords internal
#' @noRd
.shape_manipulation_summary <- function(position_matrix, storage) {
  # Resolve the radius profile when one can be derived ========================
  row_major <- identical(storage, "profile_row_major")
  radius_profile <- tryCatch(
    .shape_radius_profile(
      position_matrix = position_matrix,
      row_major = row_major,
      error_context = "manipulated geometry"
    ),
    error = function(e) NULL
  )

  # Return the manipulated geometry summary ==================================
  list(
    length = .shape_length(position_matrix = position_matrix, row_major = row_major),
    n_segments = .shape_segment_count(position_matrix, row_major = row_major),
    radius_profile = radius_profile,
    max_radius = if (is.null(radius_profile)) {
      NA_real_
    } else {
      max(radius_profile, na.rm = TRUE)
    }
  )
}

#' Update geometry metadata after a manipulation step
#' @param params Existing shape-parameter list.
#' @param summary Output of `.shape_manipulation_summary()`.
#' @param force_arbitrary Logical; whether the manipulated geometry should be
#'   relabeled as arbitrary.
#' @keywords internal
#' @noRd
.update_manipulated_shape_params <- function(params,
                                             summary,
                                             force_arbitrary = FALSE) {
  out <- params

  # Preserve the updated length and segmentation metadata =====================
  if ("length" %in% names(out)) {
    out$length <- summary$length
  }
  if ("n_segments" %in% names(out)) {
    out$n_segments <- summary$n_segments
  }

  # Refresh any stored radius-derived summaries ===============================
  if (!is.null(summary$radius_profile)) {
    if ("radius_shape" %in% names(out)) {
      out$radius_shape <- summary$radius_profile
    }
    if ("diameter_shape" %in% names(out)) {
      if (length(out$diameter_shape) > 1) {
        out$diameter_shape <- 2 * summary$radius_profile
      } else {
        out$diameter_shape <- 2 * summary$max_radius
      }
    }
    if ("radius" %in% names(out)) {
      if (length(out$radius) > 1 || all(is.na(out$radius))) {
        out$radius <- summary$radius_profile
      } else {
        out$radius <- summary$max_radius
      }
    }
    if ("diameter" %in% names(out)) {
      out$diameter <- 2 * summary$max_radius
    }
    if ("mean_radius" %in% names(out)) {
      out$mean_radius <- mean(summary$radius_profile, na.rm = TRUE)
    }
    if ("max_radius" %in% names(out)) {
      out$max_radius <- summary$max_radius
    }
    if ("semimajor_length" %in% names(out)) {
      out$semimajor_length <- summary$length / 2
    }
    if ("semiminor_length" %in% names(out)) {
      out$semiminor_length <- summary$max_radius
    }
    if ("length_radius_ratio" %in% names(out) &&
      is.finite(summary$max_radius) &&
      summary$max_radius > 0) {
      out$length_radius_ratio <- summary$length / summary$max_radius
    }
  }

  # Relabel edited canonical geometry when requested ==========================
  if (force_arbitrary && "shape" %in% names(out)) {
    out$shape <- "Arbitrary"
  }

  # Return the updated parameter list =========================================
  out
}

#' Build an Arbitrary shape after a profile-changing manipulation
#' @param position_matrix Column-major position matrix.
#' @param old_params Existing shape-parameter list.
#' @keywords internal
#' @noRd
.shape_to_arbitrary <- function(position_matrix, old_params = list()) {
  # Resolve the manipulated radius and length summaries =======================
  summary <- .shape_manipulation_summary(
    position_matrix,
    storage = "shape_column_major"
  )
  radius_profile <- summary$radius_profile

  # Return the edited geometry as an Arbitrary shape ==========================
  methods::new(
    "Arbitrary",
    position_matrix = position_matrix,
    shape_parameters = list(
      n_segments = summary$n_segments,
      length_units = old_params$length_units %||% "m",
      diameter_units = old_params$diameter_units %||%
        old_params$length_units %||%
        "m",
      radius = radius_profile %||% NA_real_,
      mean_radius = if (is.null(radius_profile)) {
        NA_real_
      } else {
        mean(radius_profile, na.rm = TRUE)
      },
      max_radius = summary$max_radius
    )
  )
}

#' Write a manipulated geometry matrix back into a Shape or Scatterer object
#' @param object Shape or Scatterer object.
#' @param info Output of `.shape_manipulation_info()`.
#' @param position_matrix Updated position matrix.
#' @param force_arbitrary Logical; whether manipulated geometry should be
#'   relabeled as arbitrary.
#' @keywords internal
#' @noRd
.shape_manipulation_apply <- function(object,
                                      info,
                                      position_matrix,
                                      force_arbitrary = FALSE) {
  # Rewrite Shape objects through their direct geometry slots =================
  if (identical(info$object_type, "Shape")) {
    if (force_arbitrary) {
      return(.shape_to_arbitrary(
        position_matrix,
        old_params = acousticTS::extract(object, "shape_parameters")
      ))
    }

    summary <- .shape_manipulation_summary(position_matrix, info$storage)
    methods::slot(object, "position_matrix") <- position_matrix
    methods::slot(object, "shape_parameters") <- .update_manipulated_shape_params(
      methods::slot(object, "shape_parameters"),
      summary,
      force_arbitrary = FALSE
    )
    return(object)
  }

  # Rewrite the requested scatterer component geometry ========================
  component_value <- methods::slot(object, info$component)
  summary <- .shape_manipulation_summary(position_matrix, info$storage)
  component_value$rpos <- position_matrix

  if ("radius" %in% names(component_value) && !is.null(summary$radius_profile)) {
    if (length(component_value$radius) > 1 || all(is.na(component_value$radius))) {
      component_value$radius <- summary$radius_profile
    } else {
      component_value$radius <- summary$max_radius
    }
  }
  if ("diameter" %in% names(component_value) && is.finite(summary$max_radius)) {
    component_value$diameter <- 2 * summary$max_radius
  }

  # Preserve any component registry mirrors ==================================
  methods::slot(object, info$component) <- component_value

  if ("components" %in% methods::slotNames(object)) {
    component_registry <- methods::slot(object, "components")
    if (info$component %in% names(component_registry)) {
      component_registry[[info$component]] <- component_value
      methods::slot(object, "components") <- component_registry
    }
  }

  # Refresh stored shape-parameter metadata ==================================
  shape_parameters <- methods::slot(object, "shape_parameters")
  if (info$component %in% names(shape_parameters) &&
    is.list(shape_parameters[[info$component]])) {
    shape_parameters[[info$component]] <- .update_manipulated_shape_params(
      shape_parameters[[info$component]],
      summary,
      force_arbitrary = force_arbitrary
    )
  } else {
    shape_parameters <- .update_manipulated_shape_params(
      shape_parameters,
      summary,
      force_arbitrary = force_arbitrary
    )
  }
  methods::slot(object, "shape_parameters") <- shape_parameters

  # Return the rewritten scatterer ============================================
  object
}

#' Re-run containment checks after an internal-component manipulation
#' @param object Scatterer object.
#' @param component Modified component name.
#' @param action Containment action.
#' @keywords internal
#' @noRd
.shape_manipulation_check_containment <- function(object,
                                                  component,
                                                  action = c(
                                                    "warn",
                                                    "error",
                                                    "ignore"
                                                  )) {
  # Resolve the requested containment action ==================================
  action <- match.arg(action)
  if (identical(action, "ignore") || !methods::is(object, "Scatterer")) {
    return(invisible(TRUE))
  }

  # Re-check swimbladder containment after editing ============================
  if (methods::is(object, "SBF") && component %in% c("body", "bladder")) {
    return(.reforge_check_internal_containment(
      acousticTS::extract(object, "body")$rpos,
      acousticTS::extract(object, "bladder")$rpos,
      component_label = "Swimbladder",
      action = action
    ))
  }

  # Re-check backbone containment after editing ===============================
  if (methods::is(object, "BBF") && component %in% c("body", "backbone")) {
    return(.reforge_check_internal_containment(
      acousticTS::extract(object, "body")$rpos,
      acousticTS::extract(object, "backbone")$rpos,
      component_label = "Backbone",
      action = action
    ))
  }

  # Return silently when no containment rule applies ==========================
  invisible(TRUE)
}

#' Translate a row-major component profile in x/z space
#' @param position_matrix Row-major profile matrix.
#' @param x_offset Along-axis translation.
#' @param z_offset Vertical translation.
#' @keywords internal
#' @noRd
.translate_profile_position_matrix <- function(position_matrix,
                                               x_offset = 0,
                                               z_offset = 0) {
  translated <- position_matrix

  # Resolve the supported row-major geometry fields ===========================
  x_idx <- .reforge_profile_row_idx(
    translated,
    .geometry_contract_schema()$profile_row_major$x
  )
  z_idx <- .reforge_profile_row_idx(
    translated,
    c("z", "z_body", "z_bladder", "z_shell", "z_fluid", "z_backbone")
  )
  zU_idx <- .reforge_profile_row_idx(
    translated,
    .geometry_contract_schema()$profile_row_major$zU
  )
  zL_idx <- .reforge_profile_row_idx(
    translated,
    .geometry_contract_schema()$profile_row_major$zL
  )

  # Apply the requested axial and vertical translations =======================
  if (length(x_idx) > 0) {
    translated[x_idx, ] <- translated[x_idx, ] + x_offset
  }
  if (length(z_idx) > 0) {
    translated[z_idx, ] <- translated[z_idx, ] + z_offset
  }
  if (length(zU_idx) > 0) {
    translated[zU_idx, ] <- translated[zU_idx, ] + z_offset
  }
  if (length(zL_idx) > 0) {
    translated[zL_idx, ] <- translated[zL_idx, ] + z_offset
  }

  # Return the translated row-major geometry ==================================
  translated
}

#' Translate a column-major shape matrix in x/y/z space
#' @param position_matrix Column-major shape matrix.
#' @param x_offset Along-axis translation.
#' @param y_offset Lateral translation.
#' @param z_offset Vertical translation.
#' @keywords internal
#' @noRd
.translate_column_position_matrix <- function(position_matrix,
                                              x_offset = 0,
                                              y_offset = 0,
                                              z_offset = 0) {
  translated <- position_matrix

  # Resolve the supported column-major geometry fields ========================
  x_candidates <- .geometry_contract_schema()$shape_column_major$x
  y_candidates <- c(
    "y", "y_body", "y_bladder", "y_shell", "y_fluid", "y_backbone"
  )
  z_candidates <- c(
    "z", "z_body", "z_bladder", "z_shell", "z_fluid", "z_backbone"
  )
  zU_candidates <- .geometry_contract_schema()$shape_column_major$zU
  zL_candidates <- .geometry_contract_schema()$shape_column_major$zL

  shift_one <- function(candidates, offset) {
    idx <- match(candidates, colnames(translated), nomatch = 0)
    idx <- idx[idx > 0]
    if (length(idx) > 0) {
      translated[, idx[1]] <<- translated[, idx[1]] + offset
      return(TRUE)
    }
    FALSE
  }

  # Apply the requested translations to available axes ========================
  shift_one(x_candidates, x_offset)
  has_y <- shift_one(y_candidates, y_offset)
  shift_one(z_candidates, z_offset)
  shift_one(zU_candidates, z_offset)
  shift_one(zL_candidates, z_offset)

  # Warn when no lateral centerline is stored ================================
  if (!has_y && !isTRUE(all.equal(y_offset, 0))) {
    warning(
      "This geometry does not store an explicit lateral centerline; 'y_offset' ",
      "was ignored.",
      call. = FALSE
    )
  }

  # Return the translated column-major geometry ===============================
  translated
}

#' Translate one supported geometry matrix
#' @param position_matrix Numeric position matrix.
#' @param storage Geometry storage label.
#' @param x_offset Along-axis translation.
#' @param y_offset Lateral translation.
#' @param z_offset Vertical translation.
#' @keywords internal
#' @noRd
.shape_translate_matrix <- function(position_matrix,
                                    storage,
                                    x_offset = 0,
                                    y_offset = 0,
                                    z_offset = 0) {
  # Dispatch the translation by stored geometry layout ========================
  if (identical(storage, "profile_row_major")) {
    if (!isTRUE(all.equal(y_offset, 0))) {
      warning(
        "Profile-style geometry does not store a lateral centerline; ",
        "'y_offset' was ignored.",
        call. = FALSE
      )
    }
    return(.translate_profile_position_matrix(
      position_matrix,
      x_offset = x_offset,
      z_offset = z_offset
    ))
  }

  # Return the translated column-major geometry ===============================
  .translate_column_position_matrix(
    position_matrix,
    x_offset = x_offset,
    y_offset = y_offset,
    z_offset = z_offset
  )
}

#' Re-anchor one supported geometry matrix along the x axis
#' @param position_matrix Numeric position matrix.
#' @param storage Geometry storage label.
#' @param anchor Anchor location.
#' @param at Target x position.
#' @keywords internal
#' @noRd
.shape_reanchor_matrix <- function(position_matrix,
                                   storage,
                                   anchor = c(
                                     "nose",
                                     "center",
                                     "tail",
                                     "max_x",
                                     "min_x"
                                   ),
                                   at = 0) {
  # Resolve the requested anchor location =====================================
  anchor <- match.arg(anchor)
  row_major <- identical(storage, "profile_row_major")
  x_vals <- .shape_x(position_matrix, row_major = row_major)
  anchor_value <- switch(anchor,
    nose = max(x_vals, na.rm = TRUE),
    max_x = max(x_vals, na.rm = TRUE),
    tail = min(x_vals, na.rm = TRUE),
    min_x = min(x_vals, na.rm = TRUE),
    center = mean(range(x_vals, na.rm = TRUE))
  )

  # Translate the geometry to the requested anchor position ===================
  .shape_translate_matrix(
    position_matrix,
    storage = storage,
    x_offset = at - anchor_value
  )
}

#' Flip one supported geometry matrix across the requested axis
#' @param position_matrix Numeric position matrix.
#' @param storage Geometry storage label.
#' @param axis Either `"x"` or `"z"`.
#' @keywords internal
#' @noRd
.shape_flip_matrix <- function(position_matrix,
                               storage,
                               axis = c("x", "z")) {
  # Resolve the requested reflection axis =====================================
  axis <- match.arg(axis)

  # Flip row-major component profiles in place ================================
  if (identical(storage, "profile_row_major")) {
    flipped <- position_matrix
    x_idx <- .reforge_profile_row_idx(
      flipped,
      .geometry_contract_schema()$profile_row_major$x
    )
    z_idx <- .reforge_profile_row_idx(
      flipped,
      c("z", "z_body", "z_bladder", "z_shell", "z_fluid", "z_backbone")
    )
    zU_idx <- .reforge_profile_row_idx(
      flipped,
      .geometry_contract_schema()$profile_row_major$zU
    )
    zL_idx <- .reforge_profile_row_idx(
      flipped,
      .geometry_contract_schema()$profile_row_major$zL
    )

    if (identical(axis, "x")) {
      # Reverse all non-x profiles to swap nose and tail ======================
      non_x <- seq_len(nrow(flipped))
      if (length(x_idx) > 0) {
        non_x <- setdiff(non_x, x_idx)
      }
      flipped[non_x, ] <- flipped[non_x, rev(seq_len(ncol(flipped))), drop = FALSE]
      return(flipped)
    }

    # Mirror the vertical fields about z = 0 ==================================
    if (length(z_idx) > 0) {
      flipped[z_idx, ] <- -flipped[z_idx, ]
    }
    if (length(zU_idx) > 0 && length(zL_idx) > 0) {
      old_zU <- flipped[zU_idx, ]
      old_zL <- flipped[zL_idx, ]
      flipped[zU_idx, ] <- -old_zL
      flipped[zL_idx, ] <- -old_zU
    }

    return(flipped)
  }

  # Flip column-major shape matrices in place =================================
  flipped <- position_matrix
  x_candidates <- .geometry_contract_schema()$shape_column_major$x
  z_candidates <- c(
    "z", "z_body", "z_bladder", "z_shell", "z_fluid", "z_backbone"
  )
  zU_candidates <- .geometry_contract_schema()$shape_column_major$zU
  zL_candidates <- .geometry_contract_schema()$shape_column_major$zL

  if (identical(axis, "x")) {
    # Reverse all non-x columns to swap nose and tail =========================
    x_idx <- match(x_candidates, colnames(flipped), nomatch = 0)
    x_idx <- x_idx[x_idx > 0]
    non_x <- seq_len(ncol(flipped))
    if (length(x_idx) > 0) {
      non_x <- setdiff(non_x, x_idx[1])
    }
    flipped[, non_x] <- flipped[rev(seq_len(nrow(flipped))), non_x, drop = FALSE]
    return(flipped)
  }

  # Mirror the vertical fields and any stored radius column ===================
  z_idx <- match(z_candidates, colnames(flipped), nomatch = 0)
  z_idx <- z_idx[z_idx > 0]
  zU_idx <- match(zU_candidates, colnames(flipped), nomatch = 0)
  zU_idx <- zU_idx[zU_idx > 0]
  zL_idx <- match(zL_candidates, colnames(flipped), nomatch = 0)
  zL_idx <- zL_idx[zL_idx > 0]

  if (length(z_idx) > 0) {
    flipped[, z_idx[1]] <- -flipped[, z_idx[1]]
  }
  if (length(zU_idx) > 0 && length(zL_idx) > 0) {
    old_zU <- flipped[, zU_idx[1]]
    old_zL <- flipped[, zL_idx[1]]
    flipped[, zU_idx[1]] <- -old_zL
    flipped[, zL_idx[1]] <- -old_zU
  }

  if ("a" %in% colnames(flipped)) {
    flipped[, "a"] <- abs(flipped[, "a"])
  }

  # Return the flipped column-major geometry ==================================
  flipped
}

#' Resample one supported geometry matrix to a new segment count
#' @param position_matrix Numeric position matrix.
#' @param storage Geometry storage label.
#' @param n_segments Desired number of intervals.
#' @keywords internal
#' @noRd
.shape_resample_matrix <- function(position_matrix,
                                   storage,
                                   n_segments) {
  # Validate the requested segment count ======================================
  if (!is.numeric(n_segments) || length(n_segments) != 1 || n_segments < 1) {
    stop("'n_segments' must be a single positive integer.", call. = FALSE)
  }

  # Dispatch resampling by stored geometry layout =============================
  n_points <- as.integer(n_segments) + 1L
  if (identical(storage, "profile_row_major")) {
    return(.resample_rpos_rows(position_matrix, n_points))
  }

  # Return the resampled column-major geometry ================================
  .resample_rpos(position_matrix, n_points)
}

#' Build a localized window used for profile inflation/tapering
#' @param x Along-axis coordinates.
#' @param x_range Length-2 numeric range.
#' @param profile Window profile.
#' @keywords internal
#' @noRd
.shape_window <- function(x,
                          x_range,
                          profile = c("cosine", "linear", "box")) {
  # Resolve the requested local-window profile ================================
  profile <- match.arg(profile)
  if (is.null(x_range)) {
    x_range <- range(x, na.rm = TRUE)
  }
  if (!is.numeric(x_range) || length(x_range) != 2) {
    stop("'x_range' must be NULL or a numeric vector of length two.", call. = FALSE)
  }

  # Validate the finite axial interval ========================================
  x0 <- min(x_range, na.rm = TRUE)
  x1 <- max(x_range, na.rm = TRUE)
  if (!is.finite(x0) || !is.finite(x1) || x1 < x0) {
    stop("'x_range' must define a finite axial interval.", call. = FALSE)
  }
  if (identical(x0, x1)) {
    # Return a point window for zero-width intervals ==========================
    idx <- which.min(abs(x - x0))
    out <- rep(0, length(x))
    out[idx] <- 1
    return(out)
  }

  # Return the requested localized taper window ===============================
  out <- rep(0, length(x))
  inside <- x >= x0 & x <= x1
  t <- (x[inside] - x0) / (x1 - x0)
  out[inside] <- switch(profile,
    box = 1,
    linear = 1 - abs(2 * t - 1),
    cosine = sin(pi * t)
  )
  out
}

#' Apply localized inflation/tapering to one supported geometry matrix
#' @param position_matrix Numeric position matrix.
#' @param storage Geometry storage label.
#' @param x_range Optional axial interval.
#' @param scale Local scale factor.
#' @param axis Local axis to scale.
#' @param profile Window profile.
#' @keywords internal
#' @noRd
.shape_local_scale_matrix <- function(position_matrix,
                                      storage,
                                      x_range = NULL,
                                      scale = 1,
                                      axis = c("radius", "width", "height"),
                                      profile = c("cosine", "linear", "box")) {
  # Resolve the requested local-scaling recipe ================================
  axis <- match.arg(axis)
  profile <- match.arg(profile)
  if (!is.numeric(scale) || length(scale) != 1 || scale <= 0) {
    stop("'scale' must be a single positive number.", call. = FALSE)
  }

  # Build the axial scaling window ============================================
  row_major <- identical(storage, "profile_row_major")
  x_vals <- .shape_x(position_matrix, row_major = row_major)
  window <- .shape_window(x_vals, x_range = x_range, profile = profile)
  factors <- 1 + (scale - 1) * window

  # Apply the localized scaling to row-major component profiles ===============
  if (row_major) {
    scaled <- position_matrix
    fields <- .profile_fields(scaled)
    width_idx <- .reforge_profile_row_idx(
      scaled,
      .geometry_contract_schema()$profile_row_major$w
    )
    z_idx <- .reforge_profile_row_idx(
      scaled,
      c("z", "z_body", "z_bladder", "z_shell", "z_fluid", "z_backbone")
    )
    zU_idx <- .reforge_profile_row_idx(
      scaled,
      .geometry_contract_schema()$profile_row_major$zU
    )
    zL_idx <- .reforge_profile_row_idx(
      scaled,
      .geometry_contract_schema()$profile_row_major$zL
    )

    if (axis %in% c("radius", "width") && length(width_idx) > 0 &&
      !all(abs(fields$w) <= sqrt(.Machine$double.eps), na.rm = TRUE)) {
      scaled[width_idx, ] <- scaled[width_idx, ] * factors
    }

    if (axis %in% c("radius", "height") && length(zU_idx) > 0 && length(zL_idx) > 0) {
      center <- .reforge_profile_centerline(fields)
      half_height <- .reforge_profile_half_height(fields) * factors
      scaled[zU_idx, ] <- center + half_height
      scaled[zL_idx, ] <- center - half_height
      if (length(z_idx) > 0) {
        scaled[z_idx, ] <- center
      }
    }

    return(scaled)
  }

  # Apply the localized scaling to column-major shape matrices ================
  scaled <- position_matrix
  width_idx <- match(
    .geometry_contract_schema()$shape_column_major$w,
    colnames(scaled),
    nomatch = 0
  )
  width_idx <- width_idx[width_idx > 0]
  z_idx <- match(
    c("z", "z_body", "z_bladder", "z_shell", "z_fluid", "z_backbone"),
    colnames(scaled),
    nomatch = 0
  )
  z_idx <- z_idx[z_idx > 0]
  zU_idx <- match(
    .geometry_contract_schema()$shape_column_major$zU,
    colnames(scaled),
    nomatch = 0
  )
  zU_idx <- zU_idx[zU_idx > 0]
  zL_idx <- match(
    .geometry_contract_schema()$shape_column_major$zL,
    colnames(scaled),
    nomatch = 0
  )
  zL_idx <- zL_idx[zL_idx > 0]
  radius_idx <- match(
    .geometry_contract_schema()$shape_column_major$radius,
    colnames(scaled),
    nomatch = 0
  )
  radius_idx <- radius_idx[radius_idx > 0]

  if (axis %in% c("radius", "width") && length(width_idx) > 0) {
    scaled[, width_idx[1]] <- scaled[, width_idx[1]] * factors
  }

  if (axis %in% c("radius", "height") && length(zU_idx) > 0 && length(zL_idx) > 0) {
    center <- if (length(z_idx) > 0) {
      scaled[, z_idx[1]]
    } else {
      (scaled[, zU_idx[1]] + scaled[, zL_idx[1]]) / 2
    }
    half_height <- ((scaled[, zU_idx[1]] - scaled[, zL_idx[1]]) / 2) * factors
    scaled[, zU_idx[1]] <- center + half_height
    scaled[, zL_idx[1]] <- center - half_height
    if (length(z_idx) > 0) {
      scaled[, z_idx[1]] <- center
    }
    if (length(radius_idx) > 0) {
      scaled[, radius_idx[1]] <- abs(half_height)
    }
  } else if (identical(axis, "radius") && length(radius_idx) > 0) {
    scaled[, radius_idx[1]] <- scaled[, radius_idx[1]] * factors
  }

  # Return the locally rescaled geometry ======================================
  scaled
}

#' Smooth a numeric profile using a centered moving average
#' @param values Numeric profile.
#' @param span Window span.
#' @param preserve_ends Logical; whether to preserve the original endpoints.
#' @keywords internal
#' @noRd
.smooth_numeric_profile <- function(values,
                                    span = 5,
                                    preserve_ends = TRUE) {
  # Validate and normalize the smoothing span =================================
  if (!is.numeric(span) || length(span) != 1 || span < 3) {
    stop("'span' must be a single integer >= 3.", call. = FALSE)
  }
  span <- as.integer(span)
  if (span %% 2 == 0) {
    span <- span + 1L
  }

  # Apply the centered moving-average smoother ================================
  weights <- rep(1 / span, span)
  smoothed <- as.numeric(stats::filter(values, weights, sides = 2))
  missing <- is.na(smoothed)
  smoothed[missing] <- values[missing]

  # Preserve the original endpoints when requested ============================
  if (isTRUE(preserve_ends)) {
    smoothed[1] <- values[1]
    smoothed[length(smoothed)] <- values[length(values)]
  }

  # Return the smoothed numeric profile =======================================
  smoothed
}

#' Smooth one supported geometry matrix
#' @param position_matrix Numeric position matrix.
#' @param storage Geometry storage label.
#' @param span Smoother span.
#' @param preserve_ends Logical; whether to preserve the original endpoints.
#' @keywords internal
#' @noRd
.shape_smooth_matrix <- function(position_matrix,
                                 storage,
                                 span = 5,
                                 preserve_ends = TRUE) {
  # Smooth row-major component profiles =======================================
  if (identical(storage, "profile_row_major")) {
    smoothed <- position_matrix
    fields <- .profile_fields(smoothed)
    width_idx <- .reforge_profile_row_idx(
      smoothed,
      .geometry_contract_schema()$profile_row_major$w
    )
    z_idx <- .reforge_profile_row_idx(
      smoothed,
      c("z", "z_body", "z_bladder", "z_shell", "z_fluid", "z_backbone")
    )
    zU_idx <- .reforge_profile_row_idx(
      smoothed,
      .geometry_contract_schema()$profile_row_major$zU
    )
    zL_idx <- .reforge_profile_row_idx(
      smoothed,
      .geometry_contract_schema()$profile_row_major$zL
    )

    if (length(width_idx) > 0 && !all(is.na(fields$w))) {
      smoothed[width_idx, ] <- .smooth_numeric_profile(
        fields$w,
        span = span,
        preserve_ends = preserve_ends
      )
    }

    if (length(zU_idx) > 0 && length(zL_idx) > 0) {
      center <- .smooth_numeric_profile(
        .reforge_profile_centerline(fields),
        span = span,
        preserve_ends = preserve_ends
      )
      half_height <- .smooth_numeric_profile(
        pmax(.reforge_profile_half_height(fields), 0),
        span = span,
        preserve_ends = preserve_ends
      )
      smoothed[zU_idx, ] <- center + half_height
      smoothed[zL_idx, ] <- center - half_height
      if (length(z_idx) > 0) {
        smoothed[z_idx, ] <- center
      }
    }

    return(smoothed)
  }

  # Smooth column-major shape matrices ========================================
  smoothed <- position_matrix
  width_idx <- match(
    .geometry_contract_schema()$shape_column_major$w,
    colnames(smoothed),
    nomatch = 0
  )
  width_idx <- width_idx[width_idx > 0]
  z_idx <- match(
    c("z", "z_body", "z_bladder", "z_shell", "z_fluid", "z_backbone"),
    colnames(smoothed),
    nomatch = 0
  )
  z_idx <- z_idx[z_idx > 0]
  zU_idx <- match(
    .geometry_contract_schema()$shape_column_major$zU,
    colnames(smoothed),
    nomatch = 0
  )
  zU_idx <- zU_idx[zU_idx > 0]
  zL_idx <- match(
    .geometry_contract_schema()$shape_column_major$zL,
    colnames(smoothed),
    nomatch = 0
  )
  zL_idx <- zL_idx[zL_idx > 0]
  radius_idx <- match(
    .geometry_contract_schema()$shape_column_major$radius,
    colnames(smoothed),
    nomatch = 0
  )
  radius_idx <- radius_idx[radius_idx > 0]

  if (length(width_idx) > 0) {
    smoothed[, width_idx[1]] <- .smooth_numeric_profile(
      smoothed[, width_idx[1]],
      span = span,
      preserve_ends = preserve_ends
    )
  }

  if (length(zU_idx) > 0 && length(zL_idx) > 0) {
    center <- if (length(z_idx) > 0) {
      smoothed[, z_idx[1]]
    } else {
      (smoothed[, zU_idx[1]] + smoothed[, zL_idx[1]]) / 2
    }
    center <- .smooth_numeric_profile(
      center,
      span = span,
      preserve_ends = preserve_ends
    )
    half_height <- .smooth_numeric_profile(
      pmax((smoothed[, zU_idx[1]] - smoothed[, zL_idx[1]]) / 2, 0),
      span = span,
      preserve_ends = preserve_ends
    )
    smoothed[, zU_idx[1]] <- center + half_height
    smoothed[, zL_idx[1]] <- center - half_height
    if (length(z_idx) > 0) {
      smoothed[, z_idx[1]] <- center
    }
    if (length(radius_idx) > 0) {
      smoothed[, radius_idx[1]] <- abs(half_height)
    }
  } else if (length(radius_idx) > 0) {
    smoothed[, radius_idx[1]] <- .smooth_numeric_profile(
      smoothed[, radius_idx[1]],
      span = span,
      preserve_ends = preserve_ends
    )
  }

  # Return the smoothed column-major geometry =================================
  smoothed
}

#' Translate a stored shape or scatterer component
#'
#' @description
#' Shift a `Shape` object or one geometry-bearing component of a `Scatterer`
#' without rebuilding it from scratch. This is most useful for re-centering
#' profiles, aligning a stored body to a preferred axial origin, or nudging an
#' internal component before model comparisons.
#'
#' @param object Shape or Scatterer object.
#' @param x_offset Along-axis translation (m).
#' @param y_offset Lateral translation (m). This is ignored for profile-style
#'   geometries that do not store an explicit lateral centerline.
#' @param z_offset Vertical translation (m).
#' @param component Optional component name for scatterers. Defaults to the
#'   primary geometry (`"body"` for most scatterers and `"shell"` for `ESS`).
#' @param containment Containment policy used when a moved swimbladder or
#'   backbone is checked against its body: `"warn"`, `"error"`, or `"ignore"`.
#'
#' @return The modified object, returned as the same broad object type.
#'
#' @examples
#' shape_obj <- cylinder(length_body = 0.05, radius_body = 0.003, n_segments = 40)
#' moved_shape <- translate_shape(shape_obj, x_offset = 0.01)
#' range(extract(moved_shape, c("position_matrix", "x")))
#'
#' @seealso [reanchor_shape()], [offset_component()], [brake()], [reforge()],
#'   [extract()]
#' @keywords shape manipulator
#' @export
translate_shape <- function(object,
                            x_offset = 0,
                            y_offset = 0,
                            z_offset = 0,
                            component = NULL,
                            containment = c("warn", "error", "ignore")) {
  # Resolve the requested geometry and containment policy =====================
  containment <- match.arg(containment)
  info <- .shape_manipulation_info(object, component = component)
  translated <- .shape_translate_matrix(
    info$position_matrix,
    storage = info$storage,
    x_offset = x_offset,
    y_offset = y_offset,
    z_offset = z_offset
  )

  # Apply the translated geometry and re-check containment ====================
  object <- .shape_manipulation_apply(object, info, translated)
  .shape_manipulation_check_containment(
    object,
    component = info$component %||% .default_shape_component(object),
    action = containment
  )
  # Return the translated object ==============================================
  object
}

#' Re-anchor a stored shape or scatterer component along the x axis
#'
#' @description
#' Translate a shape or scatterer component so that its nose, center, or tail
#' lies at a specified x position.
#'
#' @inheritParams translate_shape
#' @param anchor Anchor location used to define the translation target. `nose`
#'   is treated as the maximum stored x position and `tail` as the minimum.
#' @param at Target x location (m).
#'
#' @return The modified object, returned as the same broad object type.
#'
#' @examples
#' shape_obj <- prolate_spheroid(length_body = 0.04, radius_body = 0.004)
#' centered_shape <- reanchor_shape(shape_obj, anchor = "center", at = 0)
#' range(extract(centered_shape, c("position_matrix", "x")))
#'
#' @seealso [translate_shape()], [brake()], [reforge()], [extract()]
#' @keywords shape manipulator
#' @export
reanchor_shape <- function(object,
                           anchor = c("nose", "center", "tail", "max_x", "min_x"),
                           at = 0,
                           component = NULL,
                           containment = c("warn", "error", "ignore")) {
  # Resolve the requested geometry and anchor target ==========================
  containment <- match.arg(containment)
  info <- .shape_manipulation_info(object, component = component)
  anchored <- .shape_reanchor_matrix(
    info$position_matrix,
    storage = info$storage,
    anchor = anchor,
    at = at
  )

  # Apply the anchored geometry and re-check containment ======================
  object <- .shape_manipulation_apply(object, info, anchored)
  .shape_manipulation_check_containment(
    object,
    component = info$component %||% .default_shape_component(object),
    action = containment
  )
  # Return the anchored object ================================================
  object
}

#' Offset an internal scatterer component without rebuilding the object
#'
#' @description
#' Translate an internal geometry-bearing component such as a swimbladder or
#' backbone while leaving the outer body unchanged.
#'
#' @param object Scatterer object containing an internal component.
#' @param component Internal component name. Currently most useful for
#'   `"bladder"` and `"backbone"`.
#' @param x_offset Along-axis translation (m).
#' @param z_offset Vertical translation (m).
#' @param containment Containment policy used when the shifted component is
#'   checked against the body: `"warn"`, `"error"`, or `"ignore"`.
#'
#' @return The modified scatterer object.
#'
#' @examples
#' fish <- sbf_generate(
#'   x_body = c(0, 0.1),
#'   w_body = c(0.006, 0.008),
#'   zU_body = c(0.001, 0.002),
#'   zL_body = c(-0.001, -0.002),
#'   x_bladder = c(0.02, 0.08),
#'   w_bladder = c(0, 0),
#'   zU_bladder = c(0.0012, 0.0012),
#'   zL_bladder = c(-0.0012, -0.0012),
#'   density_body = 1040,
#'   density_bladder = 1.2,
#'   sound_speed_body = 1500,
#'   sound_speed_bladder = 340
#' )
#' shifted_fish <- offset_component(fish, component = "bladder", x_offset = 0.003)
#' min(extract(shifted_fish, c("bladder", "rpos", "x_bladder")))
#'
#' @seealso [translate_shape()], [reanchor_shape()], [reforge()]
#' @keywords shape manipulator
#' @export
offset_component <- function(object,
                             component = c("bladder", "backbone"),
                             x_offset = 0,
                             z_offset = 0,
                             containment = c("warn", "error", "ignore")) {
  # Validate the scatterer-only manipulation pathway ==========================
  if (!methods::is(object, "Scatterer")) {
    stop("`offset_component()` expects a Scatterer object.", call. = FALSE)
  }

  # Resolve the requested internal component ==================================
  component <- match.arg(component)
  # Dispatch through translate_shape() to preserve one code path ==============
  translate_shape(
    object,
    x_offset = x_offset,
    z_offset = z_offset,
    component = component,
    containment = containment
  )
}

#' Locally widen, pinch, or taper a stored shape profile
#'
#' @description
#' Apply a localized scaling window to a stored profile. Values of `scale > 1`
#' inflate the selected region, while `scale < 1` pinch or taper it. For
#' canonical `Shape` objects, the result is returned as an `Arbitrary` shape
#' because the edited profile is no longer guaranteed to preserve the canonical
#' class geometry.
#'
#' @inheritParams translate_shape
#' @param x_range Optional axial interval over which the local manipulation is
#'   applied.
#' @param scale Positive local scale factor.
#' @param axis Which profile dimension to scale.
#' @param profile Local window profile. `"cosine"` gives a smooth bump centered
#'   inside `x_range`, `"linear"` gives a triangular bump, and `"box"` applies a
#'   uniform factor inside the interval.
#'
#' @return The modified object. `Shape` inputs that are locally reshaped are
#'   returned as `Arbitrary` shapes.
#'
#' @examples
#' shape_obj <- cylinder(length_body = 0.05, radius_body = 0.003, n_segments = 80)
#' pinched_shape <- inflate_shape(
#'   shape_obj,
#'   x_range = c(0.015, 0.035),
#'   scale = 0.7
#' )
#' max(extract(pinched_shape, c("shape_parameters", "max_radius")))
#'
#' @seealso [smooth_shape()], [resample_shape()], [flip_shape()], [brake()],
#'   [reforge()]
#' @keywords shape manipulator
#' @export
inflate_shape <- function(object,
                          x_range = NULL,
                          scale = 1,
                          component = NULL,
                          axis = c("radius", "width", "height"),
                          profile = c("cosine", "linear", "box"),
                          containment = c("warn", "error", "ignore")) {
  # Resolve the requested geometry edit =======================================
  containment <- match.arg(containment)
  info <- .shape_manipulation_info(object, component = component)
  scaled <- .shape_local_scale_matrix(
    info$position_matrix,
    storage = info$storage,
    x_range = x_range,
    scale = scale,
    axis = axis,
    profile = profile
  )

  # Apply the edited profile and relabel canonical shapes when needed =========
  object <- .shape_manipulation_apply(
    object,
    info,
    scaled,
    force_arbitrary = TRUE
  )
  # Re-check internal containment after the local profile edit ================
  .shape_manipulation_check_containment(
    object,
    component = info$component %||% .default_shape_component(object),
    action = containment
  )
  # Return the edited object ==================================================
  object
}

#' Smooth a stored shape or scatterer profile
#'
#' @description
#' Smooth the stored geometry using a centered moving-average filter applied to
#' profile coordinates. This is useful for cleaning digitized outlines before
#' canonicalization or model runs. For canonical `Shape` objects, the result is
#' returned as an `Arbitrary` shape because the edited profile is no longer
#' guaranteed to preserve the canonical class geometry.
#'
#' @inheritParams translate_shape
#' @param span Centered moving-average span. Even values are rounded up to the
#'   next odd integer.
#' @param preserve_ends Logical; whether to keep the first and last profile
#'   points fixed.
#'
#' @return The modified object. `Shape` inputs that are smoothed are returned as
#'   `Arbitrary` shapes.
#'
#' @examples
#' shape_obj <- arbitrary(
#'   x_body = c(0, 0.01, 0.02, 0.03, 0.04),
#'   zU_body = c(0, 0.003, 0.006, 0.0035, 0),
#'   zL_body = c(0, -0.0025, -0.0055, -0.003, 0)
#' )
#' smoothed_shape <- smooth_shape(shape_obj, span = 3)
#' extract(smoothed_shape, c("shape_parameters", "n_segments"))
#'
#' @seealso [inflate_shape()], [resample_shape()], [brake()], [reforge()]
#' @keywords shape manipulator
#' @export
smooth_shape <- function(object,
                         span = 5,
                         component = NULL,
                         preserve_ends = TRUE,
                         containment = c("warn", "error", "ignore")) {
  # Resolve the requested smoothing pathway ===================================
  containment <- match.arg(containment)
  info <- .shape_manipulation_info(object, component = component)
  smoothed <- .shape_smooth_matrix(
    info$position_matrix,
    storage = info$storage,
    span = span,
    preserve_ends = preserve_ends
  )

  # Apply the smoothed profile and relabel canonical shapes when needed =======
  object <- .shape_manipulation_apply(
    object,
    info,
    smoothed,
    force_arbitrary = TRUE
  )
  # Re-check internal containment after smoothing =============================
  .shape_manipulation_check_containment(
    object,
    component = info$component %||% .default_shape_component(object),
    action = containment
  )
  # Return the smoothed object ================================================
  object
}

#' Resample a stored shape or scatterer profile to a new segment count
#'
#' @description
#' Re-discretize a shape or scatterer component along its x axis without
#' rebuilding it from scratch.
#'
#' @inheritParams translate_shape
#' @param n_segments New number of intervals used to discretize the profile.
#'
#' @return The modified object, returned as the same broad object type.
#'
#' @examples
#' shape_obj <- sphere(radius_body = 0.01, n_segments = 20)
#' shape_fine <- resample_shape(shape_obj, n_segments = 80)
#' extract(shape_fine, c("shape_parameters", "n_segments"))
#'
#' @seealso [smooth_shape()], [inflate_shape()], [reforge()], [extract()]
#' @keywords shape manipulator
#' @export
resample_shape <- function(object,
                           n_segments,
                           component = NULL,
                           containment = c("warn", "error", "ignore")) {
  # Resolve the requested resampling pathway ==================================
  containment <- match.arg(containment)
  info <- .shape_manipulation_info(object, component = component)
  resampled <- .shape_resample_matrix(
    info$position_matrix,
    storage = info$storage,
    n_segments = n_segments
  )

  # Apply the resampled geometry and re-check containment =====================
  object <- .shape_manipulation_apply(object, info, resampled)
  .shape_manipulation_check_containment(
    object,
    component = info$component %||% .default_shape_component(object),
    action = containment
  )
  # Return the resampled object ===============================================
  object
}

#' Flip a stored shape or scatterer profile across the x or z axis
#'
#' @description
#' Reverse axial orientation (`axis = "x"`) or mirror the geometry
#' dorsoventrally (`axis = "z"`) without manually editing the position matrix.
#'
#' @inheritParams translate_shape
#' @param axis Flip axis. Use `"x"` to reverse nose/tail orientation while
#'   keeping the existing x grid, or `"z"` to mirror the profile vertically.
#'
#' @return The modified object, returned as the same broad object type.
#'
#' @examples
#' shape_obj <- arbitrary(
#'   x_body = c(0, 0.01, 0.02, 0.03),
#'   radius_body = c(0, 0.004, 0.002, 0)
#' )
#' flipped_shape <- flip_shape(shape_obj, axis = "x")
#' head(extract(flipped_shape, c("position_matrix", "zU")))
#'
#' @seealso [translate_shape()], [reanchor_shape()], [inflate_shape()],
#'   [smooth_shape()]
#' @keywords shape manipulator
#' @export
flip_shape <- function(object,
                       axis = c("x", "z"),
                       component = NULL,
                       containment = c("warn", "error", "ignore")) {
  # Resolve the requested reflection pathway ==================================
  containment <- match.arg(containment)
  info <- .shape_manipulation_info(object, component = component)
  flipped <- .shape_flip_matrix(
    info$position_matrix,
    storage = info$storage,
    axis = axis
  )

  # Apply the flipped geometry and re-check containment =======================
  object <- .shape_manipulation_apply(object, info, flipped)
  .shape_manipulation_check_containment(
    object,
    component = info$component %||% .default_shape_component(object),
    action = containment
  )
  # Return the flipped object =================================================
  object
}
