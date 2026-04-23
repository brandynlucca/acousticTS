################################################################################
# Transition matrix method (TMM) bistatic-summary helpers
################################################################################

# Convert spherical angles to a 3D unit vector.
#' @noRd
.tmm_spherical_to_cartesian <- function(theta, phi) {
  # Return the 3D unit-vector components =======================================
  c(
    sin(theta) * cos(phi),
    sin(theta) * sin(phi),
    cos(theta)
  )
}

# Convert one 3D direction vector back into spherical angles.
#' @noRd
.tmm_cartesian_to_spherical <- function(vec) {
  # Normalize the input vector =================================================
  vec <- as.numeric(vec)
  vec <- vec / sqrt(sum(vec^2))
  # Convert the normalized vector to spherical angles ==========================
  theta <- acos(max(-1, min(1, vec[3])))
  phi <- atan2(vec[2], vec[1])
  if (phi < 0) {
    phi <- phi + 2 * pi
  }
  c(theta = theta, phi = phi)
}

# Compute the cross product of two 3D vectors.
#' @noRd
.tmm_cross_product <- function(a, b) {
  # Evaluate the 3D cross product ==============================================
  c(
    a[2] * b[3] - a[3] * b[2],
    a[3] * b[1] - a[1] * b[3],
    a[1] * b[2] - a[2] * b[1]
  )
}

# Build one orthonormal basis around the forward-scattering direction so that
# local great-circle slices can be defined independently of global coordinates.
#' @noRd
.tmm_forward_basis <- function(theta_body, phi_body) {
  # Build the forward-scattering direction =====================================
  forward <- .tmm_spherical_to_cartesian(theta_body, phi_body)
  body_axis <- c(1, 0, 0)
  e1 <- body_axis - sum(body_axis * forward) * forward
  if (sqrt(sum(e1^2)) < 1e-12) {
    ref <- if (abs(forward[3]) < 0.95) c(0, 0, 1) else c(0, 1, 0)
    e1 <- ref - sum(ref * forward) * forward
  }
  e1 <- e1 / sqrt(sum(e1^2))
  e2 <- .tmm_cross_product(forward, e1)
  e2 <- e2 / sqrt(sum(e2^2))

  list(forward = forward, e1 = e1, e2 = e2)
}

# Build one full incident-centered local scattering grid using the shared
# forward/body-axis basis.
#' @noRd
.tmm_local_grid_directions <- function(theta_body,
                                       phi_body,
                                       psi_scatter,
                                       alpha_scatter) {
  basis <- .tmm_forward_basis(theta_body, phi_body)
  psi_mesh <- rep(psi_scatter, times = length(alpha_scatter))
  alpha_mesh <- rep(alpha_scatter, each = length(psi_scatter))

  dirs <- vapply(
    seq_along(psi_mesh),
    function(i) {
      basis$forward * cos(psi_mesh[i]) +
        (basis$e1 * cos(alpha_mesh[i]) + basis$e2 * sin(alpha_mesh[i])) *
          sin(psi_mesh[i])
    },
    numeric(3)
  )
  angles <- apply(dirs, 2, .tmm_cartesian_to_spherical)

  data.frame(
    psi_scatter = psi_mesh,
    alpha_scatter = alpha_mesh,
    theta_scatter = as.numeric(angles["theta", ]),
    phi_scatter = as.numeric(angles["phi", ])
  )
}

# Build one great-circle scattering slice in local coordinates centered on the
# forward-scattering direction.
#' @noRd
.tmm_local_slice_directions <- function(theta_body,
                                        phi_body,
                                        psi_scatter,
                                        alpha = 0) {
  slice <- .tmm_local_grid_directions(
    theta_body = theta_body,
    phi_body = phi_body,
    psi_scatter = psi_scatter,
    alpha_scatter = alpha
  )

  data.frame(
    psi_scatter = slice$psi_scatter,
    theta_scatter = slice$theta_scatter,
    phi_scatter = slice$phi_scatter
  )
}

# Compute the angular separation from the forward-scattering direction for one
# scattering grid.
#' @noRd
.tmm_forward_separation_matrix <- function(theta_body,
                                           phi_body,
                                           theta_scatter,
                                           phi_scatter) {
  # Build the forward direction for the incident geometry ======================
  forward <- .tmm_spherical_to_cartesian(theta_body, phi_body)
  # Compute the forward-centered angular separation grid =======================
  psi <- outer(
    theta_scatter,
    phi_scatter,
    Vectorize(function(theta_val, phi_val) {
      scat <- .tmm_spherical_to_cartesian(theta_val, phi_val)
      acos(max(-1, min(1, sum(forward * scat))))
    })
  )
  psi
}

# Measure the -drop_dB width of the local backscatter lobe on one great-circle
# slice expressed in forward-centered scattering angle.
#' @noRd
.tmm_lobe_width <- function(psi_scatter, sigma_scat_dB, center_psi = pi,
                            drop_dB = 3) {
  # Validate the requested dB drop threshold ===================================
  if (!is.numeric(drop_dB) || length(drop_dB) != 1 || !is.finite(drop_dB) ||
    drop_dB <= 0) {
    stop("'drop_dB' must be a single positive numeric value.", call. = FALSE)
  }

  idx_center <- which.min(abs(psi_scatter - center_psi))
  peak_level <- sigma_scat_dB[idx_center]
  if (!is.finite(peak_level)) {
    return(NA_real_)
  }
  # Identify the contiguous lobe around the center angle =======================
  threshold <- peak_level - drop_dB
  keep <- which(is.finite(sigma_scat_dB) & sigma_scat_dB >= threshold)
  if (!(idx_center %in% keep)) {
    return(NA_real_)
  }

  left <- idx_center
  while (left > 1 && ((left - 1) %in% keep)) {
    left <- left - 1
  }
  right <- idx_center
  while (right < length(psi_scatter) && ((right + 1) %in% keep)) {
    right <- right + 1
  }

  psi_scatter[right] - psi_scatter[left]
}
# Build one named local great-circle slice relative to the forward direction.
#' @noRd
.tmm_local_slice <- function(model_params,
                             frequency_idx,
                             shape_parameters,
                             theta_body,
                             phi_body,
                             psi_scatter,
                             alpha,
                             name) {
  # Build the requested local slice directions =================================
  slice_dirs <- .tmm_local_slice_directions(
    theta_body = theta_body,
    phi_body = phi_body,
    psi_scatter = psi_scatter,
    alpha = alpha
  )
  n_eval <- nrow(slice_dirs)
  # Evaluate scattering along the local slice ==================================
  f_scat <- .tmm_scattering_points(
    model_params = model_params,
    frequency_idx = frequency_idx,
    shape_parameters = shape_parameters,
    theta_body = rep(theta_body, n_eval),
    phi_body = rep(phi_body, n_eval),
    theta_scatter = slice_dirs$theta_scatter,
    phi_scatter = slice_dirs$phi_scatter
  )
  # Return the labeled slice data ==============================================
  data.frame(
    slice = name,
    psi_scatter = slice_dirs$psi_scatter,
    theta_scatter = slice_dirs$theta_scatter,
    phi_scatter = slice_dirs$phi_scatter,
    f_scat = f_scat,
    sigma_scat = .sigma_bs(f_scat),
    sigma_scat_dB = db(.sigma_bs(f_scat))
  )
}

# Default coarse angular sectors used by the bistatic summary helper.
#' @noRd
.tmm_default_sectors <- function() {
  data.frame(
    sector = c("forward_sector", "oblique_sector", "backward_sector"),
    psi_min = c(0, pi / 3, 2 * pi / 3),
    psi_max = c(pi / 3, 2 * pi / 3, pi),
    stringsAsFactors = FALSE
  )
}

# Validate one user-supplied angular-sector table.
#' @noRd
.tmm_validate_sectors <- function(sectors) {
  # Fall back to the package-default sectors ===================================
  if (is.null(sectors)) {
    return(.tmm_default_sectors())
  }
  # Validate the required columns and bounds ===================================
  required_cols <- c("sector", "psi_min", "psi_max")
  missing_cols <- setdiff(required_cols, names(sectors))
  if (length(missing_cols)) {
    stop(
      "'sectors' must contain the columns: ",
      paste(required_cols, collapse = ", "),
      ". Missing: ",
      paste(missing_cols, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  if (!is.character(sectors$sector) || any(!nzchar(sectors$sector))) {
    stop("'sectors$sector' must contain non-empty names.", call. = FALSE)
  }
  if (!is.numeric(sectors$psi_min) || !is.numeric(sectors$psi_max) ||
    any(!is.finite(sectors$psi_min)) || any(!is.finite(sectors$psi_max)) ||
    any(sectors$psi_min < 0) || any(sectors$psi_max > pi) ||
    any(sectors$psi_min >= sectors$psi_max)) {
    stop(
      "'sectors' must define finite angular bounds in radians with 0 <= ",
      "psi_min < psi_max <= pi.",
      call. = FALSE
    )
  }

  sectors
}

# Compute solid-angle weights for one theta-phi scattering grid.
#' @noRd
.tmm_grid_solid_angle <- function(theta_scatter, phi_scatter) {
  # Resolve the implied grid-cell edges ========================================
  theta_edges <- .tmm_grid_edges(theta_scatter, lower = 0, upper = pi)
  phi_edges <- .tmm_grid_edges(phi_scatter, lower = 0, upper = 2 * pi)
  # Convert the theta-phi cells to solid-angle weights =========================
  dphi <- diff(phi_edges)
  theta_band <- cos(theta_edges[-length(theta_edges)]) - cos(theta_edges[-1])
  tcrossprod(theta_band, dphi)
}

# Resolve and validate the core control inputs for `tmm_bistatic_summary()`.
#' @noRd
.tmm_bistatic_summary_inputs <- function(model_params,
                                         frequency,
                                         theta_body,
                                         phi_body,
                                         n_psi,
                                         sectors) {
  # Reject unsupported retained cylindrical bistatic workflows ================
  if (model_params$parameters$coordinate_system %in%
    c("cylindrical", "axisymmetric", "cylinder_native")) {
    .tmm_stop_cylinder_bistatic_public("tmm_bistatic_summary()")
  }

  # Validate the local-slice resolution and resolve stored defaults ============
  if (!is.numeric(n_psi) || length(n_psi) != 1 || !is.finite(n_psi) ||
    n_psi < 3 || n_psi %% 1 != 0) {
    stop("'n_psi' must be a single integer >= 3.", call. = FALSE)
  }

  acoustics <- model_params$parameters$acoustics
  defaults <- model_params$body

  list(
    idx = .tmm_plot_frequency_index(frequency, acoustics$frequency),
    theta_body = .tmm_scalar_angle(theta_body, defaults$theta_body, "theta_body"),
    phi_body = .tmm_scalar_angle(phi_body, defaults$phi_body %||% pi, "phi_body"),
    sectors = .tmm_validate_sectors(sectors)
  )
}

# Build the retained scattering grid and angular quadrature helpers used by the
# bistatic summary.
#' @noRd
.tmm_bistatic_summary_grid <- function(object,
                                       acoustics,
                                       idx,
                                       theta_body,
                                       phi_body,
                                       n_theta,
                                       n_phi) {
  # Reuse the stored scattering grid at the requested incident direction =======
  grid <- suppressWarnings(tmm_scattering_grid(
    object = object,
    frequency = acoustics$frequency[idx],
    theta_body = theta_body,
    phi_body = phi_body,
    n_theta = n_theta,
    n_phi = n_phi
  ))
  psi_grid <- .tmm_forward_separation_matrix(
    theta_body = theta_body,
    phi_body = phi_body,
    theta_scatter = grid$theta_scatter,
    phi_scatter = grid$phi_scatter
  )

  list(
    grid = grid,
    psi_grid = psi_grid,
    solid_angle = .tmm_grid_solid_angle(grid$theta_scatter, grid$phi_scatter)
  )
}

# Build the forward-centered diagnostic slices used by the bistatic summary.
#' @noRd
.tmm_bistatic_summary_slices <- function(model_params,
                                         idx,
                                         shape_parameters,
                                         theta_body,
                                         phi_body,
                                         n_psi) {
  # Construct the shared forward-angle grid for both local slices =============
  psi_scatter <- seq(0, pi, length.out = n_psi)

  list(
    psi_scatter = psi_scatter,
    forward_slice = .tmm_local_slice(
      model_params = model_params,
      frequency_idx = idx,
      shape_parameters = shape_parameters,
      theta_body = theta_body,
      phi_body = phi_body,
      psi_scatter = psi_scatter,
      alpha = 0,
      name = "forward_scatter_slice"
    ),
    dorsal_slice = .tmm_local_slice(
      model_params = model_params,
      frequency_idx = idx,
      shape_parameters = shape_parameters,
      theta_body = theta_body,
      phi_body = phi_body,
      psi_scatter = psi_scatter,
      alpha = pi / 2,
      name = "dorsal_ventral_slice"
    )
  )
}

# Summarize the peak-scattering cell together with the exact forward and
# monostatic evaluation points.
#' @noRd
.tmm_bistatic_summary_points <- function(object,
                                         grid,
                                         psi_grid,
                                         theta_body,
                                         phi_body) {
  # Locate the peak retained scattering cell ==================================
  peak_idx <- arrayInd(which.max(grid$sigma_scat), dim(grid$sigma_scat))
  peak <- list(
    theta = grid$theta_scatter[peak_idx[1]],
    phi = grid$phi_scatter[peak_idx[2]],
    sigma = grid$sigma_scat[peak_idx[1], peak_idx[2]],
    psi = psi_grid[peak_idx[1], peak_idx[2]]
  )

  # Evaluate the exact forward and monostatic receive directions ==============
  list(
    peak = peak,
    forward_point = suppressWarnings(tmm_scattering(
      object = object,
      theta_body = theta_body,
      phi_body = phi_body,
      theta_scatter = theta_body,
      phi_scatter = phi_body
    )),
    backscatter_point = suppressWarnings(tmm_scattering(
      object = object,
      theta_body = theta_body,
      phi_body = phi_body,
      theta_scatter = pi - theta_body,
      phi_scatter = phi_body + pi
    ))
  )
}

# Build one logical mask for a named forward-separation sector.
#' @noRd
.tmm_bistatic_sector_mask <- function(psi_grid, sectors, i) {
  # Keep the upper bound closed only for the final sector ======================
  if (i == nrow(sectors)) {
    return(psi_grid >= sectors$psi_min[i] & psi_grid <= sectors$psi_max[i])
  }

  psi_grid >= sectors$psi_min[i] & psi_grid < sectors$psi_max[i]
}

# Integrate user-defined angular sectors over the retained bistatic grid.
#' @noRd
.tmm_bistatic_sector_integrals <- function(sectors, psi_grid, sigma_scat, solid_angle) {
  # Integrate each requested sector over the retained grid =====================
  do.call(
    rbind,
    lapply(
      seq_len(nrow(sectors)),
      function(i) {
        keep <- .tmm_bistatic_sector_mask(psi_grid, sectors, i)
        data.frame(
          sector = sectors$sector[i],
          psi_min = sectors$psi_min[i],
          psi_max = sectors$psi_max[i],
          integrated_sigma_scat = sum(sigma_scat[keep] * solid_angle[keep])
        )
      }
    )
  )
}

#' Summarize bistatic products from a stored TMM object
#'
#' @description
#' Reuses the stored T-matrix blocks at one stored frequency to compute a
#' higher-level bistatic summary, including forward- and cross-plane slices,
#' peak-scattering direction, backscatter-lobe width, and integrated scattering
#' over coarse angular sectors. In the current package build, this helper is
#' available for the spherical and spheroidal stored branches. Stored cylinders
#' are intentionally outside the public scope of this helper and stop with an
#' explicit error.
#'
#' @param object Scatterer-object previously evaluated with
#'   `target_strength(..., model = "TMM", store_t_matrix = TRUE)`.
#' @param frequency Stored frequency (Hz) to summarize. Required when the object
#'   contains more than one stored frequency.
#' @param theta_body Incident polar angle (radians). Defaults to the stored TMM
#'   incident angle.
#' @param phi_body Incident azimuth angle (radians). Defaults to the stored TMM
#'   incident angle.
#' @param n_theta Number of receive polar-angle samples used by the summary
#' grid.
#' @param n_phi Number of receive azimuth samples used by the summary grid.
#' @param n_psi Number of forward-centered angular samples used for the local
#'   great-circle slices.
#' @param sectors Optional data frame with columns `sector`, `psi_min`, and
#'   `psi_max` (radians). When omitted, three coarse forward/oblique/backward
#'   sectors are used.
#' @param drop_dB Positive dB drop used to define the backscatter-lobe width on
#'   the forward-centered slice.
#' @param include_grid Logical; include the underlying scattering grid in the
#'   returned list.
#'
#' @return A list containing scalar summary metrics, the named slice data
#'   frames, sector integrals, and optionally the underlying scattering grid.
#'
#' @seealso \code{\link{tmm_scattering_grid}},
#'   \code{\link{tmm_products}}
#' @export
tmm_bistatic_summary <- function(object,
                                 frequency = NULL,
                                 theta_body = NULL,
                                 phi_body = NULL,
                                 n_theta = 91,
                                 n_phi = 181,
                                 n_psi = 181,
                                 sectors = NULL,
                                 drop_dB = 3,
                                 include_grid = FALSE) {
  # Recover the stored TMM state ===============================================
  model_params <- .tmm_require_stored_blocks(object)
  .tmm_warn_exploratory_cylinder_blocks(object, model_params)
  acoustics <- model_params$parameters$acoustics
  shape_parameters <- acousticTS::extract(object, "shape_parameters")
  summary_inputs <- .tmm_bistatic_summary_inputs(
    model_params = model_params,
    frequency = frequency,
    theta_body = theta_body,
    phi_body = phi_body,
    n_psi = n_psi,
    sectors = sectors
  )
  idx <- summary_inputs$idx
  theta_body <- summary_inputs$theta_body
  phi_body <- summary_inputs$phi_body
  sectors <- summary_inputs$sectors
  # Build the scattering grid and local slices =================================
  grid_state <- .tmm_bistatic_summary_grid(
    object = object,
    acoustics = acoustics,
    idx = idx,
    theta_body = theta_body,
    phi_body = phi_body,
    n_theta = n_theta,
    n_phi = n_phi
  )
  grid <- grid_state$grid
  psi_grid <- grid_state$psi_grid
  solid_angle <- grid_state$solid_angle
  # Build the forward-centered diagnostic slices ===============================
  slice_state <- .tmm_bistatic_summary_slices(
    model_params = model_params,
    idx = idx,
    shape_parameters = shape_parameters,
    theta_body = theta_body,
    phi_body = phi_body,
    n_psi = n_psi
  )
  psi_scatter <- slice_state$psi_scatter
  forward_slice <- slice_state$forward_slice
  dorsal_slice <- slice_state$dorsal_slice
  # Summarize the peak scattering direction and monostatic points ==============
  point_state <- .tmm_bistatic_summary_points(
    object = object,
    grid = grid,
    psi_grid = psi_grid,
    theta_body = theta_body,
    phi_body = phi_body
  )
  peak_theta <- point_state$peak$theta
  peak_phi <- point_state$peak$phi
  peak_sigma <- point_state$peak$sigma
  peak_psi <- point_state$peak$psi
  forward_point <- point_state$forward_point
  backscatter_point <- point_state$backscatter_point
  # Integrate the requested angular sectors ====================================
  sector_integrals <- .tmm_bistatic_sector_integrals(
    sectors = sectors,
    psi_grid = psi_grid,
    sigma_scat = grid$sigma_scat,
    solid_angle = solid_angle
  )
  # Assemble the scalar metrics and requested outputs ==========================
  metrics <- data.frame(
    frequency = grid$frequency,
    forward_sigma_scat = forward_point$sigma_scat[idx],
    forward_sigma_scat_dB = forward_point$sigma_scat_dB[idx],
    sigma_bs = backscatter_point$sigma_scat[idx],
    TS = db(backscatter_point$sigma_scat[idx]),
    peak_theta_scatter = peak_theta,
    peak_phi_scatter = peak_phi,
    peak_psi_scatter = peak_psi,
    peak_sigma_scat = peak_sigma,
    peak_sigma_scat_dB = db(peak_sigma),
    backscatter_lobe_width = .tmm_lobe_width(
      psi_scatter = forward_slice$psi_scatter,
      sigma_scat_dB = forward_slice$sigma_scat_dB,
      center_psi = pi,
      drop_dB = drop_dB
    )
  )

  out <- list(
    metrics = metrics,
    slices = list(
      forward_scatter = forward_slice,
      dorsal_ventral = dorsal_slice
    ),
    sector_integrals = sector_integrals
  )
  if (isTRUE(include_grid)) {
    out$grid <- grid
  }
  # Return the bistatic summary bundle =========================================
  out
}

#' Collect multiple post-processed products from one stored TMM solve
#'
#' @description
#' Provides a higher-level convenience interface for the stored-block TMM
#' workflow. One call can return the monostatic scattering spectrum, an
#' orientation average, and a higher-level bistatic summary without rebuilding
#' the underlying T-matrix solve. For stored cylinders, the currently supported
#' products are the exact monostatic scattering spectrum and the corresponding
#' orientation-averaged monostatic outputs; cylinder bistatic summaries are
#' intentionally outside the public scope.
#'
#' @param object Scatterer-object previously evaluated with
#'   `target_strength(..., model = "TMM", store_t_matrix = TRUE)`.
#' @param frequency Stored frequency (Hz) used when `bistatic_summary = TRUE`.
#' @param theta_body Incident polar angle (radians) for the monostatic and
#'   bistatic products. Defaults to the stored TMM incident angle.
#' @param phi_body Incident azimuth angle (radians) for the monostatic and
#'   bistatic products. Defaults to the stored TMM incident angle.
#' @param orientation Optional orientation distribution created by
#'   \code{\link{tmm_orientation_distribution}}.
#' @param bistatic_summary Logical; include the output of
#'   \code{\link{tmm_bistatic_summary}}.
#' @param include_grid Logical; keep the 2D scattering grid inside the bistatic
#'   summary output.
#' @param n_theta Number of receive polar-angle samples for the bistatic
#'   summary.
#' @param n_phi Number of receive azimuth samples for the bistatic summary.
#' @param n_psi Number of forward-centered angular samples used by the local
#'   summary slices.
#' @param sectors Optional angular-sector table passed to
#'   \code{\link{tmm_bistatic_summary}}.
#' @param drop_dB Positive dB drop used to define the backscatter-lobe width.
#'
#' @return A named list containing the requested post-processed TMM products.
#'
#' @seealso \code{\link{tmm_scattering}},
#'   \code{\link{tmm_average_orientation}},
#'   \code{\link{tmm_bistatic_summary}}
#' @export
tmm_products <- function(object,
                         frequency = NULL,
                         theta_body = NULL,
                         phi_body = NULL,
                         orientation = NULL,
                         bistatic_summary = FALSE,
                         include_grid = FALSE,
                         n_theta = 91,
                         n_phi = 181,
                         n_psi = 181,
                         sectors = NULL,
                         drop_dB = 3) {
  # Recover the stored TMM state and defaults ==================================
  model_params <- .tmm_require_stored_blocks(object)
  .tmm_warn_exploratory_cylinder_blocks(object, model_params)
  defaults <- model_params$body
  theta_body <- .tmm_scalar_angle(theta_body, defaults$theta_body, "theta_body")
  phi_body <- .tmm_scalar_angle(phi_body, defaults$phi_body %||% pi, "phi_body")
  # Always return the exact stored monostatic scattering =======================
  out <- list(
    monostatic = suppressWarnings(tmm_scattering(
      object = object,
      theta_body = theta_body,
      phi_body = phi_body,
      theta_scatter = pi - theta_body,
      phi_scatter = phi_body + pi
    ))
  )
  # Add the requested orientation average ======================================
  if (!is.null(orientation)) {
    out$orientation_average <- suppressWarnings(tmm_average_orientation(
      object = object,
      distribution = orientation
    ))
  }
  # Add the requested bistatic summary =========================================
  if (isTRUE(bistatic_summary)) {
    out$bistatic_summary <- suppressWarnings(tmm_bistatic_summary(
      object = object,
      frequency = frequency,
      theta_body = theta_body,
      phi_body = phi_body,
      n_theta = n_theta,
      n_phi = n_phi,
      n_psi = n_psi,
      sectors = sectors,
      drop_dB = drop_dB,
      include_grid = include_grid
    ))
  }
  # Return the collected TMM products ==========================================
  out
}
