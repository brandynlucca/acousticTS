################################################################################
# Transition matrix method (TMM) spheroidal-coordinate branch helpers
################################################################################

# Convert the public/world spherical-angle convention used by the shared TMM
# interface into the internal prolate spheroidal convention whose polar axis is
# aligned with the body-fixed x-axis.
#' @noRd
.tmm_public_to_spheroidal_angles <- function(theta, phi) {
  x <- sin(theta) * cos(phi)
  y <- sin(theta) * sin(phi)
  z <- cos(theta)

  theta_internal <- acos(pmax(-1, pmin(1, x)))
  phi_internal <- atan2(z, y)
  phi_internal <- .tmm_wrap_angle_2pi(phi_internal)

  on_axis <- abs(sin(theta_internal)) < 1e-10
  phi_internal[on_axis] <- 0

  list(theta = theta_internal, phi = phi_internal)
}

# Convert the internal prolate spheroidal convention whose polar axis is aligned
# with the body-fixed x-axis back into the shared public/world spherical-angle
# convention used by the rest of the TMM interface.
#' @noRd
.tmm_spheroidal_to_public_angles <- function(theta, phi) {
  x <- cos(theta)
  y <- sin(theta) * cos(phi)
  z <- sin(theta) * sin(phi)

  theta_public <- acos(pmax(-1, pmin(1, z)))
  phi_public <- atan2(y, x)
  phi_public <- .tmm_wrap_angle_2pi(phi_public)

  on_axis <- abs(sin(theta_public)) < 1e-10
  phi_public[on_axis] <- 0

  list(theta = theta_public, phi = phi_public)
}

# Convert the Furusawa-style retained modal sum into the far-field scattering
# amplitude used elsewhere in TMM. The spheroidal backend returns the series
# before the -2i / k prefactor from the published far-field expression.
#' @noRd
.tmm_spheroidal_sum_to_amplitude <- function(f_sum, k_sw) {
  # Apply the far-field prefactor used by the public TMM outputs ===============
  (-2i / k_sw) * f_sum
}

# Route the prolate branch through the exact scalar spheroidal modal-series
# backend. For the current single-target scope, this is the geometry-matched
# T-matrix-equivalent path.
#' @noRd
.tmm_run_spheroidal_branch <- function(object,
                                       acoustics,
                                       medium,
                                       parameters) {
  # Recover the stored prolate body parameters =================================
  body <- acousticTS::extract(object, "model_parameters")$TMM$body
  body_internal <- body
  incident_internal <- .tmm_public_to_spheroidal_angles(
    theta = body$theta_body,
    phi = body$phi_body
  )
  body_internal$theta_body <- incident_internal$theta[[1L]]
  body_internal$phi_body <- incident_internal$phi[[1L]]
  use_simplified_gas <- identical(parameters$boundary, "gas_filled")
  # Penetrable prolate cases are numerically stiffer, so they inherit the same
  # conservative PSMS settings used in the benchmarked modal-series workflow.
  if (!isTRUE(parameters$store_t_matrix)) {
    # Reuse the exact PSMS backend when no retained blocks are needed ==========
    psms_object <- psms_initialize(
      object = object,
      frequency = acoustics$frequency,
      phi_body = body_internal$phi_body,
      boundary = parameters$boundary,
      adaptive = FALSE,
      precision = parameters$precision,
      n_integration = if (is.na(parameters$n_integration)) {
        NULL
      } else {
        parameters$n_integration
      },
      # The full gas-filled spheroidal kernel remains numerically unstable in
      # the benchmark regime, so the prolate TMM gas branch must stay aligned
      # with the stabilized public PSMS pathway.
      simplify_Amn = use_simplified_gas,
      sound_speed_sw = medium$sound_speed,
      density_sw = medium$density
    )
    psms_object <- PSMS(psms_object)
    psms_parameters <- acousticTS::extract(
      psms_object,
      "model_parameters"
    )$PSMS$parameters
    f_bs <- .tmm_spheroidal_sum_to_amplitude(
      psms_object@model$PSMS$f_bs,
      acoustics$k_sw
    )
    sigma_bs <- .sigma_bs(f_bs)

    return(
      list(
        model = data.frame(
          frequency = acoustics$frequency,
          f_bs = f_bs,
          sigma_bs = sigma_bs,
          TS = db(sigma_bs),
          n_max = psms_parameters$acoustics$n_max
        ),
        n_max = psms_parameters$acoustics$n_max,
        t_matrix = NULL
      )
    )
  }
  # Build the quadrature rule for the retained spheroidal solve ================
  quad_pts <- gauss_legendre(
    n = parameters$n_integration %||% 96L,
    a = -1,
    b = 1
  )
  # Evaluate the retained spheroidal T-matrix blocks ===========================
  tmm_detail <- prolate_spheroid_tmatrix_cpp(
    acoustics = acoustics,
    body = body_internal,
    medium = medium,
    integration_pts = quad_pts,
    precision = parameters$precision,
    Amn_method = switch(parameters$boundary,
      liquid_filled = "Amn_fluid",
      gas_filled = "Amn_fluid_simplify",
      fixed_rigid = "Amn_fixed_rigid",
      pressure_release = "Amn_pressure_release"
    )
  )
  # Convert the retained modal sums into the public TMM outputs ================
  f_bs <- .tmm_spheroidal_sum_to_amplitude(
    tmm_detail$f_scat,
    acoustics$k_sw
  )
  sigma_bs <- .sigma_bs(f_bs)
  # Return the monostatic spectrum and stored spheroidal blocks ================
  list(
    model = data.frame(
      frequency = acoustics$frequency,
      f_bs = f_bs,
      sigma_bs = sigma_bs,
      TS = db(sigma_bs),
      n_max = acoustics$n_max
    ),
    n_max = acoustics$n_max,
    t_matrix = tmm_detail$t_matrix
  )
}
