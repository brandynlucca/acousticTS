################################################################################
# Transition matrix method (TMM) for single-target axisymmetric scatterers
################################################################################
#' Single-target transition matrix method (TMM)
#'
#' @description
#' Computes monostatic backscatter from a single axisymmetric target using a
#' transition-matrix formulation. The current implementation targets smooth
#' bodies of revolution and finite cylinders already represented in the
#' package as a `Sphere`, `OblateSpheroid`, `ProlateSpheroid`, or `Cylinder`.
#' It also supports spherical fluid shells and spherical elastic
#' shells carried by `ESS` objects. An explicit experimental hybrid-reference
#' path is also available for elastic-shelled prolate spheroids under axial
#' incidence. The current public boundaries therefore cover rigid,
#' pressure-release, and homogeneous penetrable fluid/gas interiors for
#' homogeneous bodies, plus `shelled_pressure_release`, `shelled_liquid`, and
#' `shelled_gas` for shell spheres, and `elastic_shelled` for spherical elastic
#' shells plus the experimental axial prolate-shell branch.
#'
#' @details
#' This implementation is intentionally scoped to **single targets** and the
#' monostatic backscatter quantity used by `target_strength()`.
#'
#' For spheres and oblate spheroids, the current implementation uses a
#' spherical-wave basis with an explicit projected boundary solve over the
#' target surface. Spherical fluid shells instead stay on the exact layered
#' sphere modal path already represented by `SPHMS`, but are retained and
#' post-processed through the same stored `TMM` workflow. Spherical elastic
#' shells analogously stay on the exact `ESSMS` modal path while exposing the
#' same retained `TMM` angular post-processing interface. For prolate
#' spheroids, it instead uses a
#' spheroidal-coordinate transition-matrix-equivalent backend, which is the
#' more natural coordinate system for that geometry and is consistent with the
#' scalar spheroidal T-matrix literature for single-target scattering. For
#' finite cylinders, the default monostatic branch uses a
#' cylindrical-coordinate modal/T-matrix-equivalent backend so that the
#' backscatter benchmark remains aligned with the exact finite-cylinder family.
#' When \code{store_t_matrix = TRUE}, sharp cylinders keep that exact
#' finite-cylinder monostatic family on the stored model table and retain only
#' the exact monostatic reuse path needed by the standard monostatic
#' post-processing helpers. Public general-angle cylinder bistatic helpers are
#' intentionally outside the current scope. Cylinder calls therefore still emit
#' a warning by default while this narrower validation scope is being
#' maintained; see \code{options(acousticTS.warn_tmm_cylinder = FALSE)} to
#' silence it in controlled test or benchmarking workflows. Elastic-shelled
#' prolate spheroids currently expose an experimental retained-grid route backed
#' by the coupled shell-fluid hybrid solver in the sibling
#' `acousticTSValidation` repository. That branch currently supports axial
#' incidence only (`theta_body_deg = 0` or `180`) and is enabled explicitly
#' with `elastic_shell_backend = "hybrid_reference"`.
#'
#' The present solver is therefore a practical single-target acoustic
#' T-matrix method motivated by the classic transition-matrix literature, but
#' it is **not yet** a full implementation of the far-field two-part
#' T-matrix workflow of Ganesh and Hawkins (2008, 2022), which assumes an
#' external far-field solver as the first stage.
#'
#' @section Usage:
#' This model is accessed via:
#' \preformatted{
#' target_strength(
#'   ...,
#'   model = "TMM",
#'   boundary,
#'   sound_speed_sw,
#'   density_sw,
#'   n_max,
#'   cylinder_backend,
#'   store_t_matrix,
#'   cylinder_endcap_fraction,
#'   elastic_shell_backend,
#'   theta_body_deg,
#'   validation_repo
#' )
#' }
#'
#' @section Arguments:
#' \describe{
#'   \item{\code{boundary}}{Boundary condition at the target surface. One of
#'   \code{"fixed_rigid"}, \code{"pressure_release"},
#'   \code{"liquid_filled"}, \code{"gas_filled"},
#'   \code{"shelled_pressure_release"}, \code{"shelled_liquid"}, or
#'   \code{"shelled_gas"}, or \code{"elastic_shelled"}. The fluid-shell
#'   boundaries are currently restricted to spherical `ESS` objects. The
#'   `elastic_shelled` boundary supports spherical `ESS` objects directly and
#'   also exposes an experimental axial-incidence prolate-shell route when
#'   `elastic_shell_backend = "hybrid_reference"` is supplied.}
#'   \item{\code{sound_speed_sw}}{Surrounding-medium sound speed
#'   (\eqn{m~s^{-1}}).}
#'   \item{\code{density_sw}}{Surrounding-medium density (\eqn{kg~m^{-3}}).}
#'   \item{\code{n_max}}{Optional truncation limit. For spheres and oblate
#'   spheroids, this is the maximum spherical-wave degree used in the truncated
#'   T-matrix solve. For the default monostatic cylinder branch, it is the
#'   cylindrical modal cutoff used in the geometry-matched backend. When left
#'   as \code{NULL}, a geometry-aware rule is used frequency-by-frequency. This
#'   argument is currently ignored for prolate spheroids, which use the
#'   spheroidal-coordinate branch.}
#'   \item{\code{cylinder_backend}}{Cylinder-specific backend selection. One of
#'   \code{"auto"}, \code{"legacy"}, or \code{"retained"}. The default
#'   \code{"auto"} keeps plain sharp-cylinder monostatic runs on the
#'   geometry-matched cylindrical branch, but can switch stored-cylinder
#'   workflows and tapered or explicitly endcap-smoothed cylinders onto the
#'   retained internal branch used for exact monostatic reuse. This argument is
#'   ignored for non-cylinder shapes.}
#'   \item{\code{cylinder_endcap_fraction}}{Optional non-negative smoothing
#'   fraction used only by the retained axisymmetric cylinder backend. When
#'   positive, this rounds each sharp cylinder endcap over that fraction of the
#'   half-length using a short spheroidal cap profile. The default \code{NULL}
#'   keeps sharp cylinders sharp. This argument is ignored for non-cylinder
#'   shapes and for the default legacy cylinder backend.}
#'   \item{\code{store_t_matrix}}{Logical flag controlling whether the
#'   frequency-specific retained state is stored under
#'   \code{object@model_parameters$TMM$parameters$t_matrix}. The default is
#'   \code{FALSE} to avoid large object sizes. Explicit block retention is
#'   available for the spherical and spheroidal branches. For cylinders, sharp
#'   stored runs keep the geometry-matched monostatic family on the model table
#'   while retaining only the exact monostatic reuse path; public cylinder
#'   bistatic and grid helpers remain outside scope. For the experimental
#'   elastic-shelled prolate branch, the retained state is the externally
#'   generated hybrid scattering grid.}
#'   \item{\code{elastic_shell_backend}}{Elastic-shelled prolate backend
#'   selection. The default \code{"default"} keeps that branch disabled.
#'   Supply \code{"hybrid_reference"} to use the coupled shell-fluid hybrid
#'   solver from the sibling \code{acousticTSValidation} repository. Ignored
#'   for all non-prolate or non-elastic-shell calls.}
#'   \item{\code{theta_body_deg}}{Optional body-axis incidence angle in
#'   degrees for the experimental elastic-shelled prolate branch. Only
#'   \code{0} and \code{180} are supported at present. Ignored otherwise.}
#'   \item{\code{validation_repo}}{Optional path to the sibling
#'   \code{acousticTSValidation} repository used by the experimental
#'   elastic-shelled prolate hybrid backend.}
#' }
#'
#' @section Theory:
#' For a single target, the incident and scattered fields are expanded in
#' regular and outgoing modal bases, respectively:
#'
#' \deqn{
#' p^{inc} = \sum_{\nu} a_{\nu} \, \psi_{\nu}^{(1)}, \qquad
#' p^{sca} = \sum_{\nu} f_{\nu} \, \psi_{\nu}^{(3)},
#' }
#'
#' where the transition matrix \eqn{\mathbf{T}} maps incident coefficients to
#' scattered coefficients:
#'
#' \deqn{
#' \mathbf{f} = \mathbf{T}\mathbf{a}.
#' }
#'
#' For the axisymmetric single-target case used here, the azimuthal orders
#' decouple and each block is recovered by enforcing the boundary conditions on
#' the target surface for the retained modal basis functions. The
#' backscatter amplitude is then obtained by evaluating the outgoing expansion
#' in the monostatic receive direction opposite to the incident plane wave.
#'
#' @seealso
#' \code{\link{target_strength}}, \code{\link{FLS}}, \code{\link{GAS}},
#' \code{\link{ESS}}, \code{\link{Sphere}}, \code{\link{OblateSpheroid}},
#' \code{\link{ProlateSpheroid}}, \code{\link{Cylinder}},
#' \code{\link{sphere}}, \code{\link{oblate_spheroid}},
#' \code{\link{prolate_spheroid}}, \code{\link{cylinder}}
#'
#' @references
#'
#' Waterman, P. C. (1969). New formulation of acoustic scattering. *The Journal
#' of the Acoustical Society of America*, **45**, 1417-1429.
#'
#' Varadan, V. K., Varadan, V. V., Bringi, V. N., and Waterman, P. C. (1982).
#' Computation of rigid body scattering by prolate spheroids using the T-matrix
#' approach. *The Journal of the Acoustical Society of America*, **71**, 22-25.
#'
#' Hackman, R. H. (1984). An application of the spheroidal-coordinate-based
#' transition matrix: The acoustic scattering from high aspect ratio solids.
#' *The Journal of the Acoustical Society of America*, **76**, 1058-1070.

#' Waterman, P. C. (2009). T-matrix methods in acoustic scattering. *The
#' Journal of the Acoustical Society of America*, **125**, 42-51.
#'
#' Ganesh, M., and Hawkins, S. C. (2008). A far-field based T-matrix method for
#' three dimensional acoustic scattering. *Wave Motion*, **45**, 1441-1460.
#'
#' Ganesh, M., and Hawkins, S. C. (2022). A numerically stable T-matrix method
#' for acoustic scattering by nonspherical particles with large aspect ratios
#' and size parameters. *The Journal of the Acoustical Society of America*,
#' **151**, 1978-1988.
#'
#' @name TMM
#' @aliases tmm TMM
#' @docType data
#' @keywords models acoustics internal
NULL

#' Initialize object for the transition matrix method
#'
#' @param object Scatterer-class object.
#' @param frequency Frequency vector (Hz).
#' @param boundary Boundary condition at the target surface.
#' @param sound_speed_sw Surrounding-medium sound speed (m/s).
#' @param density_sw Surrounding-medium density (kg/m^3).
#' @param n_max Optional truncation degree for the spherical-wave basis.
#' @param cylinder_backend Optional cylinder backend selector.
#' @param cylinder_endcap_fraction Optional smoothing fraction for the retained
#'   axisymmetric cylinder branch.
#' @param store_t_matrix Logical; retain the frequency-specific T-matrix blocks.
#' @keywords internal
#' @noRd
tmm_initialize <- function(object,
                           frequency,
                           boundary = NULL,
                           sound_speed_sw = .SEAWATER_SOUND_SPEED_DEFAULT,
                           density_sw = .SEAWATER_DENSITY_DEFAULT,
                           n_max = NULL,
                           cylinder_backend = "auto",
                           cylinder_endcap_fraction = NULL,
                           store_t_matrix = FALSE,
                           elastic_shell_backend = c("default", "hybrid_reference"),
                           validation_repo = NULL,
                           python = Sys.getenv("ACOUSTICTS_PYTHON", unset = "python"),
                           theta_body_deg = NULL,
                           n_eta = 129L,
                           mesh_n_u = 12L,
                           mesh_n_v = 16L,
                           n_theta = 91L,
                           n_phi = 181L) {
  # Enforce the current homogeneous-fluid scatterer scope ======================
  .tmm_validate_object_scope(object)
  .tmm_validate_store_t_matrix(store_t_matrix)
  elastic_shell_backend <- match.arg(elastic_shell_backend)
  cylinder_endcap_fraction <- .tmm_resolve_cylinder_endcap_fraction(
    cylinder_endcap_fraction
  )
  boundary <- .tmm_resolve_boundary(object, boundary)
  shape_parameters <- acousticTS::extract(object, "shape_parameters")
  use_elastic_shell_prolate_branch <- .tmm_is_elastic_shell_prolate_branch(
    object = object,
    shape_parameters = shape_parameters,
    boundary = boundary
  )
  if (methods::is(object, "ESS") &&
    !.tmm_is_sphere_modal_branch(
      object = object,
      shape_parameters = shape_parameters,
      boundary = boundary
    ) &&
    !use_elastic_shell_prolate_branch) {
    stop(
      "The current TMM shell branch supports spherical fluid shells plus ",
      "spherical elastic shells, along with the experimental axial-incidence ",
      "elastic-shelled prolate branch. Supported ESS boundaries are ",
      "'shelled_pressure_release', 'shelled_liquid', 'shelled_gas', and ",
      "'elastic_shelled'.",
      call. = FALSE
    )
  }
  branch_flags <- .tmm_branch_flags(
    shape_parameters = shape_parameters,
    boundary = boundary,
    store_t_matrix = store_t_matrix,
    cylinder_backend = cylinder_backend,
    cylinder_endcap_fraction = cylinder_endcap_fraction
  )
  use_spheroidal_branch <- branch_flags$use_spheroidal_branch
  use_cylindrical_branch <- branch_flags$use_cylindrical_branch
  use_shell_sphere_branch <- .tmm_is_sphere_modal_branch(
    object = object,
    shape_parameters = shape_parameters,
    boundary = boundary
  )
  resolved_cylinder_backend <- branch_flags$cylinder_backend

  # Resolve the boundary condition and validate the storage controls ===========
  n_max <- .tmm_branch_n_max(n_max, use_spheroidal_branch)
  body <- .tmm_prepare_body(object, sound_speed_sw, density_sw, boundary)
  if (use_elastic_shell_prolate_branch) {
    if (!identical(elastic_shell_backend, "hybrid_reference")) {
      stop(
        paste(
          "Elastic-shelled prolate TMM currently requires",
          "'elastic_shell_backend = \"hybrid_reference\"'."
        ),
        call. = FALSE
      )
    }
    if (is.null(theta_body_deg)) {
      stop(
        paste(
          "Elastic-shelled prolate TMM currently needs an explicit axial",
          "'theta_body_deg' override of 0 or 180."
        ),
        call. = FALSE
      )
    }
    theta_body_deg <- as.numeric(theta_body_deg)[1]
    if (!is.finite(theta_body_deg) ||
      !(isTRUE(all.equal(theta_body_deg %% 180, 0, tolerance = 1e-8)))) {
      stop(
        "Elastic-shelled prolate TMM currently supports axial 'theta_body_deg' only (0 or 180).",
        call. = FALSE
      )
    }
    body$theta <- pi / 2
    body$phi_body <- if ((theta_body_deg %% 360) < 90 || (theta_body_deg %% 360) > 270) 0 else pi
  }

  # Build the shared acoustics table for the requested frequencies =============
  acoustics_info <- .tmm_prepare_acoustics(
    frequency = frequency,
    sound_speed_sw = sound_speed_sw,
    body = body,
    boundary = boundary,
    shape_parameters = shape_parameters,
    use_spheroidal_branch = use_spheroidal_branch,
    use_cylindrical_branch = use_cylindrical_branch,
    n_max = n_max,
    cylinder_backend = resolved_cylinder_backend,
    cylinder_endcap_fraction = cylinder_endcap_fraction,
    store_t_matrix = store_t_matrix
  )
  acoustics <- acoustics_info$acoustics
  geometry <- acoustics_info$geometry

  # Assemble the stored body/geometry metadata for downstream reuse ============
  body_params <- .tmm_body_parameters(body, geometry)

  # Initialize the TMM model slots and optional retained-block storage =========
  .init_model_slots(
    object = object,
    model_name = "TMM",
    frequency = frequency,
    model_parameters = list(
      parameters = list(
        acoustics = acoustics,
        boundary = boundary,
        cylinder_backend = resolved_cylinder_backend,
        coordinate_system = .tmm_coordinate_system(
          use_spheroidal_branch,
          use_cylindrical_branch,
          use_shell_sphere_branch = use_shell_sphere_branch,
          use_elastic_shell_prolate_branch = use_elastic_shell_prolate_branch,
          shape_parameters = shape_parameters,
          cylinder_backend = resolved_cylinder_backend
        ),
        precision = .tmm_precision_label(use_spheroidal_branch, boundary),
        n_integration = .tmm_n_integration_label(
          use_spheroidal_branch,
          boundary
        ),
        elastic_shell_backend = elastic_shell_backend,
        validation_repo = validation_repo,
        python = python,
        theta_body_deg = theta_body_deg,
        n_eta = as.integer(n_eta),
        mesh_n_u = as.integer(mesh_n_u),
        mesh_n_v = as.integer(mesh_n_v),
        n_theta = as.integer(n_theta),
        n_phi = as.integer(n_phi),
        cylinder_endcap_fraction = cylinder_endcap_fraction,
        store_t_matrix = store_t_matrix,
        t_matrix = if (isTRUE(store_t_matrix)) {
          vector("list", length(frequency))
        } else {
          NULL
        }
      ),
      medium = .init_medium_params(sound_speed_sw, density_sw),
      body = body_params
    ),
    result_cols = c("f_bs", "sigma_bs", "TS", "n_max")
  )
}

# Build the stored hybrid-grid branch used for the current elastic-shelled
# prolate TMM rebuild route.
#' @noRd
.tmm_run_elastic_shell_prolate_hybrid_branch <- function(acoustics,
                                                         body,
                                                         parameters) {
  theta_body_deg <- parameters$theta_body_deg %||% 0
  t_store <- if (isTRUE(parameters$store_t_matrix)) {
    lapply(
      acoustics$frequency,
      function(freq_i) {
        .espsms_hybrid_reference_grid_from_body(
          body = body,
          frequency = freq_i,
          theta_body_deg = theta_body_deg,
          n_eta = parameters$n_eta %||% 129L,
          mesh_n_u = parameters$mesh_n_u %||% 12L,
          mesh_n_v = parameters$mesh_n_v %||% 16L,
          n_theta = parameters$n_theta %||% 91L,
          n_phi = parameters$n_phi %||% 181L,
          validation_repo = parameters$validation_repo,
          python = parameters$python %||% Sys.getenv("ACOUSTICTS_PYTHON", unset = "python")
        )
      }
    )
  } else {
    NULL
  }

  if (isTRUE(parameters$store_t_matrix)) {
    monostatic_index <- function(store_i) {
      theta_target <- pi - store_i$theta_body
      phi_target <- .tmm_wrap_angle_2pi(store_i$phi_body + pi)
      theta_idx <- which.min(abs(store_i$theta_scatter - theta_target))
      phi_idx <- which.min(abs(store_i$phi_scatter - phi_target))
      c(theta_idx = theta_idx, phi_idx = phi_idx)
    }
    f_bs <- vapply(
      t_store,
      function(store_i) {
        idx <- monostatic_index(store_i)
        store_i$f_scat[idx[[1]], idx[[2]]]
      },
      complex(1)
    )
  } else {
    mono <- .espsms_hybrid_reference_monostatic_from_body(
      body = body,
      frequency = acoustics$frequency,
      theta_body_deg = theta_body_deg,
      n_eta = parameters$n_eta %||% 129L,
      mesh_n_u = parameters$mesh_n_u %||% 12L,
      mesh_n_v = parameters$mesh_n_v %||% 16L,
      validation_repo = parameters$validation_repo,
      python = parameters$python %||% Sys.getenv("ACOUSTICTS_PYTHON", unset = "python")
    )
    f_bs <- complex(real = mono$f_bs_real, imaginary = mono$f_bs_imag)
  }

  sigma_bs <- .sigma_bs(f_bs)
  list(
    model = data.frame(
      frequency = acoustics$frequency,
      f_bs = f_bs,
      sigma_bs = sigma_bs,
      TS = db(sigma_bs),
      n_max = acoustics$n_max
    ),
    t_matrix = t_store
  )
}

# Build the exact spherical shell modal coefficients retained by the TMM shell
# branch so bistatic scattering can be reconstructed directly from the
# layered-fluid modal series.
#' @noRd
.tmm_store_sphere_modal_branch <- function(acoustics, body, boundary) {
  Am <- .sphms_modal_coefficients(
    k1a = acoustics$k_sw * body$radius_shell,
    k2a = acoustics$k_shell * body$radius_shell,
    k2b = acoustics$k_shell * body$radius_fluid,
    k3a = acoustics$k_fluid * body$radius_fluid,
    k3b = acoustics$k_fluid * body$radius_fluid,
    g21 = body$g21,
    g31 = body$g31,
    g32 = body$g32,
    h21 = body$h21,
    h31 = body$h31,
    h32 = body$h32,
    m_limit = acoustics$n_max,
    Bm_method = .sphms_Bm_method(boundary)
  )

  lapply(
    seq_len(nrow(acoustics)),
    function(i) {
      n_seq <- 0:as.integer(acoustics$n_max[i])
      list(
        family = "sphere_modal",
        n_seq = n_seq,
        A_n = as.vector(Am[seq_along(n_seq), i])
      )
    }
  )
}

# Store the exact elastic-shell sphere modal coefficients in the same retained
# spherical-modal structure used by the shell-sphere fluid branch.
#' @noRd
.tmm_store_elastic_shell_modal_branch <- function(acoustics, body) {
  sound_speed_longitudinal <- sqrt(
    (body$shell_lambda + 2 * body$shell_G) / body$shell_density
  )
  sound_speed_transversal <- sqrt(body$shell_G / body$shell_density)
  ka_matrix <- .calculate_ka_matrix(
    frequency = acoustics$frequency,
    sound_speed_sw = body$medium_sound_speed,
    sound_speed_fluid = body$fluid_sound_speed,
    sound_speed_longitudinal = sound_speed_longitudinal,
    sound_speed_transversal = sound_speed_transversal,
    radius_shell = body$radius_shell,
    radius_fluid = body$radius_fluid
  )
  Am <- elastic_shell_boundary_conditions(
    ka_matrix = ka_matrix,
    m_limit = acoustics$m_limit,
    lambda = body$shell_lambda,
    mu = body$shell_G,
    rho_ratio_sw = body$medium_density / body$shell_density,
    rho_ratio_fl = body$fluid_density / body$shell_density
  )

  lapply(
    seq_len(nrow(acoustics)),
    function(i) {
      n_seq <- 0:as.integer(acoustics$m_limit[i])
      list(
        family = "sphere_modal",
        n_seq = n_seq,
        A_n = as.vector(Am[i, seq_along(n_seq)])
      )
    }
  )
}

# Evaluate the exact stored spherical-modal TMM branches used by shell spheres.
#' @noRd
.tmm_run_shell_sphere_branch <- function(acoustics, body, boundary) {
  t_store <- if (identical(boundary, "elastic_shelled")) {
    .tmm_store_elastic_shell_modal_branch(
      acoustics = acoustics,
      body = body
    )
  } else {
    .tmm_store_sphere_modal_branch(
      acoustics = acoustics,
      body = body,
      boundary = boundary
    )
  }
  sph_bm <- vapply(
    seq_along(t_store),
    function(i) {
      sum(
        .sphms_modal_weights(acoustics$n_max[i]) * t_store[[i]]$A_n,
        na.rm = TRUE
      )
    },
    complex(1)
  )
  f_bs <- -1i / acoustics$k_sw * sph_bm
  sigma_bs <- .sigma_bs(f_bs)

  list(
    model = data.frame(
      frequency = acoustics$frequency,
      f_bs = f_bs,
      sigma_bs = sigma_bs,
      TS = db(sigma_bs),
      n_max = acoustics$n_max
    ),
    t_matrix = t_store
  )
}


#' Single-target transition matrix method (TMM)
#'
#' @param object Scatterer-object initialized for TMM.
#' @noRd
TMM <- function(object) {
  # Recover the initialization recipe prepared by `tmm_initialize()` ===========
  model_params <- acousticTS::extract(object, "model_parameters")$TMM
  parameters <- model_params$parameters
  acoustics <- parameters$acoustics
  medium <- model_params$medium
  body <- model_params$body
  shape_parameters <- acousticTS::extract(object, "shape_parameters")

  if (identical(parameters$coordinate_system, "sphere_modal")) {
    shell_result <- .tmm_run_shell_sphere_branch(
      acoustics = acoustics,
      body = body,
      boundary = parameters$boundary
    )
    methods::slot(object, "model")$TMM <- shell_result$model
    if (isTRUE(parameters$store_t_matrix)) {
      methods::slot(object, "model_parameters")$TMM$parameters$t_matrix <-
        shell_result$t_matrix
    }

    return(object)
  }

  if (identical(parameters$coordinate_system, "espsms_hybrid_grid")) {
    hybrid_result <- .tmm_run_elastic_shell_prolate_hybrid_branch(
      acoustics = acoustics,
      body = body,
      parameters = parameters
    )
    methods::slot(object, "model")$TMM <- hybrid_result$model
    if (isTRUE(parameters$store_t_matrix)) {
      methods::slot(object, "model_parameters")$TMM$parameters$t_matrix <- hybrid_result$t_matrix
    }
    return(object)
  }

  # Route prolate targets through the geometry-matched spheroidal backend ======
  if (.tmm_is_spheroidal_branch(shape_parameters, parameters$boundary)) {
    spheroidal_result <- .tmm_run_spheroidal_branch(
      object = object,
      acoustics = acoustics,
      medium = medium,
      parameters = parameters
    )

    methods::slot(object, "model")$TMM <- spheroidal_result$model
    methods::slot(object, "model_parameters")$TMM$parameters$acoustics$n_max <-
      spheroidal_result$n_max
    if (isTRUE(parameters$store_t_matrix)) {
      methods::slot(object, "model_parameters")$TMM$parameters$t_matrix <-
        spheroidal_result$t_matrix
    }

    return(object)
  }

  # Route cylinders through the finite-cylinder exact-family bookkeeping =======
  if (identical(parameters$coordinate_system, "cylindrical")) {
    cylindrical_result <- .tmm_run_cylindrical_branch(
      shape_parameters = shape_parameters,
      acoustics = acoustics,
      body = body,
      boundary = parameters$boundary
    )
    methods::slot(object, "model")$TMM <- cylindrical_result$model
    if (isTRUE(parameters$store_t_matrix)) {
      methods::slot(object, "model_parameters")$TMM$parameters$t_matrix <-
        .tmm_store_cylindrical_branch(acoustics)
    }
    return(object)
  }

  # Stored cylinder workflows now use the dedicated cylinder-native retained
  # state rather than the shared spherical retained operator. Sharp cylinders
  # still report the exact finite-cylinder monostatic family on the model
  # table, while tapered or explicitly smoothed cylinders use the native
  # retained evaluator for their monostatic outputs too.
  if (identical(parameters$coordinate_system, "cylinder_native")) {
    is_tapered <- "taper_order" %in% names(shape_parameters) &&
      is.finite(as.numeric(shape_parameters$taper_order)[1])
    is_sharp_cylinder <- !is_tapered &&
      (.tmm_resolve_cylinder_endcap_fraction(
        parameters$cylinder_endcap_fraction
      ) %||% 0) <= 0

    if (isTRUE(parameters$store_t_matrix)) {
      methods::slot(object, "model_parameters")$TMM$parameters$t_matrix <-
        .tmm_store_cylindrical_branch(
          acoustics = acoustics,
          family = "cylinder_native"
        )
    }

    if (is_sharp_cylinder) {
      exact_cyl <- .tmm_run_cylindrical_branch(
        shape_parameters = shape_parameters,
        acoustics = acoustics,
        body = body,
        boundary = parameters$boundary
      )
      f_bs_model <- exact_cyl$model$f_bs
    } else {
      native_params <- model_params
      native_params$parameters$t_matrix <- .tmm_store_cylindrical_branch(
        acoustics = acoustics,
        family = "cylinder_native"
      )
      f_bs_model <- .tmm_scattering_cylinder_native(
        model_params = native_params,
        shape_parameters = shape_parameters,
        theta_body = body$theta_body,
        phi_body = body$phi_body %||% pi,
        theta_scatter = pi - body$theta_body,
        phi_scatter = (body$phi_body %||% pi) + pi
      )
    }

    sigma_bs <- .sigma_bs(f_bs_model)
    methods::slot(object, "model")$TMM <- data.frame(
      frequency = acoustics$frequency,
      f_bs = f_bs_model,
      sigma_bs = sigma_bs,
      TS = db(sigma_bs),
      n_max = acoustics$n_max
    )

    return(object)
  }

  # Use the fast compiled spherical branch when block storage is disabled ======
  use_compiled_spherical <- !(
    identical(.tmm_shape_name(shape_parameters), "Cylinder") &&
      identical(parameters$cylinder_backend, "retained")
  ) && !identical(parameters$boundary, "elastic_shelled")

  if (!isTRUE(parameters$store_t_matrix) && use_compiled_spherical) {
    f_bs <- tmm_backscatter_cpp(
      frequency = acoustics$frequency,
      theta_body = body$theta_body,
      shape = .tmm_shape_name(shape_parameters),
      shape_values = .tmm_shape_values(shape_parameters),
      boundary = parameters$boundary,
      sound_speed_sw = medium$sound_speed,
      density_sw = medium$density,
      density_body = body$density %||% NA_real_,
      sound_speed_body = body$sound_speed %||% NA_real_,
      n_max = acoustics$n_max
    )
  } else {
    # Retain the per-frequency spherical blocks for post-processing ============
    f_bs <- complex(length.out = nrow(acoustics))
    t_store <- vector("list", nrow(acoustics))

    for (i in seq_len(nrow(acoustics))) {
      tmm_i <- .tmm_single_frequency_spherical(
        k_sw = acoustics$k_sw[i],
        k_body = acoustics$k_body[i],
        theta_body = body$theta_body,
        boundary = parameters$boundary,
        shape_parameters = shape_parameters,
        rho_sw = medium$density,
        rho_body = body$density,
        n_max = acoustics$n_max[i],
        cylinder_endcap_fraction = parameters$cylinder_endcap_fraction,
        store_t_matrix = parameters$store_t_matrix,
        shell_body = if (identical(parameters$boundary, "elastic_shelled")) {
          body
        } else {
          NULL
        },
        frequency_hz = acoustics$frequency[i]
      )

      f_bs[i] <- tmm_i$f_bs
      t_store[[i]] <- tmm_i$blocks
    }

    methods::slot(object, "model_parameters")$TMM$parameters$t_matrix <- t_store
  }

  # Store the monostatic backscatter outputs on the scatterer object ===========
  sigma_bs <- .sigma_bs(f_bs)
  methods::slot(object, "model")$TMM <- data.frame(
    frequency = acoustics$frequency,
    f_bs = f_bs,
    sigma_bs = sigma_bs,
    TS = db(sigma_bs),
    n_max = acoustics$n_max
  )

  object
}
