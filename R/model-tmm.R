################################################################################
# Transition matrix method (TMM) for single-target axisymmetric scatterers
################################################################################
#' Single-target transition matrix method (TMM)
#'
#' @description
#' Computes monostatic backscatter from a single axisymmetric target using a
#' transition-matrix formulation. The current implementation targets smooth
#' bodies of revolution and finite cylinders already represented in the
#' package as a `Sphere`, `OblateSpheroid`, `ProlateSpheroid`, or `Cylinder`,
#' and supports rigid, pressure-release, and homogeneous penetrable fluid/gas
#' interiors.
#'
#' @details
#' This implementation is intentionally scoped to **single targets**. The
#' default `target_strength()` result is monostatic backscatter, and retained
#' T-matrix state can also be reused for angular post-processing where the
#' active branch supports it.
#'
#' For spheres and oblate spheroids, the current implementation uses a
#' spherical-wave T-matrix: each retained azimuthal order is solved as an
#' incident-to-outgoing modal coefficient map after projecting the boundary
#' conditions on the target surface. For prolate spheroids, the implementation
#' uses a geometry-matched spheroidal-basis T-matrix formulation, which is the
#' natural coordinate system for that geometry and is consistent with the
#' scalar spheroidal transition-matrix literature for single-target scattering.
#' For finite cylinders, the default monostatic branch uses a
#' cylindrical-coordinate modal T-matrix-style backend so that the backscatter
#' benchmark remains aligned with the exact finite-cylinder family.
#' When \code{store_t_matrix = TRUE}, cylinders retain lightweight
#' cylindrical-family state that supports exact monostatic reuse and
#' orientation-averaged monostatic products, while general-angle cylinder
#' bistatic post-processing remains outside the current validated scope. Because
#' of that narrower validation status, cylinder calls emit a warning by default;
#' see \code{options(acousticTS.warn_tmm_cylinder = FALSE)} to silence it in
#' controlled test or benchmarking workflows.
#'
#' The sphere, oblate, prolate, and shell-sphere branches are therefore
#' single-target acoustic T-matrix methods in the modal coefficient-map sense:
#' they represent the target response as a map from incident modal amplitudes to
#' scattered modal amplitudes and reuse that retained state for supported
#' post-processing. The cylinder branch is more limited and should be treated as
#' a guarded monostatic modal T-matrix-style branch until a validated
#' general-angle cylindrical operator is available.
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
#'   store_t_matrix
#' )
#' }
#'
#' @section Arguments:
#' \describe{
#'   \item{\code{boundary}}{Boundary condition at the target surface. One of
#'   \code{"fixed_rigid"}, \code{"pressure_release"},
#'   \code{"liquid_filled"}, or \code{"gas_filled"}.}
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
#'   \item{\code{store_t_matrix}}{Logical flag controlling whether the
#'   frequency-specific retained state is stored under
#'   \code{object@model_parameters$TMM$parameters$t_matrix}. The default is
#'   \code{FALSE} to avoid large object sizes. Explicit block retention is
#'   available for the spherical and spheroidal branches. For cylinders, the
#'   stored state keeps the geometry-matched cylindrical monostatic family
#'   available for exact monostatic reuse and orientation-averaged monostatic
#'   products; full general-angle cylinder bistatic post-processing is not yet
#'   provided.}
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
#' decouple. Each retained block is recovered in the basis used by the active
#' geometry branch: spherical for spheres and oblates, spheroidal for prolates,
#' and exact spherical modal coefficients for supported shell spheres. The
#' backscatter amplitude is obtained by evaluating the outgoing expansion in the
#' monostatic receive direction opposite to the incident plane wave. When the
#' retained state is stored, the same coefficient map can be reused for
#' supported angular, bistatic, diagnostic, and orientation-averaged products.
#'
#' @seealso
#' \code{\link{target_strength}}, \code{\link{FLS}}, \code{\link{GAS}},
#' \code{\link{Sphere}}, \code{\link{OblateSpheroid}},
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
#' @param store_t_matrix Logical; retain the frequency-specific T-matrix blocks.
#' @keywords internal
#' @noRd
tmm_initialize <- function(object,
                           frequency,
                           boundary = NULL,
                           sound_speed_sw = .SEAWATER_SOUND_SPEED_DEFAULT,
                           density_sw = .SEAWATER_DENSITY_DEFAULT,
                           n_max = NULL,
                           store_t_matrix = FALSE) {
  # Enforce the current homogeneous-fluid scatterer scope ======================
  .tmm_validate_object_scope(object)
  shape_parameters <- acousticTS::extract(object, "shape_parameters")
  boundary <- .tmm_resolve_boundary(object, boundary)
  if (methods::is(object, "ESS") &&
    !.tmm_is_sphere_modal_branch(
      object = object,
      shape_parameters = shape_parameters,
      boundary = boundary
    )) {
    stop(
      "The current TMM shell branch supports spherical fluid shells plus ",
      "spherical elastic shells only. Supported ESS boundaries are ",
      "'shelled_pressure_release', 'shelled_liquid', 'shelled_gas', and ",
      "'elastic_shelled'.",
      call. = FALSE
    )
  }
  branch_flags <- .tmm_branch_flags(shape_parameters, boundary = boundary)
  use_spheroidal_branch <- branch_flags$use_spheroidal_branch
  use_cylindrical_branch <- branch_flags$use_cylindrical_branch
  use_shell_sphere_branch <- .tmm_is_sphere_modal_branch(
    object = object,
    shape_parameters = shape_parameters,
    boundary = boundary
  )

  # Resolve the boundary condition and validate the storage controls ===========
  .tmm_validate_store_t_matrix(store_t_matrix)
  n_max <- .tmm_branch_n_max(n_max, use_spheroidal_branch)
  body <- .tmm_prepare_body(object, sound_speed_sw, density_sw, boundary)

  # Build the shared acoustics table for the requested frequencies =============
  acoustics_info <- .tmm_prepare_acoustics(
    frequency = frequency,
    sound_speed_sw = sound_speed_sw,
    body = body,
    boundary = boundary,
    shape_parameters = shape_parameters,
    use_spheroidal_branch = use_spheroidal_branch,
    use_cylindrical_branch = use_cylindrical_branch,
    n_max = n_max
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
        coordinate_system = .tmm_coordinate_system(
          use_spheroidal_branch,
          use_cylindrical_branch,
          use_shell_sphere_branch = use_shell_sphere_branch
        ),
        precision = .tmm_precision_label(use_spheroidal_branch, boundary),
        n_integration = .tmm_n_integration_label(
          use_spheroidal_branch,
          boundary
        ),
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

  # Use the fast compiled spherical branch when block storage is disabled ======
  if (!isTRUE(parameters$store_t_matrix)) {
    f_bs <- tmm_backscatter_cpp(
      frequency = acoustics$frequency,
      theta_body = body$theta_body,
      shape = shape_parameters$shape,
      shape_values = .tmm_shape_values(shape_parameters),
      boundary = parameters$boundary,
      sound_speed_sw = medium$sound_speed,
      density_sw = medium$density,
      density_body = body$density %||% NA_real_,
      sound_speed_body = body$sound_speed %||% NA_real_,
      n_max = acoustics$n_max
    )
  } else {
    # Retain spherical blocks while sharing setup work within n_max buckets ===
    stored_spherical <- .tmm_spherical_stored_frequency_sweep(
      acoustics = acoustics,
      theta_body = body$theta_body,
      boundary = parameters$boundary,
      shape_parameters = shape_parameters,
      rho_sw = medium$density,
      rho_body = body$density,
      store_t_matrix = parameters$store_t_matrix
    )
    f_bs <- stored_spherical$f_bs
    t_store <- stored_spherical$t_matrix

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
