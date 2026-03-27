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
#' This implementation is intentionally scoped to **single targets** and the
#' monostatic backscatter quantity used by `target_strength()`.
#'
#' For spheres and oblate spheroids, the current implementation uses a
#' spherical-wave basis with an explicit projected boundary solve over the
#' target surface. For prolate spheroids, it instead uses a
#' spheroidal-coordinate transition-matrix-equivalent backend, which is the
#' more natural coordinate system for that geometry and is consistent with the
#' scalar spheroidal T-matrix literature for single-target scattering. For
#' finite cylinders, the default monostatic branch uses a
#' cylindrical-coordinate modal/T-matrix-equivalent backend so that the
#' backscatter benchmark remains aligned with the exact finite-cylinder family.
#' When \code{store_t_matrix = TRUE}, cylinders retain lightweight
#' cylindrical-family state that supports exact monostatic reuse and
#' orientation-averaged monostatic products, while general-angle cylinder
#' bistatic post-processing remains outside the current validated scope. Because
#' of that narrower validation status, cylinder calls emit a warning by default;
#' see \code{options(acousticTS.warn_tmm_cylinder = FALSE)} to silence it in
#' controlled test or benchmarking workflows.
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
#' decouple and each block is recovered by enforcing the boundary conditions on
#' the target surface for the retained modal basis functions. The
#' backscatter amplitude is then obtained by evaluating the outgoing expansion
#' in the monostatic receive direction opposite to the incident plane wave.
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
  branch_flags <- .tmm_branch_flags(shape_parameters)
  use_spheroidal_branch <- branch_flags$use_spheroidal_branch
  use_cylindrical_branch <- branch_flags$use_cylindrical_branch

  # Resolve the boundary condition and validate the storage controls ===========
  boundary <- .tmm_resolve_boundary(object, boundary)
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
          use_cylindrical_branch
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

  # Route prolate targets through the geometry-matched spheroidal backend ======
  if (.tmm_is_spheroidal_branch(shape_parameters)) {
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
        store_t_matrix = parameters$store_t_matrix
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
