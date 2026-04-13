################################################################################
# Elastic-shell prolate spheroidal modal series (ESPSMS)
################################################################################
#' Elastic-shell prolate spheroidal model (ESPSMS)
#'
#' @description
#' The original experimental `ESPSMS` branch used a projected thin-shell
#' impedance approximation on the prolate surface. That route did not validate
#' against the shell-only literature or against independent shell validation and
#' is now disabled while a paper-backed elastic-shelled prolate formulation is
#' rebuilt.
#'
#' @section Usage:
#' By default this model identifier fails fast. During the rebuild, an explicit
#' experimental backend may be requested:
#' \preformatted{
#' target_strength(
#'   ...,
#'   model = "espsms",
#'   sound_speed_sw,
#'   density_sw,
#'   n_max,
#'   solver_backend = "hybrid_reference",
#'   theta_body_deg = 0,
#'   validation_repo = "path/to/acousticTSValidation"
#' )
#' }
#'
#' @section Arguments:
#' \describe{
#'   \item{\code{sound_speed_sw}}{Surrounding-medium sound speed
#'   (\eqn{m~s^{-1}}).}
#'   \item{\code{density_sw}}{Surrounding-medium density (\eqn{kg~m^{-3}}).}
#'   \item{\code{n_max}}{Reserved for the future native shell-only rebuild.}
#'   \item{\code{solver_backend}}{Use \code{"disabled"} to keep the rebuild
#'   branch off, or \code{"hybrid_reference"} to use the external coupled
#'   shell-fluid reference backend.}
#'   \item{\code{theta_body_deg}}{Optional orientation override for the current
#'   experimental backend. Only \code{0} and \code{180} are supported at
#'   present.}
#'   \item{\code{validation_repo}}{Optional path to the sibling
#'   \code{acousticTSValidation} repository when using the hybrid-reference
#'   backend.}
#' }
#'
#' @section Implementation route:
#' The next supported implementation follows the elastic-shelled prolate
#' literature directly, using shell-only modal/numerical equations from the
#' prolate shell vibration/scattering literature together with an independent
#' hybrid shell validator before any `TMM` wrapping is restored.
#'
#' The current rebuild keeps the old proxy branch disabled by default. An
#' explicit experimental backend, `solver_backend = "hybrid_reference"`, is
#' available for monostatic axial-incidence runs when the sibling
#' `acousticTSValidation` repository and the Python hybrid helper are available.
#'
#' @seealso
#' \code{\link{target_strength}}, \code{\link{ESS}},
#' \code{\link{ProlateSpheroid}}, \code{\link{prolate_spheroid}}
#'
#' @references
#'
#' DiMaggio, F. L., and Rand, R. (1966). Axisymmetric vibrations of prolate
#' spheroidal shells. *The Journal of the Acoustical Society of America*,
#' **40**, 179-186.
#'
#' Rand, R., and DiMaggio, F. L. (1967). Vibrations of fluid-filled spherical
#' and spheroidal shells. *The Journal of the Acoustical Society of America*,
#' **42**, 1278-1286.
#'
#' Hayek, S. I., and Boisvert, J. E. (2003). Vibration of prolate spheroidal
#' shells with shear deformation and rotatory inertia: Axisymmetric case.
#' *The Journal of the Acoustical Society of America*, **114**, 2799-2811.
#'
#' @name ESPSMS
#' @aliases espsms ESPSMS epsms EPSMS
#' @docType data
#' @keywords models acoustics internal
NULL

#' Validate the current ESPSMS shape and shell scope.
#' @noRd
.espsms_validate_object <- function(object) {
  if (!methods::is(object, "ESS")) {
    stop(
      "ESPSMS requires an 'ESS' scatterer object.",
      call. = FALSE
    )
  }

  shape_parameters <- acousticTS::extract(object, "shape_parameters")
  shape_name <- as.character(shape_parameters$shape)[1]
  if (!identical(shape_name, "ProlateSpheroid")) {
    stop(
      "ESPSMS currently supports prolate-spheroidal ESS objects only.",
      call. = FALSE
    )
  }

  shell <- acousticTS::extract(object, "shell")
  if (!acousticTS:::.tmm_has_elastic_shell(shell)) {
    stop(
      "ESPSMS requires elastic shell properties on the ESS shell component.",
      call. = FALSE
    )
  }

  fluid <- acousticTS::extract(object, "fluid")
  if (is.null(fluid$density) || is.null(fluid$sound_speed)) {
    stop(
      "ESPSMS requires an explicit interior fluid density and sound speed.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}

#' Build the shared thin-shell prolate body state used by ESPSMS.
#' @noRd
.elastic_shell_prolate_body <- function(object,
                                        sound_speed_sw,
                                        density_sw) {
  .espsms_validate_object(object)

  shape_parameters <- acousticTS::extract(object, "shape_parameters")
  shell_raw <- acousticTS::extract(object, "shell")
  fluid_raw <- acousticTS::extract(object, "fluid")
  shell <- .extract_material_props(shell_raw, sound_speed_sw, density_sw)
  fluid <- .extract_material_props(fluid_raw, sound_speed_sw, density_sw)
  shell[names(.complete_elastic_moduli(
    E = shell$E,
    G = shell$G,
    K = shell$K,
    nu = shell$nu
  ))] <- .complete_elastic_moduli(
    E = shell$E,
    G = shell$G,
    K = shell$K,
    nu = shell$nu
  )

  list(
    theta = shell_raw$theta,
    phi_body = pi,
    density = fluid$density,
    sound_speed = fluid$sound_speed,
    fluid_density = fluid$density,
    fluid_sound_speed = fluid$sound_speed,
    shell_density = shell$density,
    shell_E = shell$E,
    shell_G = shell$G,
    shell_K = shell$K,
    shell_nu = shell$nu,
    shell_lambda = shell$lambda,
    shell_thickness = shell_raw$shell_thickness,
    medium_density = density_sw,
    medium_sound_speed = sound_speed_sw,
    semimajor_length = as.numeric(shape_parameters$shell$semimajor_length)[1],
    semiminor_length = as.numeric(shape_parameters$shell$semiminor_length)[1]
  )
}

#' Build the open meridional eta-grid used by the prolate shell rebuild.
#' @noRd
.espsms_uniform_eta_grid <- function(n_eta = 129L,
                                     pole_offset = 1e-4) {
  if (!is.numeric(n_eta) || length(n_eta) != 1L || !is.finite(n_eta) ||
    n_eta < 5) {
    stop("'n_eta' must be one integer >= 5.", call. = FALSE)
  }
  if (!is.numeric(pole_offset) || length(pole_offset) != 1L ||
    !is.finite(pole_offset) || pole_offset <= 0 || pole_offset >= 0.5) {
    stop("'pole_offset' must lie in (0, 0.5).", call. = FALSE)
  }

  seq(
    from = -1 + pole_offset,
    to = 1 - pole_offset,
    length.out = as.integer(n_eta)
  )
}

#' Build second-order finite-difference matrices on a uniform eta-grid.
#' @noRd
.espsms_fd_matrices <- function(eta) {
  eta <- as.numeric(eta)
  if (!is.numeric(eta) || any(!is.finite(eta)) || length(eta) < 5) {
    stop("'eta' must be a finite numeric vector with length >= 5.", call. = FALSE)
  }

  h <- diff(eta)
  if (!all(abs(h - h[[1]]) <= 1e-12 + 1e-10 * abs(h[[1]]))) {
    stop("Only uniform eta-grids are supported.", call. = FALSE)
  }

  step <- h[[1]]
  n <- length(eta)
  d1 <- matrix(0, nrow = n, ncol = n)
  d2 <- matrix(0, nrow = n, ncol = n)

  for (i in 2:(n - 1)) {
    d1[i, i - 1] <- -0.5 / step
    d1[i, i + 1] <- 0.5 / step
    d2[i, i - 1] <- 1 / (step * step)
    d2[i, i] <- -2 / (step * step)
    d2[i, i + 1] <- 1 / (step * step)
  }

  d1[1, 1:3] <- c(-3, 4, -1) / (2 * step)
  d1[n, (n - 2):n] <- c(1, -4, 3) / (2 * step)
  d2[1, 1:4] <- c(2, -5, 4, -1) / (step * step)
  d2[n, (n - 3):n] <- c(-1, 4, -5, 2) / (step * step)

  list(D1 = d1, D2 = d2)
}

#' Return the prolate shell geometry/state used by the ESPSMS rebuild core.
#' @noRd
.espsms_shell_geometry_state <- function(body,
                                         eta) {
  semimajor_length <- as.numeric(body$semimajor_length)[1]
  semiminor_length <- as.numeric(body$semiminor_length)[1]
  shell_thickness <- as.numeric(body$shell_thickness)[1]
  if (!all(is.finite(c(semimajor_length, semiminor_length, shell_thickness)))) {
    stop(
      "ESPSMS shell geometry requires finite semimajor length, semiminor length, ",
      "and shell thickness.",
      call. = FALSE
    )
  }
  if (semimajor_length <= semiminor_length) {
    stop("ESPSMS rebuild requires a prolate shell with semimajor > semiminor.", call. = FALSE)
  }

  focal_radius <- sqrt(semimajor_length^2 - semiminor_length^2)
  a_shape <- semimajor_length / focal_radius
  interfocal_distance <- 2 * focal_radius
  bending_epsilon <- (shell_thickness / semimajor_length)^2 / 12

  A <- rep(a_shape^2 - 1, length(eta))
  B <- 1 - eta^2
  C <- a_shape^2 - eta^2
  D <- sqrt(A / C)

  list(
    semimajor_length = semimajor_length,
    semiminor_length = semiminor_length,
    shell_thickness = shell_thickness,
    focal_radius = focal_radius,
    interfocal_distance = interfocal_distance,
    aspect_ratio = semimajor_length / semiminor_length,
    a_shape = a_shape,
    bending_epsilon = bending_epsilon,
    eta = eta,
    A = A,
    B = B,
    C = C,
    D = D
  )
}

#' Assemble the axisymmetric nontorsional shell dynamic system for ESPSMS.
#' @noRd
.espsms_axisymmetric_shell_system <- function(body,
                                              frequency_hz,
                                              n_eta = 129L,
                                              pole_offset = 1e-4) {
  eta <- .espsms_uniform_eta_grid(n_eta = n_eta, pole_offset = pole_offset)
  fd <- .espsms_fd_matrices(eta)
  geom <- .espsms_shell_geometry_state(body, eta)

  shell_nu <- as.numeric(body$shell_nu)[1]
  shell_density <- as.numeric(body$shell_density)[1]
  shell_E <- as.numeric(body$shell_E)[1]
  if (!all(is.finite(c(shell_nu, shell_density, shell_E)))) {
    stop(
      "ESPSMS shell system requires finite shell Young's modulus, density, ",
      "and Poisson ratio.",
      call. = FALSE
    )
  }

  cp <- sqrt(shell_E / (shell_density * (1 - shell_nu^2)))
  omega_hat <- 2 * pi * as.numeric(frequency_hz)[1] * geom$semimajor_length / cp

  a <- geom$a_shape
  A <- geom$A
  B <- geom$B
  C <- geom$C
  D <- geom$D
  eps <- geom$bending_epsilon
  F <- (1 - shell_nu) / 2
  gamma <- 3 / 7
  kappa <- pi^2 / 12

  diagm <- function(x) diag(as.numeric(x), nrow = length(x))
  block <- function(d2_coef = NULL, d1_coef = NULL, diag_coef = NULL) {
    n <- length(eta)
    out <- matrix(0, nrow = n, ncol = n)
    if (!is.null(d2_coef)) out <- out + diagm(d2_coef) %*% fd$D2
    if (!is.null(d1_coef)) out <- out + diagm(d1_coef) %*% fd$D1
    if (!is.null(diag_coef)) out <- out + diagm(diag_coef)
    out
  }

  kuu <- block(
    d2_coef = D * B - eps * (a^4) * (B^2) * D / (C^3),
    d1_coef = -eta * D * (1 + D^2) -
      eps * (a^4) * B * D * eta * (3 - 7 * D^2) / (C^3),
    diag_coef = -(eta^2) * D / B -
      shell_nu * (a^2) * D / C -
      F * kappa * (a^2) * (D^3) / C +
      eps * (
        -(a^4) * (eta^2) * D / (A * C^2) +
          F * gamma * kappa * (a^6) * (D^3) * B / (C^4)
      )
  )

  kuw <- block(
    d1_coef = -a * D * sqrt(B) * ((1 + F * kappa) * D^2 + shell_nu) / sqrt(C) -
      eps * (a^5) * A * (B^(3 / 2)) * (1 + F * kappa * gamma) / (C^(9 / 2)),
    diag_coef = a * eta * sqrt(B) * (4 * A + B) / (C^(5 / 2)) +
      eps * (a^5) * A * sqrt(B) * eta *
      (1 / (A^2) - 6 / (C^2) + 9 * A / (C^3)) / (C^(5 / 2))
  )

  kub <- block(
    d2_coef = eps * (a^2) * (B^2) / (C^2),
    d1_coef = -4 * eps * (a^2) * eta * A * B / (C^3),
    diag_coef = F * kappa * D^2 +
      eps * ((a^2) * (eta^2) / (C^2) - F * gamma * kappa * (a^4) * B / (C^4))
  )

  kwu <- block(
    d1_coef = -a * sqrt(B) * (A * (1 + F * kappa) + shell_nu * C) / (C^(3 / 2)) +
      eps * (a^5) * A * (B^(3 / 2)) * (1 + F * kappa * gamma) / (C^(9 / 2)),
    diag_coef = a * eta *
      (1 + (shell_nu + F * kappa) * D^2 - 3 * F * kappa * A * B / (C^2)) /
      sqrt(B * C) +
      eps * eta * (a^5) * sqrt(B) *
      (1 + 3 * F * kappa * gamma * D^4 * (2 - 3 * D^2)) /
      (A * C^(5 / 2))
  )

  kww <- block(
    d2_coef = F * kappa * D * B -
      eps * F * kappa * gamma * (a^4) * D * (B^2) / (C^3),
    d1_coef = -F * kappa * D * (1 + D^2) * eta -
      eps * F * kappa * gamma * (a^4) * D * B * eta * (3 - 7 * D^2) / (C^3),
    diag_coef = -(a^2) * D * (1 + D^4 + 2 * shell_nu * D^2) / A +
      eps * (a^6) * (D^3) * B / C * (1 / (C^3) - 1 / (A^3))
  )

  kwb <- block(
    d1_coef = F * kappa * sqrt(A * B) / a -
      eps * (a^3) * sqrt(A * B) * B * (F * kappa * gamma + 1) / (C^3),
    diag_coef = -F * kappa * eta * sqrt(A / B) / a -
      eps * (a^3) * sqrt(B) * eta *
      (1 + 3 * F * kappa * gamma * D^2 * (1 - 2 * D^2)) /
      (C^2 * sqrt(A))
  )

  kbu <- block(
    d2_coef = eps * (a^2) * (B^2) / (C^2),
    d1_coef = -4 * eps * (a^2) * A * B * eta / (C^3),
    diag_coef = F * kappa * D^2 +
      eps * ((a^2) * (eta^2) / (C^2) - F * kappa * gamma * (a^4) * A * B / (C^4))
  )

  kbw <- block(
    d1_coef = -F * kappa * sqrt(A * B) / a +
      eps * (a^3) * sqrt(B) * (1 + F * kappa * gamma) * A * B / (C^3 * sqrt(A)),
    diag_coef = eps * (a^3) * sqrt(B) * eta * C *
      (-1 + 3 * D^2 * (1 - 2 * D^2)) / (C^3 * sqrt(A))
  )

  kbb <- block(
    d2_coef = eps * D * (B^2),
    d1_coef = -eps * D * eta * (1 + D^2),
    diag_coef = -F * kappa * D * C / (a^2) +
      eps * D * ((eta^2) / B - (a^2) / C + F * kappa * gamma * (a^2) * B / (C^2))
  )

  zmat <- matrix(0, nrow = length(eta), ncol = length(eta))
  muu <- diagm(1 + eps * (a^4) / (C^2))
  mub <- diagm(eps * (1 + D^2))
  mww <- diagm(sqrt(A * C) / (a^2) * (1 + eps * (a^4) / (C^2)))
  mbu <- diagm(eps * (1 + D^2))
  mbb <- diagm(eps * sqrt(A * C) / (a^2))

  structural_matrix <- rbind(
    cbind(kuu, kuw, kub),
    cbind(kwu, kww, kwb),
    cbind(kbu, kbw, kbb)
  )
  mass_matrix <- rbind(
    cbind(muu, zmat, mub),
    cbind(zmat, mww, zmat),
    cbind(mbu, zmat, mbb)
  )
  dynamic_matrix <- structural_matrix + (omega_hat^2) * mass_matrix

  load_scale_q <- -(1 - shell_nu^2) * (geom$semimajor_length^2) *
    sqrt(A * C) / (shell_E * geom$shell_thickness * a^2)
  load_scale_m <- -(1 - shell_nu^2) * sqrt(A * C) /
    (shell_E * geom$shell_thickness * a^2)

  list(
    eta = eta,
    geometry = geom,
    shell = list(
      density = shell_density,
      E = shell_E,
      nu = shell_nu,
      cp = cp,
      F = F,
      gamma = gamma,
      kappa = kappa
    ),
    structural_blocks = list(
      K_uu = kuu,
      K_uw = kuw,
      K_u_beta = kub,
      K_wu = kwu,
      K_ww = kww,
      K_w_beta = kwb,
      K_beta_u = kbu,
      K_beta_w = kbw,
      K_beta_beta = kbb
    ),
    mass_blocks = list(
      M_uu = muu,
      M_uw = zmat,
      M_u_beta = mub,
      M_wu = zmat,
      M_ww = mww,
      M_w_beta = zmat,
      M_beta_u = mbu,
      M_beta_w = zmat,
      M_beta_beta = mbb
    ),
    dynamic_matrix = dynamic_matrix,
    load_scales = list(
      q_u = load_scale_q,
      q_w = load_scale_q,
      m_beta_eta = load_scale_m
    ),
    nondimensional_frequency = omega_hat
  )
}

#' Resolve the sibling validation repository used for the ESPSMS rebuild.
#' @noRd
.espsms_validation_repo_root <- function(validation_repo = NULL) {
  if (!is.null(validation_repo)) {
    root <- normalizePath(validation_repo, winslash = "/", mustWork = FALSE)
    if (!dir.exists(root)) {
      stop(
        "The requested ESPSMS validation repository does not exist: ",
        validation_repo,
        call. = FALSE
      )
    }
    return(root)
  }

  env_root <- Sys.getenv("ACOUSTICTS_VALIDATION_REPO", unset = NA_character_)
  if (is.character(env_root) && length(env_root) == 1L &&
    !is.na(env_root) && nzchar(env_root) && dir.exists(env_root)) {
    return(normalizePath(env_root, winslash = "/", mustWork = TRUE))
  }

  sibling_root <- normalizePath(
    file.path(dirname(normalizePath(getwd(), winslash = "/", mustWork = TRUE)), "acousticTSValidation"),
    winslash = "/",
    mustWork = FALSE
  )
  if (dir.exists(sibling_root)) {
    return(sibling_root)
  }

  stop(
    paste(
      "Could not locate the 'acousticTSValidation' repository needed for the",
      "ESPSMS hybrid reference bridge. Set 'validation_repo' explicitly or",
      "define ACOUSTICTS_VALIDATION_REPO."
    ),
    call. = FALSE
  )
}

#' Convert an ESPSMS body description into the canonical hybrid-spec schema.
#' @noRd
.espsms_validation_spec <- function(body, frequency_hz, case_id = "espsms_package_bridge") {
  list(
    case_id = case_id,
    family = "ESPSMS",
    workflow = "package_bridge_hybrid_reference",
    recommended_backend = "hybrid_shell_fem_bem",
    geometry = list(
      type = "prolate_spheroid",
      semimajor_length_m = as.numeric(body$semimajor_length)[1],
      semiminor_length_m = as.numeric(body$semiminor_length)[1],
      total_length_m = 2 * as.numeric(body$semimajor_length)[1],
      aspect_ratio = as.numeric(body$semimajor_length)[1] /
        as.numeric(body$semiminor_length)[1],
      shell_thickness_m = as.numeric(body$shell_thickness)[1]
    ),
    media = list(
      surrounding = list(
        density_kg_m3 = as.numeric(body$medium_density)[1],
        sound_speed_m_s = as.numeric(body$medium_sound_speed)[1]
      ),
      shell = list(
        density_kg_m3 = as.numeric(body$shell_density)[1],
        sound_speed_m_s = as.numeric(body$shell_sound_speed %||% sqrt(
          as.numeric(body$shell_E)[1] /
            (as.numeric(body$shell_density)[1] * (1 - as.numeric(body$shell_nu)[1]^2))
        ))[1],
        youngs_modulus_Pa = as.numeric(body$shell_E)[1],
        poissons_ratio = as.numeric(body$shell_nu)[1]
      ),
      interior = list(
        density_kg_m3 = as.numeric(body$fluid_density)[1],
        sound_speed_m_s = as.numeric(body$fluid_sound_speed)[1]
      )
    ),
    frequency_hz = as.list(as.numeric(frequency_hz)),
    notes = list(
      "Temporary internal package bridge to the external ESPSMS hybrid reference.",
      "This is for rebuild comparison only and is not the public ESPSMS direct solver."
    )
  )
}

#' Run the external hybrid ESPSMS reference for a monostatic ladder.
#' @noRd
.espsms_hybrid_reference_monostatic_from_body <- function(body,
                                                          frequency,
                                                          theta_body_deg = 0,
                                                          n_eta = 129L,
                                                          mesh_n_u = 12L,
                                                          mesh_n_v = 16L,
                                                          validation_repo = NULL,
                                                          python = Sys.getenv("ACOUSTICTS_PYTHON", unset = "python")) {
  repo_root <- .espsms_validation_repo_root(validation_repo)
  helper_py <- file.path(
    repo_root,
    "tools",
    "implementation-figures",
    "helpers",
    "prolate_elastic_shell_hybrid_solver.py"
  )
  if (!file.exists(helper_py)) {
    stop("Missing ESPSMS hybrid helper: ", helper_py, call. = FALSE)
  }

  tmp_dir <- tempfile("espsms-hybrid-bridge-")
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(tmp_dir, recursive = TRUE, force = TRUE), add = TRUE)

  spec_path <- file.path(tmp_dir, "spec.json")
  out_path <- file.path(tmp_dir, "monostatic.csv")
  jsonlite::write_json(
    .espsms_validation_spec(body, frequency_hz = frequency),
    path = spec_path,
    auto_unbox = TRUE,
    pretty = TRUE
  )

  out <- system2(
    python,
    args = c(
      helper_py,
      "--spec", spec_path,
      "--case-id", "espsms_package_bridge",
      "--theta-body-deg", as.character(theta_body_deg),
      "--n-eta", as.character(as.integer(n_eta)),
      "--mesh-n-u", as.character(as.integer(mesh_n_u)),
      "--mesh-n-v", as.character(as.integer(mesh_n_v)),
      "--mode", "monostatic",
      "--monostatic-out", out_path
    ),
    stdout = TRUE,
    stderr = TRUE
  )
  status <- attr(out, "status")
  if (!is.null(status) && status != 0) {
    stop(
      paste(c("ESPSMS hybrid reference bridge failed:", out), collapse = "\n"),
      call. = FALSE
    )
  }
  if (!file.exists(out_path)) {
    stop(
      "ESPSMS hybrid reference bridge did not produce the monostatic output file.",
      call. = FALSE
    )
  }

  utils::read.csv(out_path, stringsAsFactors = FALSE)
}

#' Run the external hybrid ESPSMS reference for a monostatic ladder.
#' @noRd
.espsms_hybrid_reference_monostatic <- function(object,
                                                frequency,
                                                sound_speed_sw = .SEAWATER_SOUND_SPEED_DEFAULT,
                                                density_sw = .SEAWATER_DENSITY_DEFAULT,
                                                theta_body_deg = 0,
                                                n_eta = 129L,
                                                mesh_n_u = 12L,
                                                mesh_n_v = 16L,
                                                validation_repo = NULL,
                                                python = Sys.getenv("ACOUSTICTS_PYTHON", unset = "python")) {
  body <- .elastic_shell_prolate_body(
    object = object,
    sound_speed_sw = sound_speed_sw,
    density_sw = density_sw
  )

  .espsms_hybrid_reference_monostatic_from_body(
    body = body,
    frequency = frequency,
    theta_body_deg = theta_body_deg,
    n_eta = n_eta,
    mesh_n_u = mesh_n_u,
    mesh_n_v = mesh_n_v,
    validation_repo = validation_repo,
    python = python
  )
}

#' Run the external hybrid ESPSMS reference grid for one stored frequency.
#' @noRd
.espsms_hybrid_reference_grid_from_body <- function(body,
                                                    frequency,
                                                    theta_body_deg = 0,
                                                    n_eta = 129L,
                                                    mesh_n_u = 12L,
                                                    mesh_n_v = 16L,
                                                    n_theta = 91L,
                                                    n_phi = 181L,
                                                    validation_repo = NULL,
                                                    python = Sys.getenv("ACOUSTICTS_PYTHON", unset = "python")) {
  repo_root <- .espsms_validation_repo_root(validation_repo)
  helper_py <- file.path(
    repo_root,
    "tools",
    "implementation-figures",
    "helpers",
    "prolate_elastic_shell_hybrid_solver.py"
  )
  if (!file.exists(helper_py)) {
    stop("Missing ESPSMS hybrid helper: ", helper_py, call. = FALSE)
  }

  tmp_dir <- tempfile("espsms-hybrid-grid-")
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(tmp_dir, recursive = TRUE, force = TRUE), add = TRUE)

  spec_path <- file.path(tmp_dir, "spec.json")
  out_path <- file.path(tmp_dir, "grid.csv")
  jsonlite::write_json(
    .espsms_validation_spec(body, frequency_hz = frequency, case_id = "espsms_package_grid_bridge"),
    path = spec_path,
    auto_unbox = TRUE,
    pretty = TRUE
  )

  out <- system2(
    python,
    args = c(
      helper_py,
      "--spec", spec_path,
      "--case-id", "espsms_package_grid_bridge",
      "--theta-body-deg", as.character(theta_body_deg),
      "--n-eta", as.character(as.integer(n_eta)),
      "--mesh-n-u", as.character(as.integer(mesh_n_u)),
      "--mesh-n-v", as.character(as.integer(mesh_n_v)),
      "--mode", "bistatic_grid",
      "--grid-frequency-hz", as.character(as.numeric(frequency)[1]),
      "--n-theta", as.character(as.integer(n_theta)),
      "--n-phi", as.character(as.integer(n_phi)),
      "--bistatic-out", out_path
    ),
    stdout = TRUE,
    stderr = TRUE
  )
  status <- attr(out, "status")
  if (!is.null(status) && status != 0) {
    stop(
      paste(c("ESPSMS hybrid grid bridge failed:", out), collapse = "\n"),
      call. = FALSE
    )
  }
  if (!file.exists(out_path)) {
    stop(
      "ESPSMS hybrid grid bridge did not produce the bistatic output file.",
      call. = FALSE
    )
  }

  grid_df <- utils::read.csv(out_path, stringsAsFactors = FALSE)
  theta_vals <- sort(unique(as.numeric(grid_df$theta_scatter_rad)))
  phi_vals <- sort(unique(as.numeric(grid_df$phi_scatter_rad)))
  f_vals <- complex(
    real = as.numeric(grid_df$f_scat_real),
    imaginary = as.numeric(grid_df$f_scat_imag)
  )
  sigma_vals <- as.numeric(grid_df$sigma_scat)

  list(
    family = "espsms_hybrid_grid",
    theta_body = pi / 2,
    phi_body = if ((theta_body_deg %% 360) < 90 || (theta_body_deg %% 360) > 270) 0 else pi,
    frequency = as.numeric(frequency)[1],
    theta_scatter = theta_vals,
    phi_scatter = phi_vals,
    f_scat = matrix(
      f_vals,
      nrow = length(theta_vals),
      ncol = length(phi_vals),
      byrow = TRUE
    ),
    sigma_scat = matrix(
      sigma_vals,
      nrow = length(theta_vals),
      ncol = length(phi_vals),
      byrow = TRUE
    )
  )
}

#' Resolve the direct-solver profile for the current ESPSMS branch.
#' @noRd
.espsms_resolve_solver_profile <- function(frequency,
                                           sound_speed_sw,
                                           shape_parameters,
                                           theta_body,
                                           n_max = NULL) {
  if (!is.null(n_max)) {
    if (!is.numeric(n_max) || any(!is.finite(n_max)) || any(n_max < 0)) {
      stop(
        "'n_max' must be a non-negative numeric scalar or vector.",
        call. = FALSE
      )
    }
    if (length(n_max) == 1L) {
      resolved_n_max <- rep(as.integer(round(n_max)), length(frequency))
    } else if (length(n_max) != length(frequency)) {
      stop(
        "'n_max' must have length 1 or the same length as 'frequency'.",
        call. = FALSE
      )
    } else {
      resolved_n_max <- as.integer(round(n_max))
    }
  } else {
    shape_source <- shape_parameters$shell %||% shape_parameters
    major_radius <- as.numeric(
      shape_source$semimajor_length %||% (shape_source$length / 2)
    )[1]
    if (!is.finite(major_radius) || major_radius <= 0) {
      stop(
        "ESPSMS requires a finite positive prolate semimajor length.",
        call. = FALSE
      )
    }

    ka_major <- wavenumber(frequency, sound_speed_sw) * major_radius
    cos_theta <- abs(cos(theta_body))
    sin_theta <- abs(sin(theta_body))

    resolved_n_max <- ifelse(
      ka_major <= 5,
      ifelse(
        cos_theta >= 0.75,
        48L,
        ifelse(cos_theta >= 0.25, 32L, 12L)
      ),
      ifelse(sin_theta >= 0.95, 48L, 32L)
    )
  }

  resolved_n_max <- as.integer(pmin(64L, pmax(12L, resolved_n_max)))
  shape_source <- shape_parameters$shell %||% shape_parameters
  major_radius <- as.numeric(
    shape_source$semimajor_length %||% (shape_source$length / 2)
  )[1]
  ka_major <- wavenumber(frequency, sound_speed_sw) * major_radius
  cos_theta <- abs(cos(theta_body))
  sin_theta <- abs(sin(theta_body))
  high_ka <- ka_major > 5

  collocation_multiplier <- rep(as.integer(8L), length(frequency))
  svd_rel_tol <- rep(1e-10, length(frequency))
  enhanced_projection <- high_ka & sin_theta >= 0.45
  collocation_multiplier[enhanced_projection] <- 12L
  svd_rel_tol[enhanced_projection] <- 1e-5

  list(
    n_max = resolved_n_max,
    collocation_multiplier = collocation_multiplier,
    svd_rel_tol = svd_rel_tol
  )
}

#' Initialize a prolate elastic-shell object for the ESPSMS branch.
#' @param object ESS-class object.
#' @param frequency Frequency vector (Hz).
#' @param sound_speed_sw Surrounding-medium sound speed (m/s).
#' @param density_sw Surrounding-medium density (kg/m^3).
#' @param n_max Optional retained modal cutoff.
#' @param solver_backend ESPSMS backend selector. Use `"disabled"` to keep the
#'   branch off, or `"hybrid_reference"` to use the external coupled shell-fluid
#'   hybrid solver during the rebuild.
#' @param validation_repo Optional path to the `acousticTSValidation`
#'   repository when using `solver_backend = "hybrid_reference"`.
#' @param python Python executable to use for the hybrid reference backend.
#' @param theta_body_deg Optional incidence angle override for the current
#'   experimental hybrid backend. At present only `0` or `180` are supported.
#' @param n_eta Shell meridional grid size for the hybrid reference backend.
#' @param mesh_n_u Structured surface-mesh meridional count for the hybrid
#'   reference backend.
#' @param mesh_n_v Structured surface-mesh azimuthal count for the hybrid
#'   reference backend.
#' @noRd
espsms_initialize <- function(object,
                              frequency,
                              sound_speed_sw = .SEAWATER_SOUND_SPEED_DEFAULT,
                              density_sw = .SEAWATER_DENSITY_DEFAULT,
                              n_max = NULL,
                              solver_backend = c("disabled", "hybrid_reference"),
                              validation_repo = NULL,
                              python = Sys.getenv("ACOUSTICTS_PYTHON", unset = "python"),
                              theta_body_deg = NULL,
                              n_eta = 129L,
                              mesh_n_u = 12L,
                              mesh_n_v = 16L) {
  solver_backend <- match.arg(solver_backend)
  if (identical(solver_backend, "disabled")) {
    stop(
      paste(
        "ESPSMS is temporarily disabled.",
        "The previous projected thin-shell approximation did not validate as a",
        "genuine elastic-shelled prolate spheroid model and has been retired",
        "pending a shell-only modal/numerical rebuild."
      ),
      call. = FALSE
    )
  }

  body <- .elastic_shell_prolate_body(
    object = object,
    sound_speed_sw = sound_speed_sw,
    density_sw = density_sw
  )
  shape_parameters <- acousticTS::extract(object, "shape_parameters")
  acoustics <- .init_acoustics_df(
    frequency,
    k_sw = sound_speed_sw,
    k_body = body$fluid_sound_speed
  )
  solver_profile <- .espsms_resolve_solver_profile(
    frequency = frequency,
    sound_speed_sw = sound_speed_sw,
    shape_parameters = shape_parameters,
    theta_body = body$theta,
    n_max = n_max
  )
  acoustics$n_max <- solver_profile$n_max

  .init_model_slots(
    object = object,
    model_name = "ESPSMS",
    frequency = frequency,
    model_parameters = list(
      parameters = list(
        acoustics = acoustics,
        shell_theory = "thin_isotropic_love_kirchhoff",
        solver_profile = solver_profile,
        solver_backend = solver_backend,
        validation_repo = validation_repo,
        python = python,
        theta_body_deg = theta_body_deg,
        n_eta = as.integer(n_eta),
        mesh_n_u = as.integer(mesh_n_u),
        mesh_n_v = as.integer(mesh_n_v)
      ),
      medium = .init_medium_params(sound_speed_sw, density_sw),
      body = body
    ),
    result_cols = c(
      "f_bs",
      "sigma_bs",
      "TS",
      "n_max",
      "solver_backend",
      "n_surface_dof",
      "n_shell_dof",
      "shell_w_abs_max",
      "pressure_jump_abs_max"
    )
  )
}

#' Evaluate the thin-shell prolate ESPSMS branch.
#' @param object ESS-class object initialized for ESPSMS.
#' @noRd
ESPSMS <- function(object) {
  model <- acousticTS::extract(object, "model_parameters")$ESPSMS
  acoustics <- model$parameters$acoustics
  solver_profile <- model$parameters$solver_profile
  body <- model$body
  solver_backend <- model$parameters$solver_backend %||% "disabled"

  if (identical(solver_backend, "disabled")) {
    stop(
      paste(
        "ESPSMS is temporarily disabled.",
        "The previous projected thin-shell approximation did not validate as a",
        "genuine elastic-shelled prolate spheroid model and has been retired",
        "pending a shell-only modal/numerical rebuild."
      ),
      call. = FALSE
    )
  }

  if (!identical(solver_backend, "hybrid_reference")) {
    stop(
      "Unsupported ESPSMS solver backend: ", solver_backend,
      call. = FALSE
    )
  }

  theta_body_deg <- model$parameters$theta_body_deg %||% NA_real_
  if (is.na(theta_body_deg)) {
    theta_body_deg <- if (isTRUE(all.equal(abs(cos(body$theta)), 1, tolerance = 1e-8))) {
      if (cos(body$theta) >= 0) 0 else 180
    } else {
      stop(
        paste(
          "The current ESPSMS hybrid-reference backend only supports axial",
          "incidence (theta_body = 0 or pi). Supply 'theta_body_deg = 0' or",
          "'theta_body_deg = 180' explicitly when initializing the model."
        ),
        call. = FALSE
      )
    }
  }
  if (!isTRUE(all.equal(theta_body_deg %% 180, 0, tolerance = 1e-8))) {
    stop(
      paste(
        "The current ESPSMS hybrid-reference backend only supports axial",
        "incidence (theta_body = 0 or pi)."
      ),
      call. = FALSE
    )
  }
  theta_body_deg <- ifelse((theta_body_deg %% 360) >= 180, 180, 0)

  hybrid <- .espsms_hybrid_reference_monostatic(
    object = object,
    frequency = acoustics$frequency,
    sound_speed_sw = model$medium$sound_speed,
    density_sw = model$medium$density,
    theta_body_deg = theta_body_deg,
    n_eta = model$parameters$n_eta %||% 129L,
    mesh_n_u = model$parameters$mesh_n_u %||% 12L,
    mesh_n_v = model$parameters$mesh_n_v %||% 16L,
    validation_repo = model$parameters$validation_repo,
    python = model$parameters$python %||% Sys.getenv("ACOUSTICTS_PYTHON", unset = "python")
  )

  methods::slot(object, "model")$ESPSMS <- data.frame(
    frequency = acoustics$frequency,
    f_bs = complex(real = hybrid$f_bs_real, imaginary = hybrid$f_bs_imag),
    sigma_bs = hybrid$sigma_bs,
    TS = hybrid$TS,
    n_max = acoustics$n_max,
    solver_backend = solver_backend,
    n_surface_dof = hybrid$n_surface_dof,
    n_shell_dof = hybrid$n_shell_dof,
    shell_w_abs_max = hybrid$shell_w_abs_max,
    pressure_jump_abs_max = hybrid$pressure_jump_abs_max
  )

  object
}

#' Backward-compatible initializer alias for older experimental naming.
#' @noRd
epsms_initialize <- function(...) {
  espsms_initialize(...)
}

#' Backward-compatible solver alias for older experimental naming.
#' @noRd
EPSMS <- function(...) {
  ESPSMS(...)
}
