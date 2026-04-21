#' Elastic-shell prolate spheroidal model (ESPSMS)
#'
#' @description
#' `ESPSMS` is currently disabled in `acousticTS` while the elastic-shelled
#' prolate implementation is rebuilt and validated in the standalone
#' `acousticValidator` repository.
#'
#' @section Usage:
#' Calls routed through `model = "espsms"` currently fail fast.
#'
#' @section Implementation route:
#' The package retains only the shell-theory helper routines needed for the
#' ongoing rebuild. Validation solvers, external bridges, and experimental
#' hybrid workflows live outside `acousticTS`.
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
.espsms_uniform_eta_grid <- function(eta_grid_size = 129L,
                                     pole_offset = 1e-4) {
  if (!is.numeric(eta_grid_size) || length(eta_grid_size) != 1L ||
    !is.finite(eta_grid_size) || eta_grid_size < 5) {
    stop("'eta_grid_size' must be one integer >= 5.", call. = FALSE)
  }
  if (!is.numeric(pole_offset) || length(pole_offset) != 1L ||
    !is.finite(pole_offset) || pole_offset <= 0 || pole_offset >= 0.5) {
    stop("'pole_offset' must lie in (0, 0.5).", call. = FALSE)
  }

  seq(
    from = -1 + pole_offset,
    to = 1 - pole_offset,
    length.out = as.integer(eta_grid_size)
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
                                              eta_grid_size = 129L,
                                              pole_offset = 1e-4) {
  eta <- .espsms_uniform_eta_grid(
    eta_grid_size = eta_grid_size,
    pole_offset = pole_offset
  )
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

#' Initialize a prolate elastic-shell object for the ESPSMS branch.
#' @param object ESS-class object.
#' @param frequency Frequency vector (Hz).
#' @param sound_speed_sw Surrounding-medium sound speed (m/s).
#' @param density_sw Surrounding-medium density (kg/m^3).
#' @noRd
espsms_initialize <- function(object,
                              frequency,
                              sound_speed_sw = .SEAWATER_SOUND_SPEED_DEFAULT,
                              density_sw = .SEAWATER_DENSITY_DEFAULT) {
  stop(
    paste(
      "ESPSMS is temporarily disabled.",
      "The elastic-shelled prolate rebuild and all validator workflows now live",
      "in the standalone acousticValidator repository until the package-side",
      "model is restored as a validated direct implementation."
    ),
    call. = FALSE
  )
}

#' Evaluate the thin-shell prolate ESPSMS branch.
#' @param object ESS-class object initialized for ESPSMS.
#' @noRd
ESPSMS <- function(object) {
  stop(
    paste(
      "ESPSMS is temporarily disabled.",
      "The elastic-shelled prolate rebuild and all validator workflows now live",
      "in the standalone acousticValidator repository until the package-side",
      "model is restored as a validated direct implementation."
    ),
    call. = FALSE
  )
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
