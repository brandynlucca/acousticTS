################################################################################
# Viscous-elastic spherical scattering model
################################################################################
#' Viscous-elastic spherical model (VESMS)
#'
#' @description
#' Computes backscatter from a gas-filled spherical inclusion surrounded by an
#' elastic shell and an outer viscous layer, following the wideband
#' backscattering model used by Khodabandeloo et al. (2021). The model is
#' parameterized on `ESS` objects, where the `fluid` slot represents the inner
#' gas sphere and the `shell` slot represents the elastic shell. The external
#' viscous layer is supplied as model-specific arguments.
#'
#' @section Usage:
#' This model is accessed via:
#' \preformatted{
#' target_strength(
#'   ...,
#'   model = "VESMS",
#'   sound_speed_sw,
#'   density_sw,
#'   sound_speed_viscous,
#'   density_viscous,
#'   shear_viscosity_viscous,
#'   bulk_viscosity_viscous,
#'   radius_viscous,
#'   viscous_thickness,
#'   m_limit
#' )
#' }
#'
#' @section Arguments:
#' \describe{
#'   \item{\code{sound_speed_sw}}{Seawater sound speed (\eqn{m~s^{-1}}).}
#'   \item{\code{density_sw}}{Seawater density (\eqn{kg~m^{-3}}).}
#'   \item{\code{sound_speed_viscous}}{Compressional sound speed in the outer
#'   viscous layer (\eqn{m~s^{-1}}).}
#'   \item{\code{density_viscous}}{Density of the outer viscous layer
#'   (\eqn{kg~m^{-3}}).}
#'   \item{\code{shear_viscosity_viscous}}{Shear viscosity of the outer viscous
#'   layer (\eqn{kg~m^{-1}~s^{-1}}).}
#'   \item{\code{bulk_viscosity_viscous}}{Optional bulk viscosity of the outer
#'   viscous layer. When omitted, it defaults to the supplied shear viscosity to
#'   match the reference implementation.}
#'   \item{\code{radius_viscous}}{Optional outer radius of the viscous layer
#'   (\eqn{R_2}, m).}
#'   \item{\code{viscous_thickness}}{Optional thickness of the viscous layer
#'   relative to the shell outer radius. Supply either this or
#'   \code{radius_viscous}, not both.}
#'   \item{\code{m_limit}}{Optional truncation limit for the modal summation.
#'   A value of 2 retains the \eqn{m = 0, 1, 2} terms. When omitted, the model
#'   uses \eqn{\max(2, \mathrm{round}(k_1 R_2) + 10)} at each frequency.}
#' }
#'
#' @details
#' When \code{radius_viscous} and \code{viscous_thickness} are both omitted, the
#' model estimates the viscous-layer radius from the neutral-buoyancy relation
#' used by Khodabandeloo et al. (2021, Eq. 18).
#'
#' @references
#' Khodabandeloo, B., Agersted, M.D., Klevjer, T., Macaulay, G.J., and Melle,
#' W. (2021). Estimating target strength and physical characteristics of
#' gas-bearing mesopelagic fish from wideband in situ echoes using a
#' viscous-elastic scattering model. \emph{The Journal of the Acoustical
#' Society of America}, 149, 673-691.
#'
#' Feuillade, C., and Nero, R.W. (1998). A viscous-elastic swimbladder model for
#' describing enhanced-frequency resonance scattering from fish. \emph{The
#' Journal of the Acoustical Society of America}, 103, 3245-3255.
#'
#' Anson, D.S., and Chivers, R.C. (1993). An irregular frequencies solution to
#' acoustic scattering by elastic spheres. \emph{The Journal of the Acoustical
#' Society of America}, 93, 118-123.
#'
#' @name VESMS
#' @aliases vesms VESMS
#' @docType data
#' @keywords models acoustics
NULL

#' @noRd
.vesms_is_spherical <- function(object) {
  shape <- extract(object, "shape_parameters")
  shape_tag <- shape$shape %||% ""
  any(grepl("sphere", tolower(as.character(shape_tag))))
}

#' @noRd
.vesms_validate_outer_radius <- function(radius_viscous, radius_shell) {
  if (!is.numeric(radius_viscous) ||
      length(radius_viscous) != 1L ||
      !is.finite(radius_viscous) ||
      radius_viscous <= radius_shell) {
    stop(
      "VESMS requires the viscous-layer radius to be a finite scalar larger ",
      "than the shell radius."
    )
  }

  radius_viscous
}

#' @noRd
.vesms_neutral_radius <- function(radius_gas,
                                  density_sw,
                                  density_viscous,
                                  density_gas) {
  denom <- density_viscous - density_sw
  if (!is.finite(denom) || abs(denom) <= sqrt(.Machine$double.eps)) {
    stop(
      "Neutral-buoyancy estimation of the viscous-layer radius is undefined ",
      "when 'density_viscous' equals the surrounding-medium density. Supply ",
      "'radius_viscous' or 'viscous_thickness' explicitly."
    )
  }

  radius_gas * (1 + (density_sw - density_gas) / denom)^(1 / 3)
}

#' @noRd
.vesms_sph_j <- function(order, z) {
  z <- as.complex(z)
  eps <- sqrt(.Machine$double.eps)
  out <- complex(length = length(z))

  for (i in seq_along(z)) {
    zi <- z[i]
    if (Mod(zi) <= eps) {
      out[i] <- if (order == 0L) 1 + 0i else 0 + 0i
      next
    }

    j0 <- sin(zi) / zi
    if (order == 0L) {
      out[i] <- j0
      next
    }

    j1 <- sin(zi) / zi^2 - cos(zi) / zi
    if (order == 1L) {
      out[i] <- j1
      next
    }

    jm1 <- j0
    jn <- j1
    for (n in seq_len(order - 1L)) {
      jp1 <- ((2 * n + 1) / zi) * jn - jm1
      jm1 <- jn
      jn <- jp1
    }
    out[i] <- jn
  }

  out
}

#' @noRd
.vesms_sph_y <- function(order, z) {
  z <- as.complex(z)
  eps <- sqrt(.Machine$double.eps)
  out <- complex(length = length(z))

  for (i in seq_along(z)) {
    zi <- z[i]
    if (Mod(zi) <= eps) {
      out[i] <- complex(real = NaN, imaginary = NaN)
      next
    }

    y0 <- -cos(zi) / zi
    if (order == 0L) {
      out[i] <- y0
      next
    }

    y1 <- -cos(zi) / zi^2 - sin(zi) / zi
    if (order == 1L) {
      out[i] <- y1
      next
    }

    ym1 <- y0
    yn <- y1
    for (n in seq_len(order - 1L)) {
      yp1 <- ((2 * n + 1) / zi) * yn - ym1
      ym1 <- yn
      yn <- yp1
    }
    out[i] <- yn
  }

  out
}

#' @noRd
.vesms_sph_h <- function(order, z) {
  .vesms_sph_j(order, z) + 1i * .vesms_sph_y(order, z)
}

#' @noRd
.vesms_sph_jd <- function(order, z) {
  z <- as.complex(z)
  if (order == 0L) {
    return(-.vesms_sph_j(1L, z))
  }
  .vesms_sph_j(order - 1L, z) - ((order + 1) / z) * .vesms_sph_j(order, z)
}

#' @noRd
.vesms_sph_yd <- function(order, z) {
  z <- as.complex(z)
  if (order == 0L) {
    return(-.vesms_sph_y(1L, z))
  }
  .vesms_sph_y(order - 1L, z) - ((order + 1) / z) * .vesms_sph_y(order, z)
}

#' @noRd
.vesms_sph_hd <- function(order, z) {
  .vesms_sph_jd(order, z) + 1i * .vesms_sph_yd(order, z)
}

#' @noRd
.vesms_sph_jdd <- function(order, z) {
  z <- as.complex(z)
  jn <- .vesms_sph_j(order, z)
  jnd <- .vesms_sph_jd(order, z)
  ((order * (order + 1) / z^2) - 1) * jn - (2 / z) * jnd
}

#' @noRd
.vesms_sph_hdd <- function(order, z) {
  z <- as.complex(z)
  hn <- .vesms_sph_h(order, z)
  hnd <- .vesms_sph_hd(order, z)
  ((order * (order + 1) / z^2) - 1) * hn - (2 / z) * hnd
}

#' @noRd
.vesms_solve_system <- function(M, F) {
  direct <- tryCatch(
    solve(M, F, tol = 0),
    error = function(e) NULL
  )
  if (!is.null(direct)) {
    return(direct)
  }

  qr_fit <- tryCatch(
    qr.solve(M, F),
    error = function(e) NULL
  )
  if (!is.null(qr_fit)) {
    return(qr_fit)
  }

  svd_fit <- tryCatch(
    {
      dec <- svd(M)
      tol <- max(dim(M)) * .Machine$double.eps * max(dec$d)
      keep <- dec$d > tol
      if (!any(keep)) {
        stop("All singular values are below the numerical tolerance.")
      }
      inv_d <- rep(0, length(dec$d))
      inv_d[keep] <- 1 / dec$d[keep]
      dec$v %*% (diag(inv_d, nrow = length(inv_d)) %*% (Conj(t(dec$u)) %*% F))
    },
    error = function(e) NULL
  )

  if (!is.null(svd_fit)) {
    return(svd_fit)
  }

  stop(
    "VESMS was unable to solve the modal boundary system. This usually ",
    "indicates a numerically singular mode-frequency combination."
  )
}

#' @noRd
.vesms_build_mode_system <- function(order,
                                     omega,
                                     density_sw,
                                     density_viscous,
                                     density_shell,
                                     density_gas,
                                     sound_speed_sw,
                                     sound_speed_viscous,
                                     sound_speed_gas,
                                     radius_viscous,
                                     radius_shell,
                                     radius_gas,
                                     shear_viscosity_viscous,
                                     bulk_viscosity_viscous,
                                     shear_modulus_shell,
                                     lambda_shell) {
  k1 <- omega / sound_speed_sw
  phi2 <- bulk_viscosity_viscous + 4 * shear_viscosity_viscous / 3
  kc2 <- (omega / sound_speed_viscous) *
    (1 - 1i * omega * phi2 /
      (density_viscous * sound_speed_viscous^2))^(-1 / 2)
  ks2 <- (1 + 1i) * sqrt(omega * density_viscous /
    (2 * shear_viscosity_viscous))
  kc3 <- omega * sqrt(density_shell / (lambda_shell + 2 * shear_modulus_shell))
  ks3 <- omega * sqrt(density_shell / shear_modulus_shell)
  k4 <- omega / sound_speed_gas

  j_k1_R2 <- .vesms_sph_j(order, k1 * radius_viscous)
  jd_k1_R2 <- .vesms_sph_jd(order, k1 * radius_viscous)
  h_k1_R2 <- .vesms_sph_h(order, k1 * radius_viscous)
  hd_k1_R2 <- .vesms_sph_hd(order, k1 * radius_viscous)

  j_kc2_R2 <- .vesms_sph_j(order, kc2 * radius_viscous)
  jd_kc2_R2 <- .vesms_sph_jd(order, kc2 * radius_viscous)
  jdd_kc2_R2 <- .vesms_sph_jdd(order, kc2 * radius_viscous)
  h_kc2_R2 <- .vesms_sph_h(order, kc2 * radius_viscous)
  hd_kc2_R2 <- .vesms_sph_hd(order, kc2 * radius_viscous)
  hdd_kc2_R2 <- .vesms_sph_hdd(order, kc2 * radius_viscous)

  j_ks2_R2 <- .vesms_sph_j(order, ks2 * radius_viscous)
  jd_ks2_R2 <- .vesms_sph_jd(order, ks2 * radius_viscous)
  jdd_ks2_R2 <- .vesms_sph_jdd(order, ks2 * radius_viscous)
  h_ks2_R2 <- .vesms_sph_h(order, ks2 * radius_viscous)
  hd_ks2_R2 <- .vesms_sph_hd(order, ks2 * radius_viscous)
  hdd_ks2_R2 <- .vesms_sph_hdd(order, ks2 * radius_viscous)

  j_kc2_R3 <- .vesms_sph_j(order, kc2 * radius_shell)
  jd_kc2_R3 <- .vesms_sph_jd(order, kc2 * radius_shell)
  jdd_kc2_R3 <- .vesms_sph_jdd(order, kc2 * radius_shell)
  h_kc2_R3 <- .vesms_sph_h(order, kc2 * radius_shell)
  hd_kc2_R3 <- .vesms_sph_hd(order, kc2 * radius_shell)
  hdd_kc2_R3 <- .vesms_sph_hdd(order, kc2 * radius_shell)

  j_ks2_R3 <- .vesms_sph_j(order, ks2 * radius_shell)
  jd_ks2_R3 <- .vesms_sph_jd(order, ks2 * radius_shell)
  jdd_ks2_R3 <- .vesms_sph_jdd(order, ks2 * radius_shell)
  h_ks2_R3 <- .vesms_sph_h(order, ks2 * radius_shell)
  hd_ks2_R3 <- .vesms_sph_hd(order, ks2 * radius_shell)
  hdd_ks2_R3 <- .vesms_sph_hdd(order, ks2 * radius_shell)

  j_kc3_R3 <- .vesms_sph_j(order, kc3 * radius_shell)
  jd_kc3_R3 <- .vesms_sph_jd(order, kc3 * radius_shell)
  jdd_kc3_R3 <- .vesms_sph_jdd(order, kc3 * radius_shell)
  h_kc3_R3 <- .vesms_sph_h(order, kc3 * radius_shell)
  hd_kc3_R3 <- .vesms_sph_hd(order, kc3 * radius_shell)
  hdd_kc3_R3 <- .vesms_sph_hdd(order, kc3 * radius_shell)

  j_ks3_R3 <- .vesms_sph_j(order, ks3 * radius_shell)
  jd_ks3_R3 <- .vesms_sph_jd(order, ks3 * radius_shell)
  jdd_ks3_R3 <- .vesms_sph_jdd(order, ks3 * radius_shell)
  h_ks3_R3 <- .vesms_sph_h(order, ks3 * radius_shell)
  hd_ks3_R3 <- .vesms_sph_hd(order, ks3 * radius_shell)
  hdd_ks3_R3 <- .vesms_sph_hdd(order, ks3 * radius_shell)

  j_kc3_R4 <- .vesms_sph_j(order, kc3 * radius_gas)
  jd_kc3_R4 <- .vesms_sph_jd(order, kc3 * radius_gas)
  jdd_kc3_R4 <- .vesms_sph_jdd(order, kc3 * radius_gas)
  h_kc3_R4 <- .vesms_sph_h(order, kc3 * radius_gas)
  hd_kc3_R4 <- .vesms_sph_hd(order, kc3 * radius_gas)
  hdd_kc3_R4 <- .vesms_sph_hdd(order, kc3 * radius_gas)

  j_ks3_R4 <- .vesms_sph_j(order, ks3 * radius_gas)
  jd_ks3_R4 <- .vesms_sph_jd(order, ks3 * radius_gas)
  jdd_ks3_R4 <- .vesms_sph_jdd(order, ks3 * radius_gas)
  h_ks3_R4 <- .vesms_sph_h(order, ks3 * radius_gas)
  hd_ks3_R4 <- .vesms_sph_hd(order, ks3 * radius_gas)
  hdd_ks3_R4 <- .vesms_sph_hdd(order, ks3 * radius_gas)

  j_k4_R4 <- .vesms_sph_j(order, k4 * radius_gas)
  jd_k4_R4 <- .vesms_sph_jd(order, k4 * radius_gas)

  if (order == 0L) {
    M <- matrix(0 + 0i, nrow = 6, ncol = 6)
    F <- matrix(0 + 0i, nrow = 6, ncol = 1)

    M[1, 1] <- k1 * hd_k1_R2
    M[1, 2] <- -kc2 * jd_kc2_R2
    M[1, 3] <- -kc2 * hd_kc2_R2
    F[1, 1] <- -k1 * ((1i)^order) * (2 * order + 1) * jd_k1_R2

    M[2, 1] <- -1i * omega * density_sw * h_k1_R2
    M[2, 2] <- -shear_viscosity_viscous * (2 * kc2^2 - ks2^2) * j_kc2_R2 -
      2 * shear_viscosity_viscous * kc2^2 * jdd_kc2_R2
    M[2, 3] <- -shear_viscosity_viscous * (2 * kc2^2 - ks2^2) * h_kc2_R2 -
      2 * shear_viscosity_viscous * kc2^2 * hdd_kc2_R2
    F[2, 1] <- 1i * omega * density_sw * ((1i)^order) * (2 * order + 1) *
      j_k1_R2

    M[3, 2] <- kc2 * jd_kc2_R3
    M[3, 3] <- kc2 * hd_kc2_R3
    M[3, 4] <- -kc3 * jd_kc3_R3
    M[3, 5] <- -kc3 * hd_kc3_R3

    M[4, 2] <- -1i * omega * shear_viscosity_viscous * (2 * kc2^2 - ks2^2) *
      j_kc2_R3 - 2 * 1i * omega * shear_viscosity_viscous * kc2^2 * jdd_kc2_R3
    M[4, 3] <- -1i * omega * shear_viscosity_viscous * (2 * kc2^2 - ks2^2) *
      h_kc2_R3 - 2 * 1i * omega * shear_viscosity_viscous * kc2^2 * hdd_kc2_R3
    M[4, 4] <- -shear_modulus_shell * (2 * kc3^2 - ks3^2) * j_kc3_R3 -
      2 * shear_modulus_shell * kc3^2 * jdd_kc3_R3
    M[4, 5] <- -shear_modulus_shell * (2 * kc3^2 - ks3^2) * h_kc3_R3 -
      2 * shear_modulus_shell * kc3^2 * hdd_kc3_R3

    M[5, 4] <- kc3 * jd_kc3_R4
    M[5, 5] <- kc3 * hd_kc3_R4
    M[5, 6] <- -k4 * jd_k4_R4

    M[6, 4] <- shear_modulus_shell * (2 * kc3^2 - ks3^2) * j_kc3_R4 +
      2 * shear_modulus_shell * kc3^2 * jdd_kc3_R4
    M[6, 5] <- shear_modulus_shell * (2 * kc3^2 - ks3^2) * h_kc3_R4 +
      2 * shear_modulus_shell * kc3^2 * hdd_kc3_R4
    M[6, 6] <- omega^2 * density_gas * j_k4_R4
  } else {
    M <- matrix(0 + 0i, nrow = 10, ncol = 10)
    F <- matrix(0 + 0i, nrow = 10, ncol = 1)
    mm <- -order * (order + 1)
    pp <- order * (order + 1)

    M[1, 1] <- k1 * radius_viscous * hd_k1_R2
    M[1, 2] <- -kc2 * radius_viscous * jd_kc2_R2
    M[1, 3] <- -kc2 * radius_viscous * hd_kc2_R2
    M[1, 4] <- pp * j_ks2_R2
    M[1, 5] <- pp * h_ks2_R2
    F[1, 1] <- -k1 * radius_viscous * ((1i)^order) * (2 * order + 1) * jd_k1_R2

    M[2, 1] <- -1i * omega * density_sw * radius_viscous^2 * h_k1_R2
    M[2, 2] <- -shear_viscosity_viscous * radius_viscous^2 *
      (2 * kc2^2 - ks2^2) * j_kc2_R2 -
      2 * shear_viscosity_viscous * radius_viscous^2 * kc2^2 * jdd_kc2_R2
    M[2, 3] <- -shear_viscosity_viscous * radius_viscous^2 *
      (2 * kc2^2 - ks2^2) * h_kc2_R2 -
      2 * shear_viscosity_viscous * radius_viscous^2 * kc2^2 * hdd_kc2_R2
    M[2, 4] <- -2 * shear_viscosity_viscous * mm * ks2 * radius_viscous *
      jd_ks2_R2 + 2 * shear_viscosity_viscous * mm * j_ks2_R2
    M[2, 5] <- -2 * shear_viscosity_viscous * mm * ks2 * radius_viscous *
      hd_ks2_R2 + 2 * shear_viscosity_viscous * mm * h_ks2_R2
    F[2, 1] <- 1i * omega * density_sw * radius_viscous^2 * ((1i)^order) *
      (2 * order + 1) * j_k1_R2

    M[3, 2] <- 2 * kc2 * radius_viscous * jd_kc2_R2 - 2 * j_kc2_R2
    M[3, 3] <- 2 * kc2 * radius_viscous * hd_kc2_R2 - 2 * h_kc2_R2
    M[3, 4] <- -ks2^2 * radius_viscous^2 * jdd_ks2_R2 + 2 * j_ks2_R2 +
      mm * j_ks2_R2
    M[3, 5] <- -ks2^2 * radius_viscous^2 * hdd_ks2_R2 + 2 * h_ks2_R2 +
      mm * h_ks2_R2

    M[4, 2] <- kc2 * radius_shell * jd_kc2_R3
    M[4, 3] <- kc2 * radius_shell * hd_kc2_R3
    M[4, 4] <- mm * j_ks2_R3
    M[4, 5] <- mm * h_ks2_R3
    M[4, 6] <- -kc3 * radius_shell * jd_kc3_R3
    M[4, 7] <- -kc3 * radius_shell * hd_kc3_R3
    M[4, 8] <- pp * j_ks3_R3
    M[4, 9] <- pp * h_ks3_R3

    M[5, 2] <- -1i * omega * shear_viscosity_viscous *
      (2 * kc2^2 - ks2^2) * radius_shell^2 * j_kc2_R3 -
      1i * omega * 2 * shear_viscosity_viscous * kc2^2 * radius_shell^2 *
      jdd_kc2_R3
    M[5, 3] <- -1i * omega * shear_viscosity_viscous *
      (2 * kc2^2 - ks2^2) * radius_shell^2 * h_kc2_R3 -
      1i * omega * 2 * shear_viscosity_viscous * kc2^2 * radius_shell^2 *
      hdd_kc2_R3
    M[5, 4] <- -1i * omega * 2 * shear_viscosity_viscous * radius_shell *
      mm * ks2 * jd_ks2_R3 + 1i * omega * 2 * shear_viscosity_viscous *
      mm * j_ks2_R3
    M[5, 5] <- -1i * omega * 2 * shear_viscosity_viscous * radius_shell *
      mm * ks2 * hd_ks2_R3 + 1i * omega * 2 * shear_viscosity_viscous *
      mm * h_ks2_R3
    M[5, 6] <- -shear_modulus_shell * (2 * kc3^2 - ks3^2) * radius_shell^2 *
      j_kc3_R3 - 2 * shear_modulus_shell * kc3^2 * radius_shell^2 * jdd_kc3_R3
    M[5, 7] <- -shear_modulus_shell * (2 * kc3^2 - ks3^2) * radius_shell^2 *
      h_kc3_R3 - 2 * shear_modulus_shell * kc3^2 * radius_shell^2 * hdd_kc3_R3
    M[5, 8] <- -2 * shear_modulus_shell * mm * ks3 * radius_shell * jd_ks3_R3 +
      2 * shear_modulus_shell * mm * j_ks3_R3
    M[5, 9] <- -2 * shear_modulus_shell * mm * ks3 * radius_shell * hd_ks3_R3 +
      2 * shear_modulus_shell * mm * h_ks3_R3

    M[6, 2] <- -1i * omega * 2 * shear_viscosity_viscous * kc2 * radius_shell *
      jd_kc2_R3 + 1i * omega * 2 * shear_viscosity_viscous * j_kc2_R3
    M[6, 3] <- -1i * omega * 2 * shear_viscosity_viscous * kc2 * radius_shell *
      hd_kc2_R3 + 1i * omega * 2 * shear_viscosity_viscous * h_kc2_R3
    M[6, 4] <- 1i * omega * shear_viscosity_viscous * ks2^2 * radius_shell^2 *
      jdd_ks2_R3 - 1i * omega * 2 * shear_viscosity_viscous * j_ks2_R3 -
      1i * omega * shear_viscosity_viscous * mm * j_ks2_R3
    M[6, 5] <- 1i * omega * shear_viscosity_viscous * ks2^2 * radius_shell^2 *
      hdd_ks2_R3 - 1i * omega * 2 * shear_viscosity_viscous * h_ks2_R3 -
      1i * omega * shear_viscosity_viscous * mm * h_ks2_R3
    M[6, 6] <- -2 * shear_modulus_shell * kc3 * radius_shell * jd_kc3_R3 +
      2 * shear_modulus_shell * j_kc3_R3
    M[6, 7] <- -2 * shear_modulus_shell * kc3 * radius_shell * hd_kc3_R3 +
      2 * shear_modulus_shell * h_kc3_R3
    M[6, 8] <- shear_modulus_shell * ks3^2 * radius_shell^2 * jdd_ks3_R3 -
      2 * shear_modulus_shell * j_ks3_R3 - shear_modulus_shell * mm * j_ks3_R3
    M[6, 9] <- shear_modulus_shell * ks3^2 * radius_shell^2 * hdd_ks3_R3 -
      2 * shear_modulus_shell * h_ks3_R3 - shear_modulus_shell * mm * h_ks3_R3

    M[7, 2] <- j_kc2_R3
    M[7, 3] <- h_kc2_R3
    M[7, 4] <- -j_ks2_R3 - ks2 * radius_shell * jd_ks2_R3
    M[7, 5] <- -h_ks2_R3 - ks2 * radius_shell * hd_ks2_R3
    M[7, 6] <- -j_kc3_R3
    M[7, 7] <- -h_kc3_R3
    M[7, 8] <- j_ks3_R3 + ks3 * radius_shell * jd_ks3_R3
    M[7, 9] <- h_ks3_R3 + ks3 * radius_shell * hd_ks3_R3

    M[8, 6] <- kc3 * radius_gas * jd_kc3_R4
    M[8, 7] <- kc3 * radius_gas * hd_kc3_R4
    M[8, 8] <- mm * j_ks3_R4
    M[8, 9] <- mm * h_ks3_R4
    M[8, 10] <- -k4 * radius_gas * jd_k4_R4

    M[9, 6] <- shear_modulus_shell * (2 * kc3^2 - ks3^2) * radius_gas^2 *
      j_kc3_R4 + 2 * shear_modulus_shell * kc3^2 * radius_gas^2 * jdd_kc3_R4
    M[9, 7] <- shear_modulus_shell * (2 * kc3^2 - ks3^2) * radius_gas^2 *
      h_kc3_R4 + 2 * shear_modulus_shell * kc3^2 * radius_gas^2 * hdd_kc3_R4
    M[9, 8] <- 2 * shear_modulus_shell * mm * ks3 * radius_gas * jd_ks3_R4 -
      2 * shear_modulus_shell * mm * j_ks3_R4
    M[9, 9] <- 2 * shear_modulus_shell * mm * ks3 * radius_gas * hd_ks3_R4 -
      2 * shear_modulus_shell * mm * h_ks3_R4
    M[9, 10] <- omega^2 * density_gas * radius_gas^2 * j_k4_R4

    M[10, 6] <- 2 * shear_modulus_shell * kc3 * radius_gas * jd_kc3_R4 -
      2 * shear_modulus_shell * j_kc3_R4
    M[10, 7] <- 2 * shear_modulus_shell * kc3 * radius_gas * hd_kc3_R4 -
      2 * shear_modulus_shell * h_kc3_R4
    M[10, 8] <- -shear_modulus_shell * ks3^2 * radius_gas^2 * jdd_ks3_R4 +
      2 * shear_modulus_shell * j_ks3_R4 + shear_modulus_shell * mm * j_ks3_R4
    M[10, 9] <- -shear_modulus_shell * ks3^2 * radius_gas^2 * hdd_ks3_R4 +
      2 * shear_modulus_shell * h_ks3_R4 + shear_modulus_shell * mm * h_ks3_R4
  }

  list(M = M, F = F)
}

#' @noRd
.vesms_single_mode <- function(order,
                               omega,
                               density_sw,
                               density_viscous,
                               density_shell,
                               density_gas,
                               sound_speed_sw,
                               sound_speed_viscous,
                               sound_speed_gas,
                               radius_viscous,
                               radius_shell,
                               radius_gas,
                               shear_viscosity_viscous,
                               bulk_viscosity_viscous,
                               shear_modulus_shell,
                               lambda_shell) {
  system <- .vesms_build_mode_system(
    order = order,
    omega = omega,
    density_sw = density_sw,
    density_viscous = density_viscous,
    density_shell = density_shell,
    density_gas = density_gas,
    sound_speed_sw = sound_speed_sw,
    sound_speed_viscous = sound_speed_viscous,
    sound_speed_gas = sound_speed_gas,
    radius_viscous = radius_viscous,
    radius_shell = radius_shell,
    radius_gas = radius_gas,
    shear_viscosity_viscous = shear_viscosity_viscous,
    bulk_viscosity_viscous = bulk_viscosity_viscous,
    shear_modulus_shell = shear_modulus_shell,
    lambda_shell = lambda_shell
  )

  sol <- .vesms_solve_system(system$M, system$F)
  sol[1, 1]
}

#' Initialize ESS-class object for the viscous-elastic spherical model.
#' @param object ESS-class object.
#' @param frequency Frequency vector (Hz).
#' @param sound_speed_sw Seawater sound speed (m/s).
#' @param density_sw Seawater density (kg/m^3).
#' @param sound_speed_viscous Compressional sound speed in the viscous layer.
#' @param density_viscous Density of the viscous layer.
#' @param shear_viscosity_viscous Shear viscosity of the viscous layer.
#' @param bulk_viscosity_viscous Optional bulk viscosity of the viscous layer.
#' @param radius_viscous Optional outer radius of the viscous layer.
#' @param viscous_thickness Optional viscous-layer thickness.
#' @param m_limit Optional modal truncation limit.
#' @noRd
vesms_initialize <- function(object,
                             frequency,
                             sound_speed_sw = .SEAWATER_SOUND_SPEED_DEFAULT,
                             density_sw = .SEAWATER_DENSITY_DEFAULT,
                             sound_speed_viscous,
                             density_viscous,
                             shear_viscosity_viscous,
                             bulk_viscosity_viscous = NULL,
                             radius_viscous = NULL,
                             viscous_thickness = NULL,
                             m_limit = NULL) {
  if (!methods::is(object, "ESS")) {
    stop(
      "VESMS currently requires an elastic-shelled scatterer ('ESS'). Input ",
      "scatterer is type '", class(object), "'."
    )
  }

  if (!.vesms_is_spherical(object)) {
    stop("VESMS currently requires a spherical ESS geometry.")
  }

  if (!is.null(radius_viscous) && !is.null(viscous_thickness)) {
    stop(
      "Specify at most one of 'radius_viscous' or 'viscous_thickness' for ",
      "VESMS."
    )
  }

  shell <- .extract_material_props(
    extract(object, "shell"),
    sound_speed_sw,
    density_sw
  )
  fluid <- .extract_material_props(
    extract(object, "fluid"),
    sound_speed_sw,
    density_sw
  )

  shell$radius <- extract(object, "shell")$radius
  shell$shell_thickness <- extract(object, "shell")$shell_thickness
  fluid$radius <- extract(object, "fluid")$radius

  if (is.null(shell$radius) || is.null(fluid$radius) ||
      !is.finite(shell$radius) || !is.finite(fluid$radius) ||
      fluid$radius >= shell$radius) {
    stop(
      "VESMS requires a valid shell radius larger than the inner gas radius."
    )
  }

  if (is.null(shell$density) || !is.finite(shell$density) ||
      is.null(shell$G) || !is.finite(shell$G)) {
    stop(
      "VESMS requires the ESS shell to have finite density and shear modulus."
    )
  }

  if (is.null(shell$lambda) || !is.finite(shell$lambda)) {
    shell$lambda <- tryCatch(
      lame(K = shell$K, E = shell$E, G = shell$G, nu = shell$nu),
      error = function(e) NA_real_
    )
  }

  if (!is.finite(shell$lambda)) {
    stop(
      "VESMS requires Lam", "\u00E9", "'s first parameter for the shell. ",
      "Provide sufficient shell elastic properties when building the ESS ",
      "object."
    )
  }

  if (is.null(fluid$density) || !is.finite(fluid$density) ||
      is.null(fluid$sound_speed) || !is.finite(fluid$sound_speed)) {
    stop(
      "VESMS requires the ESS inner fluid slot to represent a gas core with ",
      "finite density and sound speed."
    )
  }

  if (missing(sound_speed_viscous) ||
      missing(density_viscous) ||
      missing(shear_viscosity_viscous)) {
    stop(
      "VESMS requires 'sound_speed_viscous', 'density_viscous', and ",
      "'shear_viscosity_viscous'."
    )
  }

  if (!is.finite(sound_speed_viscous) || sound_speed_viscous <= 0 ||
      !is.finite(density_viscous) || density_viscous <= 0 ||
      !is.finite(shear_viscosity_viscous) || shear_viscosity_viscous <= 0) {
    stop(
      "VESMS requires positive finite viscous-layer sound speed, density, and ",
      "shear viscosity."
    )
  }

  bulk_viscosity_viscous <- if (is.null(bulk_viscosity_viscous)) {
    shear_viscosity_viscous
  } else {
    bulk_viscosity_viscous
  }

  if (!is.finite(bulk_viscosity_viscous) || bulk_viscosity_viscous < 0) {
    stop("VESMS requires a finite non-negative bulk viscosity.")
  }

  radius_viscous <- if (!is.null(radius_viscous)) {
    .vesms_validate_outer_radius(radius_viscous, shell$radius)
  } else if (!is.null(viscous_thickness)) {
    .vesms_validate_outer_radius(shell$radius + viscous_thickness, shell$radius)
  } else {
    .vesms_validate_outer_radius(
      .vesms_neutral_radius(
        radius_gas = fluid$radius,
        density_sw = density_sw,
        density_viscous = density_viscous,
        density_gas = fluid$density
      ),
      shell$radius
    )
  }

  m_limit <- if (is.null(m_limit)) {
    pmax(2L, round(wavenumber(frequency, sound_speed_sw) * radius_viscous) + 10L)
  } else {
    as.integer(m_limit)
  }

  if (length(m_limit) == 1L) {
    m_limit <- rep(m_limit, length(frequency))
  }

  if (length(m_limit) != length(frequency) ||
      any(!is.finite(m_limit)) ||
      any(m_limit < 0)) {
    stop(
      "VESMS requires 'm_limit' to be a non-negative scalar or a vector with ",
      "one value per frequency."
    )
  }

  acoustics <- data.frame(
    frequency = frequency,
    k_sw = wavenumber(frequency, sound_speed_sw),
    m_limit = m_limit
  )

  methods::slot(object, "model_parameters")$VESMS <- list(
    parameters = list(acoustics = acoustics),
    medium = data.frame(
      sound_speed = sound_speed_sw,
      density = density_sw
    ),
    viscous = list(
      sound_speed = sound_speed_viscous,
      density = density_viscous,
      shear_viscosity = shear_viscosity_viscous,
      bulk_viscosity = bulk_viscosity_viscous,
      radius = radius_viscous
    ),
    shell = shell,
    fluid = fluid
  )

  methods::slot(object, "model")$VESMS <- data.frame(
    frequency = frequency,
    sigma_bs = rep(NA_real_, length(frequency))
  )

  object
}

#' Viscous-elastic spherical scattering model.
#' @param object ESS-class object.
#' @noRd
VESMS <- function(object) {
  model <- extract(object, "model_parameters")$VESMS
  acoustics <- model$parameters$acoustics
  medium <- model$medium
  viscous <- model$viscous
  shell <- model$shell
  fluid <- model$fluid
  r_inf <- 1

  f_bs <- vapply(seq_len(nrow(acoustics)), function(i) {
    omega <- 2 * pi * acoustics$frequency[i]
    terms <- vapply(0:acoustics$m_limit[i], function(m) {
      A1 <- .vesms_single_mode(
        order = m,
        omega = omega,
        density_sw = medium$density,
        density_viscous = viscous$density,
        density_shell = shell$density,
        density_gas = fluid$density,
        sound_speed_sw = medium$sound_speed,
        sound_speed_viscous = viscous$sound_speed,
        sound_speed_gas = fluid$sound_speed,
        radius_viscous = viscous$radius,
        radius_shell = shell$radius,
        radius_gas = fluid$radius,
        shear_viscosity_viscous = viscous$shear_viscosity,
        bulk_viscosity_viscous = viscous$bulk_viscosity,
        shear_modulus_shell = shell$G,
        lambda_shell = shell$lambda
      )
      A1 * ((-1)^m) * .vesms_sph_h(m, acoustics$k_sw[i] * r_inf)
    }, FUN.VALUE = complex(1))
    sum(terms)
  }, FUN.VALUE = complex(1))

  sigma_bs <- .sigma_bs(f_bs)

  methods::slot(object, "model")$VESMS <- data.frame(
    frequency = acoustics$frequency,
    ka_viscous = acoustics$k_sw * viscous$radius,
    ka_shell = acoustics$k_sw * shell$radius,
    ka_gas = acoustics$k_sw * fluid$radius,
    f_bs = f_bs,
    sigma_bs = sigma_bs,
    TS = db(sigma_bs)
  )

  object
}
