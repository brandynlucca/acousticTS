################################################################################
# Phase-compensated distorted wave Born approximation (PCDWBA)
################################################################################
#' Phase-compensated distorted wave Born approximation (PCDWBA)
#'
#' @description
#' Computes backscatter from weakly scattering elongated fluid-like targets
#' using the phase-compensated distorted wave Born approximation (PCDWBA)
#' described by Chu and Ye (1999). This model is particularly suited to
#' uniformly bent or gently curved cylindrical bodies with tapered ends, but it
#' can also operate on arbitrary fluid-like centerline profiles stored in
#' `FLS` objects.
#'
#' @section Usage:
#' This model is accessed via:
#' \preformatted{
#' target_strength(
#'   ...,
#'   model = "PCDWBA",
#'   sound_speed_sw,
#'   density_sw,
#'   radius_curvature,
#'   radius_curvature_ratio
#' )
#' }
#'
#' @section Arguments:
#' \describe{
#'   \item{\code{sound_speed_sw}}{Seawater sound speed (\eqn{m~s^{-1}}).}
#'   \item{\code{density_sw}}{Seawater density (\eqn{kg~m^{-3}}).}
#'   \item{\code{radius_curvature}}{Optional radius of curvature (\eqn{m}) used
#'   to rebuild a uniformly bent centerline for canonical cylinders.}
#'   \item{\code{radius_curvature_ratio}}{Optional radius of curvature divided
#'   by body length (\eqn{\rho_c / L}) used to rebuild a uniformly bent
#'   centerline for canonical cylinders.}
#' }
#'
#' @references
#' Chu, D., and Ye, Z. (1999). A phase-compensated distorted wave Born
#' approximation representation of the bistatic scattering by weakly scattering
#' objects: Application to zooplankton. \emph{The Journal of the Acoustical
#' Society of America}, 106, 1732-1743.
#'
#' Stanton, T.K. (1989). Sound scattering by cylinders of finite length. III.
#' Deformed cylinders. \emph{The Journal of the Acoustical Society of America},
#' 86, 691-705.
#'
#' @name PCDWBA
#' @aliases pcdwba PCDWBA
#' @docType data
#' @keywords models acoustics
NULL

#' @noRd
.pcdwba_axis_row <- function(rpos, candidates, default = NULL) {
  if (!is.null(rownames(rpos))) {
    idx <- match(candidates, rownames(rpos), nomatch = 0L)
    idx <- idx[idx > 0L]
    if (length(idx) > 0L) {
      return(idx[1])
    }
  }

  if (is.null(default)) {
    return(NULL)
  }

  default
}

#' @noRd
.pcdwba_is_straight <- function(body) {
  rpos <- body$rpos
  z_idx <- .pcdwba_axis_row(rpos, c("z", "z_body"), default = NULL)

  if (is.null(z_idx)) {
    return(TRUE)
  }

  all(abs(rpos[z_idx, ]) <= sqrt(.Machine$double.eps))
}

#' @noRd
.pcdwba_node_property <- function(x, n_nodes, label) {
  if (length(x) == 1L) {
    return(rep(as.numeric(x), n_nodes))
  }

  if (length(x) == n_nodes) {
    return(as.numeric(x))
  }

  if (length(x) == n_nodes - 1L) {
    return(c(as.numeric(x[1]), as.numeric(x)))
  }

  stop(
    "PCDWBA requires '", label, "' to be scalar, nodewise, or segmentwise. ",
    "Received length ", length(x), " for ", n_nodes, " nodes."
  )
}

#' @noRd
.pcdwba_regular_geometry <- function(n_nodes,
                                     radius_curvature_ratio,
                                     taper_order = NA_real_) {
  gamma <- 0.5 / radius_curvature_ratio
  norm_ratio <- radius_curvature_ratio * 2
  z <- seq(-1, 1, length.out = n_nodes)
  taper <- if (is.na(taper_order)) {
    rep(1, n_nodes)
  } else {
    sqrt(pmax(1 - z^taper_order, 0))
  }

  z_curved <- sin(gamma) * z
  x_curved <- 1 - sqrt(pmax(1 - z_curved^2, 0))
  x_norm <- x_curved * norm_ratio
  z_norm <- z_curved * norm_ratio

  gamma_t <- atan2(z_norm, x_norm)
  dz <- diff(z_norm)
  dx <- diff(x_norm) + .Machine$double.eps
  alpha_t <- c(atan(dz / dx), atan(dz[length(dz)] / dx[length(dx)]))
  theta_tilt <- ifelse(alpha_t >= 0, alpha_t - pi / 2, alpha_t + pi / 2)
  dr1 <- sqrt(dx^2 + dz^2)
  dr_norm <- c(dr1[1], dr1)
  r_pos <- sqrt(x_norm^2 + z_norm^2)

  list(
    taper = as.numeric(taper),
    theta_tilt = as.numeric(theta_tilt),
    gamma_t = as.numeric(gamma_t),
    dr_norm = as.numeric(dr_norm),
    r_pos = as.numeric(r_pos)
  )
}

#' @noRd
.pcdwba_profile_geometry <- function(body) {
  rpos <- body$rpos
  n_nodes <- ncol(rpos)

  x_idx <- .pcdwba_axis_row(rpos, c("x", "x_body"), default = 1L)
  z_idx <- .pcdwba_axis_row(rpos, c("z", "z_body"), default = NULL)

  x_long <- as.numeric(rpos[x_idx, ])
  z_trans <- if (is.null(z_idx)) {
    rep(0, n_nodes)
  } else {
    as.numeric(rpos[z_idx, ])
  }

  dx_phys <- diff(x_long)
  dz_phys <- diff(z_trans)
  length_body <- sum(sqrt(dx_phys^2 + dz_phys^2))

  if (!is.finite(length_body) || length_body <= 0) {
    stop("PCDWBA could not resolve a positive body length from the scatterer geometry.")
  }

  longitudinal <- x_long - mean(range(x_long))
  transverse <- (-z_trans) - mean(-z_trans)

  x_norm <- 2 * transverse / length_body
  z_norm <- 2 * longitudinal / length_body
  taper <- .dwba_body_radius(body) / max(.dwba_body_radius(body))

  gamma_t <- atan2(z_norm, x_norm)
  dz <- diff(z_norm)
  dx <- diff(x_norm) + .Machine$double.eps
  alpha_t <- c(atan(dz / dx), atan(dz[length(dz)] / dx[length(dx)]))
  theta_tilt <- ifelse(alpha_t >= 0, alpha_t - pi / 2, alpha_t + pi / 2)
  dr1 <- sqrt(dx^2 + dz^2)
  dr_norm <- c(dr1[1], dr1)
  r_pos <- sqrt(x_norm^2 + z_norm^2)

  list(
    length_body = length_body,
    taper = as.numeric(taper),
    theta_tilt = as.numeric(theta_tilt),
    gamma_t = as.numeric(gamma_t),
    dr_norm = as.numeric(dr_norm),
    r_pos = as.numeric(r_pos)
  )
}

#' @noRd
.pcdwba_prepare_geometry <- function(object,
                                     radius_curvature = NULL,
                                     radius_curvature_ratio = NULL) {
  if (!methods::is(object, "FLS")) {
    stop(
      "PCDWBA requires a fluid-like scatterer ('FLS'). Input scatterer is type '",
      class(object), "'."
    )
  }

  if (!is.null(radius_curvature) && !is.null(radius_curvature_ratio)) {
    stop(
      "Specify at most one of 'radius_curvature' or 'radius_curvature_ratio' ",
      "for PCDWBA."
    )
  }

  object <- .as_dwba_profile(object)
  shape <- extract(object, "shape_parameters")
  body <- extract(object, "body")
  body$radius <- .dwba_body_radius(body)
  methods::slot(object, "body") <- body

  max_radius <- max(body$radius, na.rm = TRUE)
  if (!is.finite(max_radius) || max_radius <= 0) {
    stop("PCDWBA requires a positive body radius.")
  }

  x_idx <- .pcdwba_axis_row(body$rpos, c("x", "x_body"), default = 1L)
  z_idx <- .pcdwba_axis_row(body$rpos, c("z", "z_body"), default = NULL)
  z_pos <- if (is.null(z_idx)) {
    rep(0, ncol(body$rpos))
  } else {
    body$rpos[z_idx, ]
  }
  length_body_profile <- sum(sqrt(diff(body$rpos[x_idx, ])^2 + diff(z_pos)^2))

  if (is.null(length_body_profile) || !is.finite(length_body_profile) || length_body_profile <= 0) {
    length_body_profile <- diff(range(body$rpos[1, ]))
  }

  effective_ratio <- if (!is.null(radius_curvature_ratio)) {
    radius_curvature_ratio
  } else if (!is.null(radius_curvature)) {
    radius_curvature / length_body_profile
  } else if (!is.null(shape$radius_curvature_ratio) && !is.na(shape$radius_curvature_ratio)) {
    shape$radius_curvature_ratio
  } else if (!is.null(body$radius_curvature_ratio) && !is.na(body$radius_curvature_ratio)) {
    body$radius_curvature_ratio
  } else {
    NA_real_
  }

  if (!is.na(effective_ratio) &&
      identical(shape$shape, "Cylinder")) {
    geometry <- .pcdwba_regular_geometry(
      n_nodes = ncol(body$rpos),
      radius_curvature_ratio = effective_ratio,
      taper_order = if ("taper_order" %in% names(shape)) shape$taper_order else NA_real_
    )
    geometry$length_body <- shape$length
  } else {
    if (!is.na(effective_ratio) && .pcdwba_is_straight(body)) {
      object <- brake(object, effective_ratio, mode = "ratio")
      body <- extract(object, "body")
      body$radius <- .dwba_body_radius(body)
      methods::slot(object, "body") <- body
    }
    geometry <- .pcdwba_profile_geometry(extract(object, "body"))
  }

  list(
    object = object,
    geometry = geometry,
    max_radius = max_radius
  )
}

#' Initialize FLS-class object for PCDWBA modeling.
#' @noRd
pcdwba_initialize <- function(object,
                              frequency,
                              sound_speed_sw = 1500,
                              density_sw = 1026,
                              radius_curvature = NULL,
                              radius_curvature_ratio = NULL) {
  prepared <- .pcdwba_prepare_geometry(
    object = object,
    radius_curvature = radius_curvature,
    radius_curvature_ratio = radius_curvature_ratio
  )
  object <- prepared$object
  geometry <- prepared$geometry
  body <- extract(object, "body")
  contrasts <- .derive_contrasts(body, sound_speed_sw, density_sw)
  body$h <- .pcdwba_node_property(
    contrasts$h,
    n_nodes = ncol(body$rpos),
    label = "h"
  )
  body$g <- .pcdwba_node_property(
    contrasts$g,
    n_nodes = ncol(body$rpos),
    label = "g"
  )

  if (anyNA(body$h) || anyNA(body$g)) {
    stop("PCDWBA requires finite density and sound-speed contrasts.")
  }

  Cb <- (1 - body$g * body$h * body$h) / (body$g * body$h * body$h) -
    (body$g - 1) / body$g

  medium_params <- data.frame(
    sound_speed = sound_speed_sw,
    density = density_sw
  )
  .init_model_slots(
    object = object,
    model_name = "PCDWBA",
    frequency = frequency,
    model_parameters = list(
      parameters = list(
        acoustics = .init_acoustics_df(frequency, k_sw = sound_speed_sw)
      ),
      medium = medium_params,
      body = list(
        theta = body$theta,
        h = body$h,
        g = body$g,
        Cb = Cb,
        taper = geometry$taper,
        theta_tilt = geometry$theta_tilt,
        gamma_t = geometry$gamma_t,
        dr_norm = geometry$dr_norm,
        r_pos = geometry$r_pos,
        length_body = geometry$length_body,
        radius = prepared$max_radius
      )
    )
  )
}

#' Phase-compensated distorted wave Born approximation.
#' @noRd
PCDWBA <- function(object) {
  model <- extract(object, "model_parameters")$PCDWBA
  acoustics <- model$parameters$acoustics
  body <- model$body
  theta <- body$theta
  eps <- .Machine$double.eps

  f_norm <- vapply(seq_len(nrow(acoustics)), function(i) {
    k_sw <- acoustics$k_sw[i]
    X2 <- k_sw * body$radius * body$taper / body$h
    arg <- 2 * X2 * abs(cos(body$theta_tilt - theta)) + eps
    term0 <- k_sw * body$length_body * (body$r_pos / body$h)
    term1 <- body$h^2 * body$Cb * body$dr_norm / 4
    phase <- exp(1i * term0 * cos(body$gamma_t - theta))
    term2 <- (X2^2) * (besselJ(arg, 1) / arg) * phase
    sum(term2 * term1) + eps
  }, complex(1))

  f_bs <- f_norm * body$length_body
  sigma_bs <- abs(f_bs)^2

  methods::slot(object, "model")$PCDWBA <- data.frame(
    frequency = acoustics$frequency,
    ka = acoustics$k_sw * body$radius,
    f_bs = f_bs,
    sigma_bs = sigma_bs,
    TS = 10 * log10(sigma_bs)
  )

  object
}
