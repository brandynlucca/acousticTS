################################################################################
# Body-backbone composite model (BBFM)
################################################################################
#' Body-backbone composite model (BBFM)
#'
#' @description
#' Computes a coherent composite backscatter prediction for swimbladder-less
#' fish represented as a weakly scattering flesh body plus an explicit elastic
#' backbone. The flesh contribution is evaluated with the distorted-wave Born
#' approximation (DWBA), while the backbone contribution is evaluated with the
#' elastic cylinder modal-series solution (ECMS). The backbone amplitude is then
#' translated into the stored body coordinate frame using a phase factor based
#' on the backbone centroid before the two complex amplitudes are summed.
#'
#' @section Usage:
#' This model is accessed via:
#' \preformatted{
#' target_strength(
#'   ...,
#'   model = "BBFM",
#'   sound_speed_sw,
#'   density_sw,
#'   m_limit
#' )
#' }
#'
#' @section Arguments:
#' \describe{
#'   \item{\code{sound_speed_sw}}{Seawater sound speed (\eqn{m~s^{-1}}).}
#'   \item{\code{density_sw}}{Seawater density (\eqn{kg~m^{-3}}).}
#'   \item{\code{m_limit}}{Optional modal truncation limit passed to the
#'   backbone ECMS solve.}
#' }
#'
#' @section Theory:
#' The model follows the same structural decomposition used for
#' swimbladder-less mackerel by Gorska, Ona, and Korneliussen (2005): the flesh
#' is treated as a weakly scattering fluid-like body, and the backbone is
#' treated as an elastic cylindrical structure. In this implementation, the
#' total complex backscattering amplitude is
#' \deqn{
#'   f_{bs} = f_{\mathrm{flesh}} +
#'   f_{\mathrm{backbone}}
#'   \exp\left\{ 2 i k_{sw} v_c \right\},
#' }
#' where \eqn{f_{\mathrm{flesh}}} is obtained from \code{\link{DWBA}},
#' \eqn{f_{\mathrm{backbone}}} is obtained from \code{\link{ECMS}}, and
#' \eqn{v_c = x_c \cos \theta + z_c \sin \theta} is the projection of the
#' stored backbone centroid onto the backscatter direction.
#'
#' @references
#' Gorska, N., Ona, E., and Korneliussen, R. (2005). Acoustic backscattering by
#' Atlantic mackerel as being representative of fish that lack a swimbladder.
#' Backscattering by individual fish. \emph{ICES Journal of Marine Science},
#' 62: 984-995.
#'
#' Stanton, T.K. (1988). Sound scattering by cylinders of finite length. II.
#' Elastic cylinders. \emph{The Journal of the Acoustical Society of America},
#' 83: 64-67.
#'
#' Stanton, T.K. (1989). Sound scattering by cylinders of finite length. III.
#' Deformed cylinders. \emph{The Journal of the Acoustical Society of America},
#' 86: 691-705.
#'
#' Stanton, T.K., Chu, D., and Wiebe, P.H. (1998). Sound scattering by several
#' zooplankton groups. II. Scattering models. \emph{The Journal of the
#' Acoustical Society of America}, 103: 236-253.
#'
#' @name BBFM
#' @aliases bbfm BBFM
#' @docType data
#' @keywords models acoustics
NULL

#' @noRd
.bbfm_make_fls <- function(component, shape_parameters, id = "UID") {
  methods::new("FLS",
    metadata = list(ID = id),
    model_parameters = list(),
    model = list(),
    body = component,
    shape_parameters = shape_parameters
  )
}

#' @noRd
.bbfm_centroid_projection <- function(component) {
  rpos <- component$rpos
  x_center <- if ("x" %in% rownames(rpos)) {
    mean(range(rpos["x", ], na.rm = TRUE))
  } else {
    0
  }
  z_center <- if ("z" %in% rownames(rpos)) {
    mean(range(rpos["z", ], na.rm = TRUE))
  } else {
    0
  }

  x_center * cos(component$theta) + z_center * sin(component$theta)
}

#' Initialize BBF-class object for BBFM modeling.
#' @noRd
bbfm_initialize <- function(object,
                            frequency,
                            sound_speed_sw = .SEAWATER_SOUND_SPEED_DEFAULT,
                            density_sw = .SEAWATER_DENSITY_DEFAULT,
                            m_limit = NULL) {
  if (!methods::is(object, "BBF")) {
    stop(
      "BBFM requires a 'BBF' composite scatterer. Input scatterer is type '",
      class(object), "'.",
      call. = FALSE
    )
  }

  metadata <- extract(object, "metadata")
  shape <- extract(object, "shape_parameters")
  body <- extract(object, "body")
  backbone <- extract(object, "backbone")

  if (is.null(body$density) && is.null(body$g)) {
    stop(
      "BBFM requires the flesh body to carry either absolute density ",
      "('density_body') or density contrast ('g_body').",
      call. = FALSE
    )
  }
  if (is.null(body$sound_speed) && is.null(body$h)) {
    stop(
      "BBFM requires the flesh body to carry either absolute sound speed ",
      "('sound_speed_body') or sound-speed contrast ('h_body').",
      call. = FALSE
    )
  }
  if (is.null(backbone$density) ||
      is.null(backbone$sound_speed_longitudinal) ||
      is.null(backbone$sound_speed_transversal)) {
    stop(
      "BBFM requires backbone density plus longitudinal and transversal wave ",
      "speeds.",
      call. = FALSE
    )
  }

  body_object <- .bbfm_make_fls(
    component = body,
    shape_parameters = shape$body,
    id = paste0(metadata$ID, "_body")
  )
  backbone_object <- .bbfm_make_fls(
    component = backbone,
    shape_parameters = shape$backbone,
    id = paste0(metadata$ID, "_backbone")
  )

  methods::slot(object, "model_parameters")$BBFM <- list(
    parameters = list(
      frequency = frequency,
      sound_speed_sw = sound_speed_sw,
      density_sw = density_sw,
      k_sw = wavenumber(frequency, sound_speed_sw),
      m_limit = m_limit
    ),
    components = list(
      body = body_object,
      backbone = backbone_object
    )
  )

  methods::slot(object, "model")$BBFM <- data.frame(
    frequency = frequency,
    sigma_bs = rep(NA_real_, length(frequency))
  )

  object
}

#' Body-backbone composite model.
#' @noRd
BBFM <- function(object) {
  model <- extract(object, "model_parameters")$BBFM
  params <- model$parameters

  body_object <- target_strength(
    object = model$components$body,
    frequency = params$frequency,
    model = "DWBA",
    sound_speed_sw = params$sound_speed_sw,
    density_sw = params$density_sw
  )

  backbone_component <- extract(model$components$backbone, "body")
  backbone_object <- target_strength(
    object = model$components$backbone,
    frequency = params$frequency,
    model = "ECMS",
    sound_speed_sw = params$sound_speed_sw,
    density_sw = params$density_sw,
    density_body = backbone_component$density,
    sound_speed_longitudinal_body = backbone_component$sound_speed_longitudinal,
    sound_speed_transversal_body = backbone_component$sound_speed_transversal,
    m_limit = params$m_limit
  )

  body_results <- extract(body_object, "model")$DWBA
  backbone_results <- extract(backbone_object, "model")$ECMS
  phase_shift <- exp(
    2i * params$k_sw * .bbfm_centroid_projection(
      extract(backbone_object, "body")
    )
  )
  f_backbone_aligned <- backbone_results$f_bs * phase_shift
  f_bs <- body_results$f_bs + f_backbone_aligned
  sigma_body <- abs(body_results$f_bs)^2
  sigma_backbone <- abs(f_backbone_aligned)^2
  sigma_bs <- abs(f_bs)^2

  methods::slot(object, "model")$BBFM <- data.frame(
    frequency = params$frequency,
    ka_body = body_results$ka,
    ka_backbone = backbone_results$ka,
    f_body = body_results$f_bs,
    f_backbone = backbone_results$f_bs,
    f_backbone_aligned = f_backbone_aligned,
    f_bs = f_bs,
    sigma_body = sigma_body,
    sigma_backbone = sigma_backbone,
    sigma_bs = sigma_bs,
    TS_body = 10 * log10(sigma_body),
    TS_backbone = 10 * log10(sigma_backbone),
    TS = 10 * log10(sigma_bs)
  )

  object
}
