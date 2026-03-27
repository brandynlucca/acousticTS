library(devtools)

load_all("c:/Users/Brandyn/Desktop/acousticTS_upd", quiet = TRUE)

read_bempp_slice <- function(path) {
  valid_cols <- c("angle_deg", "phi_scatter_rad", "mod_f")
  parse_try <- function(skip = 0) {
    tryCatch(
      read.delim(
        path,
        sep = "\t",
        skip = skip,
        fileEncoding = "UTF-16LE",
        check.names = FALSE
      ),
      error = function(e) NULL
    )
  }

  tab <- parse_try(0)
  if (is.null(tab) || !all(valid_cols %in% names(tab))) {
    tab <- parse_try(2)
  }
  if (is.null(tab) || !all(valid_cols %in% names(tab))) {
    stop("Could not parse BEMPP slice file: ", path, call. = FALSE)
  }

  data.frame(
    angle_deg = tab$angle_deg,
    phi_scatter_rad = tab$phi_scatter_rad,
    mod_f_bempp = tab$mod_f
  )
}

sphere_modal_coefficients <- function(boundary,
                                      radius,
                                      frequency,
                                      sound_speed_sw,
                                      density_sw,
                                      sound_speed_body = NULL,
                                      density_body = NULL) {
  k_sw <- 2 * pi * frequency / sound_speed_sw
  ka <- k_sw * radius
  m_limit <- round(ka + 20)
  m <- 0:m_limit

  if (boundary == "pressure_release") {
    return(-(js(m, ka) / hs(m, ka)))
  }

  if (boundary == "fixed_rigid") {
    return(-(jsd(m, ka) / hsd(m, ka)))
  }

  if (boundary == "liquid_filled") {
    h31 <- sound_speed_body / sound_speed_sw
    g31 <- density_body / density_sw
    k_body_a <- ka / h31
    gh <- g31 * h31
    cm_num <- (
      (jsd(m, k_body_a) * ys(m, ka)) / (js(m, k_body_a) * jsd(m, ka)) -
        (gh * ysd(m, ka) / jsd(m, ka))
    )
    cm_denom <- (
      (jsd(m, k_body_a) * js(m, ka)) / (js(m, k_body_a) * jsd(m, ka)) -
        gh
    )
    cm <- cm_num / cm_denom
    return(-1 / (1 + 1i * cm))
  }

  stop("Unsupported sphere comparison boundary.", call. = FALSE)
}

exact_sphere_amplitude <- function(boundary,
                                   radius,
                                   frequency,
                                   sound_speed_sw,
                                   density_sw,
                                   incidence_deg,
                                   phi_scatter_rad,
                                   sound_speed_body = NULL,
                                   density_body = NULL) {
  k_sw <- 2 * pi * frequency / sound_speed_sw
  coeffs <- sphere_modal_coefficients(
    boundary = boundary,
    radius = radius,
    frequency = frequency,
    sound_speed_sw = sound_speed_sw,
    density_sw = density_sw,
    sound_speed_body = sound_speed_body,
    density_body = density_body
  )

  m <- 0:(length(coeffs) - 1)
  incidence_rad <- incidence_deg * pi / 180
  d <- c(cos(incidence_rad), sin(incidence_rad), 0)
  recv_dirs <- cbind(
    cos(phi_scatter_rad),
    sin(phi_scatter_rad),
    rep(0, length(phi_scatter_rad))
  )
  gamma <- acos(pmax(-1, pmin(1, recv_dirs[, 1] * d[1] + recv_dirs[, 2] * d[2])))
  P <- Pn(m, cos(gamma))

  as.vector(-(1i / k_sw) * crossprod((2 * m + 1) * coeffs, P))
}

compare_bempp_exact_sphere <- function(path,
                                       boundary,
                                       radius,
                                       frequency,
                                       sound_speed_sw,
                                       density_sw,
                                       incidence_deg,
                                       sound_speed_body = NULL,
                                       density_body = NULL) {
  bempp <- read_bempp_slice(path)
  f_exact <- exact_sphere_amplitude(
    boundary = boundary,
    radius = radius,
    frequency = frequency,
    sound_speed_sw = sound_speed_sw,
    density_sw = density_sw,
    incidence_deg = incidence_deg,
    phi_scatter_rad = bempp$phi_scatter_rad,
    sound_speed_body = sound_speed_body,
    density_body = density_body
  )

  cmp <- transform(
    bempp,
    mod_f_exact = Mod(f_exact)
  )
  cmp$delta_amp_db <- 20 * log10(
    pmax(cmp$mod_f_exact, .Machine$double.eps) /
      pmax(cmp$mod_f_bempp, .Machine$double.eps)
  )

  list(
    summary = data.frame(
      boundary = boundary,
      max_abs_delta_amp_db = max(abs(cmp$delta_amp_db)),
      mean_abs_delta_amp_db = mean(abs(cmp$delta_amp_db))
    ),
    detail = cmp
  )
}
