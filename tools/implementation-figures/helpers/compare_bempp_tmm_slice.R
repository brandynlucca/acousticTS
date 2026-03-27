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

# Convert world-frame spherical angles with the target axis along +x into the
# body-fixed angles expected by the current axisymmetric TMM post-processing
# helpers.
world_to_body_fixed_xaxis <- function(theta, phi) {
  x <- sin(theta) * cos(phi)
  y <- sin(theta) * sin(phi)
  z <- cos(theta)

  theta_body <- acos(pmax(-1, pmin(1, x)))
  phi_body <- atan2(z, y)
  phi_body <- ifelse(phi_body < 0, phi_body + 2 * pi, phi_body)

  # Azimuth is undefined exactly on the symmetry axis, so pin it to zero there.
  on_axis <- abs(sin(theta_body)) < 1e-10
  phi_body[on_axis] <- 0

  data.frame(theta = theta_body, phi = phi_body)
}

compare_bempp_tmm_slice <- function(object,
                                    bempp_path,
                                    frame = c("body_fixed", "world_xy_xaxis"),
                                    incidence_deg = NULL,
                                    theta_body = pi / 2,
                                    phi_body = pi / 2,
                                    theta_scatter = pi / 2) {
  frame <- match.arg(frame)
  bempp <- read_bempp_slice(bempp_path)

  if (identical(frame, "world_xy_xaxis")) {
    if (!is.numeric(incidence_deg) || length(incidence_deg) != 1 || !is.finite(incidence_deg)) {
      stop(
        "'incidence_deg' must be a single finite numeric angle when ",
        "'frame = \"world_xy_xaxis\"'.",
        call. = FALSE
      )
    }
    incident_angles <- world_to_body_fixed_xaxis(
      theta = pi / 2,
      phi = incidence_deg * pi / 180
    )
  } else {
    incident_angles <- data.frame(theta = theta_body, phi = phi_body)
  }

  tmm <- do.call(
    rbind,
    lapply(
      seq_len(nrow(bempp)),
      function(i) {
        if (identical(frame, "world_xy_xaxis")) {
          scatter_angles <- world_to_body_fixed_xaxis(
            theta = pi / 2,
            phi = bempp$phi_scatter_rad[i]
          )
        } else {
          scatter_angles <- data.frame(theta = theta_scatter, phi = bempp$phi_scatter_rad[i])
        }

        scat <- tmm_scattering(
          object,
          theta_body = incident_angles$theta,
          phi_body = incident_angles$phi,
          theta_scatter = scatter_angles$theta,
          phi_scatter = scatter_angles$phi
        )
        data.frame(
          angle_deg = bempp$angle_deg[i],
          phi_scatter_rad = bempp$phi_scatter_rad[i],
          theta_body_tmm = incident_angles$theta,
          phi_body_tmm = incident_angles$phi,
          theta_scatter_tmm = scatter_angles$theta,
          phi_scatter_tmm = scatter_angles$phi,
          mod_f_tmm = Mod(scat$f_scat),
          sigma_scat_dB = scat$sigma_scat_dB
        )
      }
    )
  )

  cmp <- merge(bempp, tmm, by = c("angle_deg", "phi_scatter_rad"), all = FALSE)
  cmp$delta_amp_db <- 20 * log10(
    pmax(cmp$mod_f_tmm, .Machine$double.eps) /
      pmax(cmp$mod_f_bempp, .Machine$double.eps)
  )

  list(
    summary = data.frame(
      max_abs_delta_amp_db = max(abs(cmp$delta_amp_db)),
      mean_abs_delta_amp_db = mean(abs(cmp$delta_amp_db))
    ),
    detail = cmp
  )
}
