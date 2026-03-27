source("tools/implementation-figures/helpers/common.R")
impl_load_all()
source("tools/implementation-figures/helpers/compare_bempp_tmm_slice.R")

density_sw <- 1026.8
sound_speed_sw <- 1480

read_compare_csv <- function(path) {
  x <- read.csv(path, stringsAsFactors = FALSE)
  x[order(x$angle_deg), , drop = FALSE]
}

plot_case <- function(compare_df,
                      main_title,
                      residual_ylim = NULL) {
  amp_bempp_db <- 20 * log10(pmax(compare_df$mod_f_bempp, .Machine$double.eps))
  amp_tmm_db <- 20 * log10(pmax(compare_df$mod_f_tmm, .Machine$double.eps))

  plot(
    compare_df$angle_deg,
    amp_bempp_db,
    type = "l",
    lwd = 2,
    col = "black",
    xlab = "",
    ylab = "20 log10(|f|) (dB)",
    main = main_title
  )
  lines(
    compare_df$angle_deg,
    amp_tmm_db,
    col = "#d95f02",
    lwd = 2,
    lty = 2
  )

  plot(
    compare_df$angle_deg,
    compare_df$delta_amp_db,
    type = "l",
    lwd = 2,
    col = "#1b9e77",
    xlab = expression(phi[scatter] ~ "(deg, world x-y slice)"),
    ylab = expression(Delta * " amplitude (dB)"),
    main = ""
  )
  abline(h = 0, lty = 3, col = "grey40")
  if (!is.null(residual_ylim)) {
    graphics::par(usr = c(par("usr")[1:2], residual_ylim))
  }
}

compare_0 <- read_compare_csv(file.path(
  "tools/implementation-figures/data/tmm",
  "bempp_cylinder_pr_70x10_38k_xy_181_0deg_compare.csv"
))
compare_45 <- read_compare_csv(file.path(
  "tools/implementation-figures/data/tmm",
  "bempp_cylinder_pr_70x10_38k_xy_181_45deg_compare.csv"
))
compare_90 <- read_compare_csv(file.path(
  "tools/implementation-figures/data/tmm",
  "bempp_cylinder_pr_70x10_38k_xy_181_compare.csv"
))

png(
  filename = file.path(
    "vignettes",
    "tmm",
    "tmm-bempp-cylinder-validation.png"
  ),
  width = 2400,
  height = 2400,
  res = 200
)
op <- par(no.readonly = TRUE)
on.exit({
  par(op)
  dev.off()
}, add = TRUE)
par(mfrow = c(3, 2), mar = c(4.1, 4.4, 2.6, 1.2), oma = c(0, 0, 1.4, 0))

plot_case(
  compare_df = compare_0,
  main_title = "Cylinder, 70 x 10 mm, 38 kHz, 0 deg incidence"
)
plot_case(
  compare_df = compare_45,
  main_title = "Cylinder, 70 x 10 mm, 38 kHz, 45 deg incidence"
)
plot_case(
  compare_df = compare_90,
  main_title = "Cylinder, 70 x 10 mm, 38 kHz, 90 deg incidence"
)
par(xpd = NA)
mtext(
  "BEMPP vs exploratory retained-angle cylinder TMM",
  outer = TRUE,
  cex = 1.15,
  font = 2
)
legend(
  "top",
  inset = c(0, -0.01),
  horiz = TRUE,
  bty = "n",
  legend = c("BEMPP", "Exploratory retained-angle TMM", expression(Delta * " amplitude")),
  lty = c(1, 2, 1),
  lwd = 2,
  col = c("black", "#d95f02", "#1b9e77")
)
dev.off()

check_monostatic_case <- function(freq,
                                  length_body,
                                  radius_body,
                                  incidence_deg,
                                  compare_csv = NULL,
                                  bempp_txt = NULL) {
  object <- target_strength(
    object = fls_generate(
      shape = cylinder(
        length_body = length_body,
        radius_body = radius_body,
        n_segments = 80
      ),
      g_body = 1,
      h_body = 1,
      theta_body = pi / 2
    ),
    frequency = freq,
    model = "tmm",
    boundary = "pressure_release",
    density_sw = density_sw,
    sound_speed_sw = sound_speed_sw,
    store_t_matrix = TRUE
  )

  incident <- world_to_body_fixed_xaxis(theta = pi / 2, phi = incidence_deg * pi / 180)
  backscatter_world_deg <- (incidence_deg + 180) %% 360
  receive_theta <- pi - incident$theta[[1]]
  receive_phi <- (incident$phi[[1]] + pi) %% (2 * pi)

  tmm_point <- tmm_scattering(
    object,
    theta_body = incident$theta[[1]],
    phi_body = incident$phi[[1]],
    theta_scatter = receive_theta,
    phi_scatter = receive_phi
  )

  if (!is.null(compare_csv)) {
    compare_df <- read_compare_csv(compare_csv)
    angle_vals <- compare_df$angle_deg
    mod_f_bempp <- compare_df$mod_f_bempp
  } else if (!is.null(bempp_txt)) {
    compare_df <- read_bempp_slice(bempp_txt)
    angle_vals <- compare_df$angle_deg
    mod_f_bempp <- compare_df$mod_f_bempp
  } else {
    stop("Either 'compare_csv' or 'bempp_txt' must be supplied.", call. = FALSE)
  }

  idx <- which.min(abs(angle_vals - backscatter_world_deg))

  data.frame(
    length_body_mm = length_body * 1000,
    radius_body_mm = radius_body * 1000,
    frequency_kHz = freq / 1000,
    incidence_deg = incidence_deg,
    backscatter_angle_deg = angle_vals[idx],
    mod_f_bempp = mod_f_bempp[idx],
    mod_f_tmm_monostatic = Mod(tmm_point$f_scat[[1]]),
    delta_amp_db = 20 * log10(
      max(Mod(tmm_point$f_scat[[1]]), .Machine$double.eps) /
        max(mod_f_bempp[idx], .Machine$double.eps)
    )
  )
}

mono_summary <- do.call(
  rbind,
  list(
    check_monostatic_case(
      freq = 38000,
      length_body = 0.07,
      radius_body = 0.01,
      incidence_deg = 0,
      compare_csv = file.path("scratch", "bempp_cylinder_pr_70x10_38k_xy_181_0deg_compare.csv")
    ),
    check_monostatic_case(
      freq = 38000,
      length_body = 0.07,
      radius_body = 0.01,
      incidence_deg = 45,
      compare_csv = file.path("scratch", "bempp_cylinder_pr_70x10_38k_xy_181_45deg_compare.csv")
    ),
    check_monostatic_case(
      freq = 38000,
      length_body = 0.07,
      radius_body = 0.01,
      incidence_deg = 90,
      compare_csv = file.path("scratch", "bempp_cylinder_pr_70x10_38k_xy_181_compare.csv")
    ),
    check_monostatic_case(
      freq = 70000,
      length_body = 0.07,
      radius_body = 0.01,
      incidence_deg = 90,
      bempp_txt = file.path("scratch", "bempp_cylinder_pr_70x10_70k_xy_181.txt")
    ),
    check_monostatic_case(
      freq = 38000,
      length_body = 0.05,
      radius_body = 0.008,
      incidence_deg = 90,
      bempp_txt = file.path("scratch", "bempp_cylinder_pr_50x8_38k_xy_181.txt")
    )
  )
)

write.csv(
  mono_summary,
  file.path(
    "vignettes",
    "data",
    "implementation",
    "tmm",
    "bempp_tmm_cylinder_monostatic_summary.csv"
  ),
  row.names = FALSE
)

print(mono_summary)
