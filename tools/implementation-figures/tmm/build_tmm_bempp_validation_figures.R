source("tools/implementation-figures/helpers/common.R")
impl_load_all()
source("tools/implementation-figures/helpers/compare_bempp_tmm_slice.R")

density_sw <- 1026.8
sound_speed_sw <- 1480

build_compare <- function(object,
                          bempp_path,
                          out_csv,
                          out_summary,
                          frame = "body_fixed",
                          incidence_deg = NULL,
                          theta_body = pi / 2,
                          phi_body = pi / 2,
                          theta_scatter = pi / 2) {
  cmp <- compare_bempp_tmm_slice(
    object = object,
    bempp_path = bempp_path,
    frame = frame,
    incidence_deg = incidence_deg,
    theta_body = theta_body,
    phi_body = phi_body,
    theta_scatter = theta_scatter
  )

  write.csv(cmp$detail, out_csv, row.names = FALSE)
  write.csv(cmp$summary, out_summary, row.names = FALSE)
  cmp$detail
}

plot_validation_figure <- function(detail,
                                   title_top,
                                   title_bottom,
                                   out_path) {
  detail <- detail[order(detail$angle_deg), , drop = FALSE]
  bempp_db <- 20 * log10(pmax(detail$mod_f_bempp, .Machine$double.eps))
  tmm_db <- 20 * log10(pmax(detail$mod_f_tmm, .Machine$double.eps))
  delta_db <- detail$delta_amp_db

  png(out_path, width = 1600, height = 1400, res = 180)
  on.exit(dev.off(), add = TRUE)

  op <- par(no.readonly = TRUE)
  on.exit(par(op), add = TRUE)

  layout(matrix(c(1, 2), ncol = 1), heights = c(3.2, 1.2))
  par(mar = c(3.5, 4.5, 2.8, 1.2))

  plot(
    detail$angle_deg,
    bempp_db,
    type = "l",
    lwd = 2.5,
    col = "black",
    xlab = "",
    ylab = expression(10 * log[10](sigma[scat] / (1 ~ m^2)) ~ "(dB)"),
    main = title_top
  )
  lines(detail$angle_deg, tmm_db, col = "#c75d2c", lwd = 2, lty = 2)
  points(detail$angle_deg[seq(1, nrow(detail), by = 6)],
         tmm_db[seq(1, nrow(detail), by = 6)],
         pch = 21, bg = "white", col = "#c75d2c", cex = 0.8)
  legend(
    "topright",
    legend = c("BEMPP", "TMM"),
    col = c("black", "#c75d2c"),
    lty = c(1, 2),
    lwd = c(2.5, 2),
    pch = c(NA, 21),
    pt.bg = c(NA, "white"),
    bty = "n"
  )

  par(mar = c(4.5, 4.5, 1.2, 1.2))
  plot(
    detail$angle_deg,
    delta_db,
    type = "l",
    lwd = 2,
    col = "#365f91",
    xlab = "World-plane receive angle (deg)",
    ylab = expression(Delta ~ "amplitude" ~ "(dB)"),
    main = title_bottom
  )
  abline(h = 0, col = "grey50", lty = 3)
}

sphere_object <- target_strength(
  object = fls_generate(
    shape = sphere(radius_body = 0.01, n_segments = 80),
    g_body = 1,
    h_body = 1,
    theta_body = pi / 2
  ),
  frequency = 38000,
  model = "tmm",
  boundary = "pressure_release",
  density_sw = density_sw,
  sound_speed_sw = sound_speed_sw,
  store_t_matrix = TRUE
)

oblate_object <- target_strength(
  object = fls_generate(
    shape = oblate_spheroid(length_body = 0.012, radius_body = 0.01, n_segments = 80),
    g_body = 1,
    h_body = 1,
    theta_body = pi / 2
  ),
  frequency = 38000,
  model = "tmm",
  boundary = "pressure_release",
  density_sw = density_sw,
  sound_speed_sw = sound_speed_sw,
  store_t_matrix = TRUE
)

prolate_object <- target_strength(
  object = fls_generate(
    shape = prolate_spheroid(length_body = 0.07, radius_body = 0.01, n_segments = 80),
    g_body = 1,
    h_body = 1,
    theta_body = pi / 2
  ),
  frequency = 38000,
  model = "tmm",
  boundary = "pressure_release",
  density_sw = density_sw,
  sound_speed_sw = sound_speed_sw,
  store_t_matrix = TRUE
)

sphere_detail <- build_compare(
  object = sphere_object,
  bempp_path = paste(
    "tools/implementation-figures/data/tmm",
    "bempp_sphere_pr_r10_38k_xy_181.txt",
    sep = "/"
  ),
  frame = "body_fixed",
  theta_body = pi / 2,
  phi_body = pi / 2,
  theta_scatter = pi / 2,
  out_csv = paste(
    "tools/implementation-figures/data/tmm",
    "bempp_sphere_pr_r10_38k_xy_181_tmm_compare.csv",
    sep = "/"
  ),
  out_summary = paste(
    "tools/implementation-figures/data/tmm",
    "bempp_sphere_pr_r10_38k_xy_181_tmm_summary.csv",
    sep = "/"
  )
)

oblate_detail <- build_compare(
  object = oblate_object,
  bempp_path = paste(
    "tools/implementation-figures/data/tmm",
    "bempp_oblate_pr_c6_a10_38k_xy_181.txt",
    sep = "/"
  ),
  frame = "world_xy_xaxis",
  incidence_deg = 90,
  out_csv = paste(
    "tools/implementation-figures/data/tmm",
    "bempp_oblate_pr_c6_a10_38k_xy_181_compare.csv",
    sep = "/"
  ),
  out_summary = paste(
    "tools/implementation-figures/data/tmm",
    "bempp_oblate_pr_c6_a10_38k_xy_181_summary.csv",
    sep = "/"
  )
)

prolate_detail <- build_compare(
  object = prolate_object,
  bempp_path = paste(
    "tools/implementation-figures/data/tmm",
    "bempp_prolate_pr_70x10_38k_xy_181.txt",
    sep = "/"
  ),
  frame = "world_xy_xaxis",
  incidence_deg = 90,
  out_csv = paste(
    "tools/implementation-figures/data/tmm",
    "bempp_prolate_pr_70x10_38k_xy_181_compare.csv",
    sep = "/"
  ),
  out_summary = paste(
    "tools/implementation-figures/data/tmm",
    "bempp_prolate_pr_70x10_38k_xy_181_summary.csv",
    sep = "/"
  )
)

plot_validation_figure(
  detail = sphere_detail,
  title_top = "Pressure-release sphere: TMM vs BEMPP at 38 kHz",
  title_bottom = "Amplitude residual (TMM - BEMPP)",
  out_path = "vignettes/tmm/tmm-bempp-sphere-validation.png"
)

plot_validation_figure(
  detail = oblate_detail,
  title_top = "Pressure-release oblate: TMM vs BEMPP at 38 kHz",
  title_bottom = "Amplitude residual (TMM - BEMPP)",
  out_path = "vignettes/tmm/tmm-bempp-oblate-validation.png"
)

plot_validation_figure(
  detail = prolate_detail,
  title_top = "Pressure-release prolate: TMM vs BEMPP at 38 kHz",
  title_bottom = "Amplitude residual (TMM - BEMPP)",
  out_path = "vignettes/tmm/tmm-bempp-prolate-validation.png"
)
