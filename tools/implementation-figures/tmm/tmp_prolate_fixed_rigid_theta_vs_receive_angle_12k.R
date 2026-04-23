library(devtools)
load_all("C:/Users/Brandyn/Desktop/acousticTS", quiet = TRUE)

density_sw <- 1026.8
sound_speed_sw <- 1477.3
theta_body <- pi / 2
phi_body <- pi / 2

object <- target_strength(
  object = fls_generate(
    shape = prolate_spheroid(length_body = 0.14, radius_body = 0.01, n_segments = 80),
    density_body = density_sw,
    sound_speed_body = sound_speed_sw,
    theta_body = theta_body
  ),
  frequency = 12e3,
  model = "tmm",
  boundary = "fixed_rigid",
  density_sw = density_sw,
  sound_speed_sw = sound_speed_sw,
  store_t_matrix = TRUE
)

grid <- acousticTS:::.tmm_incident_local_polar_grid(
  object = object,
  frequency = 12e3,
  theta_body = theta_body,
  phi_body = phi_body,
  n_theta = 121,
  n_phi = 241
)

z <- grid$sigma_scat_dB
alpha_deg <- grid$alpha_scatter * 180 / pi
theta_world_deg <- grid$theta_scatter_world * 180 / pi
psi_deg <- rep(grid$psi_scatter * 180 / pi, times = length(alpha_deg))
plot_df <- data.frame(
  alpha_deg = rep(alpha_deg, each = nrow(theta_world_deg)),
  theta_scatter_deg = as.vector(theta_world_deg),
  psi_scatter_deg = psi_deg,
  TS_tmm = as.vector(z)
)

csv_path <- "tools/implementation-figures/data/final/tmm/tmm_prolate_fixed_rigid_theta_vs_receive_angle_12kHz.csv"
png_path <- "tools/implementation-figures/data/final/tmm/tmm_prolate_fixed_rigid_theta_vs_receive_angle_12kHz.png"
dir.create(dirname(csv_path), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(plot_df, csv_path, row.names = FALSE)

z_finite <- plot_df$TS_tmm[is.finite(plot_df$TS_tmm)]
z_breaks <- seq(min(z_finite), max(z_finite), length.out = 129L)
palette_vals <- grDevices::hcl.colors(128, "Spectral", rev = TRUE)
col_idx <- findInterval(plot_df$TS_tmm, z_breaks, all.inside = TRUE)
pt_cols <- palette_vals[col_idx]

grDevices::png(filename = png_path, width = 1800, height = 1300, res = 220)
op <- graphics::par(no.readonly = TRUE)
on.exit({
  graphics::par(op)
  grDevices::dev.off()
}, add = TRUE)

graphics::layout(matrix(c(1, 2), nrow = 1), widths = c(1, 0.075))
graphics::par(mar = c(5.2, 5.6, 3.0, 2.0))
graphics::plot(
  NA,
  NA,
  xlim = c(0, 360),
  ylim = c(0, 180),
  xaxs = "i",
  yaxs = "i",
  axes = FALSE,
  xlab = "Receive Angle alpha (deg)",
  ylab = expression(theta[scatter] ~ "(deg)"),
  main = "Prolate fixed-rigid TMM: theta vs receive angle at 12.0 kHz"
)
graphics::axis(1, at = seq(0, 360, by = 45))
graphics::axis(2, at = seq(0, 180, by = 30), las = 1)
graphics::grid(col = "grey85", lty = 3)
graphics::points(
  x = plot_df$alpha_deg,
  y = plot_df$theta_scatter_deg,
  pch = 15,
  cex = 0.36,
  col = pt_cols
)
graphics::box()

graphics::par(mar = c(5.2, 0.8, 3.0, 3.4))
graphics::plot.new()
graphics::plot.window(xlim = c(0, 1), ylim = range(z_breaks))
for (i in seq_along(palette_vals)) {
  graphics::rect(
    xleft = 0,
    ybottom = z_breaks[i],
    xright = 1,
    ytop = z_breaks[i + 1L],
    col = palette_vals[i],
    border = NA
  )
}
graphics::axis(4, las = 1)
graphics::mtext("TS (dB)", side = 4, line = 2.1, cex = 0.9)

cat(normalizePath(csv_path, winslash = "/", mustWork = TRUE), "\n")
cat(normalizePath(png_path, winslash = "/", mustWork = TRUE), "\n")
