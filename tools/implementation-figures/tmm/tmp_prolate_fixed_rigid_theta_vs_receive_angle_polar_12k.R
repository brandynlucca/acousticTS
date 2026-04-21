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

grid <- suppressWarnings(tmm_scattering_grid(
  object = object,
  frequency = 12e3,
  theta_body = theta_body,
  phi_body = phi_body,
  n_theta = 121,
  n_phi = 241
))

basis <- acousticTS:::.tmm_forward_basis(theta_body = theta_body, phi_body = phi_body)
theta_edges <- acousticTS:::.tmm_grid_edges(grid$theta_scatter, lower = 0, upper = pi)
phi_edges <- acousticTS:::.tmm_grid_edges(grid$phi_scatter, lower = 0, upper = 2 * pi)

corner_receive_angle <- function(theta_val, phi_val) {
  scat_vec <- acousticTS:::.tmm_spherical_to_cartesian(theta_val, phi_val)
  alpha_val <- atan2(sum(scat_vec * basis$e2), sum(scat_vec * basis$e1))
  if (alpha_val < 0) {
    alpha_val <- alpha_val + 2 * pi
  }
  alpha_val
}

alpha_center <- outer(
  grid$theta_scatter,
  grid$phi_scatter,
  Vectorize(corner_receive_angle)
)

plot_df <- expand.grid(
  theta_scatter_rad = grid$theta_scatter,
  phi_scatter_rad = grid$phi_scatter
)
plot_df$receive_angle_rad <- as.vector(alpha_center)
plot_df$receive_angle_deg <- plot_df$receive_angle_rad * 180 / pi
plot_df$TS_tmm <- as.vector(grid$sigma_scat_dB)

csv_path <- "tools/implementation-figures/data/final/tmm/tmm_prolate_fixed_rigid_theta_vs_receive_angle_polar_12kHz.csv"
png_path <- "tools/implementation-figures/data/final/tmm/tmm_prolate_fixed_rigid_theta_vs_receive_angle_polar_12kHz.png"
dir.create(dirname(csv_path), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(plot_df, csv_path, row.names = FALSE)

palette_vals <- grDevices::hcl.colors(128, "Spectral", rev = TRUE)
z_vals <- grid$sigma_scat_dB
z_finite <- z_vals[is.finite(z_vals)]
z_breaks <- seq(min(z_finite), max(z_finite), length.out = length(palette_vals) + 1L)

grDevices::png(filename = png_path, width = 1800, height = 1400, res = 220)
op <- graphics::par(no.readonly = TRUE)
on.exit({
  graphics::par(op)
  grDevices::dev.off()
}, add = TRUE)

graphics::par(fig = c(0.0, 0.83, 0.0, 1.0), mar = c(2.2, 2.2, 3.0, 1.4))
graphics::plot(
  NA,
  NA,
  xlim = c(-pi, pi),
  ylim = c(-pi, pi),
  asp = 1,
  axes = FALSE,
  xlab = "",
  ylab = "",
  xaxs = "i",
  yaxs = "i",
  main = "Prolate fixed-rigid TMM: theta vs receive angle at 12.0 kHz"
)

for (j in seq_along(grid$phi_scatter)) {
  for (i in seq_along(grid$theta_scatter)) {
    z_ij <- z_vals[i, j]
    if (!is.finite(z_ij)) {
      next
    }

    theta_poly <- c(
      theta_edges[i],
      theta_edges[i + 1L],
      theta_edges[i + 1L],
      theta_edges[i]
    )
    alpha_poly <- c(
      corner_receive_angle(theta_edges[i], phi_edges[j]),
      corner_receive_angle(theta_edges[i + 1L], phi_edges[j]),
      corner_receive_angle(theta_edges[i + 1L], phi_edges[j + 1L]),
      corner_receive_angle(theta_edges[i], phi_edges[j + 1L])
    )

    x_poly <- theta_poly * cos(alpha_poly)
    y_poly <- theta_poly * sin(alpha_poly)
    col_idx <- findInterval(z_ij, z_breaks, all.inside = TRUE)

    graphics::polygon(
      x_poly,
      y_poly,
      col = palette_vals[col_idx],
      border = NA
    )
  }
}

for (r in seq(pi / 4, pi, by = pi / 4)) {
  graphics::symbols(0, 0, circles = r, inches = FALSE, add = TRUE, fg = "grey75", bg = NA)
}
for (aa in seq(0, 330, by = 30) * pi / 180) {
  graphics::segments(0, 0, pi * cos(aa), pi * sin(aa), col = "grey80", lty = 3)
}
graphics::symbols(0, 0, circles = pi, inches = FALSE, add = TRUE, fg = "black", bg = NA)

graphics::text(
  x = c(pi + 0.14, 0, -pi - 0.14, 0),
  y = c(0, pi + 0.14, 0, -pi - 0.14),
  labels = c("0", "90", "180", "270"),
  xpd = NA,
  cex = 0.95,
  font = 2
)
graphics::text(
  x = c(pi / 4, pi / 2, 3 * pi / 4, pi),
  y = 0,
  labels = c(expression(pi / 4), expression(pi / 2), expression(3 * pi / 4), expression(pi)),
  pos = 3,
  col = "grey35",
  cex = 0.82
)
graphics::mtext("Receive Angle (deg)", side = 1, line = 0.6, cex = 1.0)
graphics::mtext(expression(theta[scatter] ~ "(rad)"), side = 2, line = 0.8, cex = 1.0)

graphics::par(fig = c(0.82, 0.98, 0.08, 0.92), mar = c(1.0, 0.2, 1.0, 2.6), new = TRUE)
graphics::plot.new()
graphics::plot.window(xlim = c(0, 1), ylim = range(z_breaks))
y_edges <- seq(min(z_breaks), max(z_breaks), length.out = length(palette_vals) + 1L)
for (k in seq_along(palette_vals)) {
  graphics::rect(
    xleft = 0,
    ybottom = y_edges[k],
    xright = 1,
    ytop = y_edges[k + 1L],
    col = palette_vals[k],
    border = NA
  )
}
graphics::axis(4, las = 1)
graphics::mtext(expression(TS ~ "(dB)"), side = 4, line = 2.0, cex = 0.9)

cat(normalizePath(csv_path, winslash = "/", mustWork = TRUE), "\n")
cat(normalizePath(png_path, winslash = "/", mustWork = TRUE), "\n")
