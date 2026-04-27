library(devtools)
load_all("C:/Users/Brandyn/Desktop/acousticTS", quiet = TRUE)

density_sw <- 1026.8
sound_speed_sw <- 1477.3
theta_body <- pi / 2
phi_body <- pi / 2
frequency <- 12e3
n_theta <- 121
n_phi <- 241
n_psi <- 241

out_dir <- "tools/implementation-figures/data/final/tmm"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

world_png <- file.path(out_dir, "tmm_prolate_fixed_rigid_world_heatmap_12kHz.png")
local_png <- file.path(out_dir, "tmm_prolate_fixed_rigid_incident_local_heatmap_12kHz.png")
cut_png <- file.path(out_dir, "tmm_prolate_fixed_rigid_incident_plane_cut_12kHz.png")
cut_csv <- file.path(out_dir, "tmm_prolate_fixed_rigid_incident_plane_cut_12kHz.csv")

object <- target_strength(
  object = fls_generate(
    shape = prolate_spheroid(length_body = 0.14, radius_body = 0.01, n_segments = 80),
    density_body = density_sw,
    sound_speed_body = sound_speed_sw,
    theta_body = theta_body
  ),
  frequency = frequency,
  model = "tmm",
  boundary = "fixed_rigid",
  density_sw = density_sw,
  sound_speed_sw = sound_speed_sw,
  store_t_matrix = TRUE
)

grDevices::png(world_png, width = 1800, height = 1500, res = 220)
plot(
  object,
  type = "scattering",
  frequency = frequency,
  heatmap = TRUE,
  frame = "world",
  n_theta = n_theta,
  n_phi = n_phi
)
grDevices::dev.off()

grDevices::png(local_png, width = 1800, height = 1500, res = 220)
plot(
  object,
  type = "scattering",
  frequency = frequency,
  heatmap = TRUE,
  frame = "incident_local",
  n_theta = n_theta,
  n_phi = n_phi
)
grDevices::dev.off()

summary_obj <- tmm_bistatic_summary(
  object = object,
  frequency = frequency,
  theta_body = theta_body,
  phi_body = phi_body,
  n_theta = n_theta,
  n_phi = n_phi,
  n_psi = n_psi,
  include_grid = FALSE
)

cut_df <- summary_obj$slices$forward_scatter
utils::write.csv(cut_df, cut_csv, row.names = FALSE)

grDevices::png(cut_png, width = 1800, height = 1200, res = 220)
op <- graphics::par(no.readonly = TRUE)
on.exit({
  graphics::par(op)
  if (names(grDevices::dev.cur()) != "null device") {
    grDevices::dev.off()
  }
}, add = TRUE)
graphics::par(mar = c(5.0, 5.6, 2.5, 1.8))
graphics::plot(
  x = cut_df$psi_scatter,
  y = cut_df$sigma_scat_dB,
  type = "l",
  lwd = 3,
  col = "black",
  xlab = expression(psi[scatter] ~ "(rad)"),
  ylab = expression(10 * log[10](sigma[scat] / (1 ~ m^2)) ~ "(dB)"),
  main = "Prolate fixed-rigid TMM incident-plane cut at 12.0 kHz",
  xaxt = "n"
)
graphics::grid(col = "grey85", lty = 3)
graphics::axis(
  1,
  at = c(0, pi / 2, pi),
  labels = c(expression(0), expression(pi / 2), expression(pi))
)
grDevices::dev.off()

cat(normalizePath(world_png, winslash = "/", mustWork = TRUE), "\n")
cat(normalizePath(local_png, winslash = "/", mustWork = TRUE), "\n")
cat(normalizePath(cut_png, winslash = "/", mustWork = TRUE), "\n")
cat(normalizePath(cut_csv, winslash = "/", mustWork = TRUE), "\n")
