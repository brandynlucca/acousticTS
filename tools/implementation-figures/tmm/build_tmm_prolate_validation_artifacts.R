source("tools/implementation-figures/helpers/common.R")
impl_load_all()

dir.create("vignettes/tmm", showWarnings = FALSE, recursive = TRUE)
dir.create(
  "tools/implementation-figures/data/tmm",
  showWarnings = FALSE,
  recursive = TRUE
)

density_sw <- 1026.8
sound_speed_sw <- 1480

# ---------------------------------------------------------------------------
# 1. Equal-volume sphere-to-spheroid continuation path
# ---------------------------------------------------------------------------
continuation_object <- target_strength(
  object = fls_generate(
    shape = prolate_spheroid(length_body = 0.07, radius_body = 0.01, n_segments = 80),
    g_body = 1,
    h_body = 1,
    theta_body = pi / 2
  ),
  frequency = c(38e3, 70e3),
  model = "tmm",
  boundary = "pressure_release",
  density_sw = density_sw,
  sound_speed_sw = sound_speed_sw,
  store_t_matrix = TRUE
)
continuation_diag <- tmm_diagnostics(
  continuation_object,
  continuation_steps = 8L,
  n_theta = 21,
  n_phi = 41
)
continuation_df <- continuation_diag$continuation
write.csv(
  continuation_df,
  "tools/implementation-figures/data/tmm/tmm_sphere_to_spheroid_continuation.csv",
  row.names = FALSE
)

png(
  "vignettes/tmm/tmm-sphere-to-spheroid-continuation.png",
  width = 1800,
  height = 1200,
  res = 200
)
par(mar = c(4.4, 4.6, 2.4, 1.2))
cols <- c("#1b9e77", "#d95f02")
plot(
  NA,
  xlim = range(continuation_df$aspect_ratio),
  ylim = range(continuation_df$TS),
  xlab = "Aspect ratio",
  ylab = "Target strength (dB)",
  main = "Sphere-to-spheroid continuation: pressure-release prolate"
)
for (i in seq_along(sort(unique(continuation_df$frequency)))) {
  freq_i <- sort(unique(continuation_df$frequency))[i]
  df_i <- continuation_df[continuation_df$frequency == freq_i, ]
  lines(df_i$aspect_ratio, df_i$TS, type = "b", lwd = 2, pch = 16, col = cols[i])
}
abline(v = 1, lty = 3, col = "grey45")
legend(
  "bottomright",
  legend = paste0(sort(unique(continuation_df$frequency)) / 1e3, " kHz"),
  col = cols,
  lwd = 2,
  pch = 16,
  bty = "n"
)
dev.off()

# ---------------------------------------------------------------------------
# The earlier paper-style prolate plots were removed because they did not
# reproduce the source papers on a true apples-to-apples basis. The current
# retained-prolate validation figure is built separately in:
#   scratch/build_tmm_prolate_exact_validation_figure.R
# and compares the stored TMM angular field directly against the exact
# general-angle spheroidal solution at matched geometry, angles, and
# normalization.
# ---------------------------------------------------------------------------
