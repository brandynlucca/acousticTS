source("tools/implementation-figures/helpers/common.R")
impl_load_all()
data_dir <- impl_data_path("")

.equivalent_length_fresnel_ref <- function(k1, l, a, rho_c) {
  gamma_max <- l / (2 * rho_c)
  z_max <- rho_c * (1 - cos(gamma_max))
  A <- z_max / a
  B <- l / a

  vapply(
    k1,
    function(k_scalar) {
      integrand <- function(x) {
        exp(1i * 8 * k_scalar * A * x^2 / (B^2 * a))
      }
      complex(
        real = integrate(function(x) Re(integrand(x)), -l / 2, l / 2)$value,
        imaginary = integrate(function(x) Im(integrand(x)), -l / 2, l / 2)$value
      )
    },
    complex(1)
  )
}

frequency <- seq(12e3, 400e3, by = 2e3)
density_sw <- 1026.8
sound_speed_sw <- 1477.3
density_body <- density_sw * 1.0357
sound_speed_body <- sound_speed_sw * 1.0279
length_body <- 10.5e-3
radius_body <- 1e-3
radius_curvature_ratio <- 1.5
radius_curvature <- radius_curvature_ratio * length_body

straight_object <- fls_generate(
  shape = cylinder(
    length_body = length_body,
    radius_body = radius_body,
    n_segments = 401
  ),
  density_body = density_body,
  sound_speed_body = sound_speed_body,
  theta_body = pi / 2
)

timing_path <- file.path(data_dir, "bcms_reference_direct_timing.csv")

if (impl_should_refresh_timings() || !file.exists(timing_path)) {
  elapsed_straight <- system.time({
    straight_ref <- target_strength(
      straight_object,
      frequency = frequency,
      model = "fcms",
      boundary = "liquid_filled",
      density_sw = density_sw,
      sound_speed_sw = sound_speed_sw
    )
  })["elapsed"]

  elapsed_bent <- system.time({
    lebc <- .equivalent_length_fresnel_ref(
      k1 = 2 * pi * frequency / sound_speed_sw,
      l = length_body,
      a = radius_body,
      rho_c = radius_curvature
    )
  })["elapsed"]

  timing_reference <- data.frame(
    case = c("straight", "bent"),
    implementation = c("FCMS-reference", "Stanton-1989-reference"),
    elapsed_s = c(as.numeric(elapsed_straight), as.numeric(elapsed_bent)),
    frequency_min_khz = min(frequency) * 1e-3,
    frequency_max_khz = max(frequency) * 1e-3,
    frequency_step_khz = 2
  )
  timing_reference <- impl_round_timing_columns(timing_reference)

  write.csv(
    timing_reference,
    timing_path,
    row.names = FALSE
  )
} else {
  timing_reference <- read.csv(timing_path)
}

straight_ref <- target_strength(
  straight_object,
  frequency = frequency,
  model = "fcms",
  boundary = "liquid_filled",
  density_sw = density_sw,
  sound_speed_sw = sound_speed_sw
)

straight_f_bs <- extract(straight_ref, "model")$FCMS$f_bs
straight_ts <- extract(straight_ref, "model")$FCMS$TS

lebc <- .equivalent_length_fresnel_ref(
  k1 = 2 * pi * frequency / sound_speed_sw,
  l = length_body,
  a = radius_body,
  rho_c = radius_curvature
)

bent_f_bs <- lebc * straight_f_bs / length_body
bent_sigma_bs <- abs(bent_f_bs)^2
bent_ts <- 10 * log10(bent_sigma_bs)

write.csv(
  rbind(
    data.frame(
      case = "straight",
      frequency_hz = frequency,
      frequency_khz = frequency * 1e-3,
      ka = 2 * pi * frequency / sound_speed_sw * radius_body,
      f_bs_real = Re(straight_f_bs),
      f_bs_imag = Im(straight_f_bs),
      sigma_bs = abs(straight_f_bs)^2,
      TS = straight_ts
    ),
    data.frame(
      case = "bent",
      frequency_hz = frequency,
      frequency_khz = frequency * 1e-3,
      ka = 2 * pi * frequency / sound_speed_sw * radius_body,
      f_bs_real = Re(bent_f_bs),
      f_bs_imag = Im(bent_f_bs),
      sigma_bs = bent_sigma_bs,
      TS = bent_ts
    )
  ),
  file.path(data_dir, "bcms_reference_direct.csv"),
  row.names = FALSE
)

acoustic_df <- read.csv(file.path(data_dir, "bcms_reference_acousticts.csv"))
reference_df <- read.csv(file.path(data_dir, "bcms_reference_direct.csv"))

merged <- merge(
  acoustic_df,
  reference_df,
  by = c("case", "frequency_hz", "frequency_khz", "ka"),
  suffixes = c("_acousticts", "_reference")
)

# Reorder rows explicitly because base::merge() does not preserve the intended
# frequency sweep ordering for the bent/straight comparison table.
case_levels <- c("straight", "bent")
merged <- merged[
  order(match(merged$case, case_levels), merged$frequency_hz),
  ,
  drop = FALSE
]

write.csv(
  merged,
  file.path(data_dir, "bcms_reference_compare.csv"),
  row.names = FALSE
)

summary_df <- do.call(
  rbind,
  lapply(split(merged, merged$case), function(df) {
    data.frame(
      case = unique(df$case),
      max_abs_delta_TS_dB = max(abs(df$TS_acousticts - df$TS_reference)),
      mean_abs_delta_TS_dB = mean(abs(df$TS_acousticts - df$TS_reference)),
      frequency_at_max_delta_kHz = df$frequency_khz[
        which.max(abs(df$TS_acousticts - df$TS_reference))
      ]
    )
  })
)

timing_df <- merge(
  read.csv(file.path(data_dir, "bcms_reference_acousticts_timing.csv")),
  timing_reference,
  by = "case",
  suffixes = c("_acousticts", "_reference")
)

summary_df <- merge(summary_df, timing_df, by = "case")
summary_df <- impl_round_timing_columns(summary_df)

write.csv(
  summary_df,
  file.path(data_dir, "bcms_reference_compare_summary.csv"),
  row.names = FALSE
)
