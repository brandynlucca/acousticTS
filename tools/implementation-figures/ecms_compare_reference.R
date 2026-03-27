source("tools/implementation-figures/helpers/common.R")
impl_set_repo_root()
data_dir <- impl_data_path("")

.jc_ref <- function(m, x) {
  besselJ(x, m)
}

.yc_ref <- function(m, x) {
  besselY(x, m)
}

.jcd_ref <- function(m, x) {
  0.5 * (.jc_ref(m - 1, x) - .jc_ref(m + 1, x))
}

.ycd_ref <- function(m, x) {
  0.5 * (.yc_ref(m - 1, x) - .yc_ref(m + 1, x))
}

.safe_divide_ref <- function(numerator, denominator, eps = 1e-12) {
  denominator_adj <- ifelse(
    abs(denominator) < eps,
    ifelse(Re(denominator) < 0, -eps, eps),
    denominator
  )
  numerator / denominator_adj
}

.ecms_reference_scalar <- function(k_sw,
                                   k_l,
                                   k_t,
                                   length_body,
                                   radius_body,
                                   theta_body,
                                   density_sw,
                                   density_body,
                                   m_limit) {
  sin_theta <- abs(sin(theta_body))
  x1 <- k_l * radius_body * sin_theta
  x2 <- k_t * radius_body * sin_theta
  x3 <- k_sw * radius_body * sin_theta

  m <- 0:m_limit
  nu <- ifelse(m == 0, 1, 2)
  m2 <- m^2

  j1 <- .jc_ref(m, x1)
  j2 <- .jc_ref(m, x2)
  j3 <- .jc_ref(m, x3)
  y3 <- .yc_ref(m, x3)

  tan_alpha1 <- -.safe_divide_ref(x1 * .jcd_ref(m, x1), j1)
  tan_alpha2 <- -.safe_divide_ref(x2 * .jcd_ref(m, x2), j2)
  tan_alpha3 <- -.safe_divide_ref(x3 * .jcd_ref(m, x3), j3)
  tan_beta3 <- -.safe_divide_ref(x3 * .ycd_ref(m, x3), y3)
  tan_delta3 <- -.safe_divide_ref(j3, y3)

  x2_half <- x2^2 / 2
  denom_a2 <- m2 - x2_half + tan_alpha2

  term1 <- .safe_divide_ref(tan_alpha1, tan_alpha1 + 1)
  term2 <- .safe_divide_ref(m2, denom_a2)
  term3 <- .safe_divide_ref(m2 - x2_half + tan_alpha1, tan_alpha1 + 1)
  term4 <- .safe_divide_ref(m2 * (tan_alpha2 + 1), denom_a2)

  tan_zeta <- .safe_divide_ref(
    (-x2_half) * (term1 - term2),
    term3 - term4
  )
  tan_phi <- -(density_sw / density_body) * tan_zeta
  tan_eta <- tan_delta3 * .safe_divide_ref(
    tan_phi + tan_alpha3,
    tan_phi + tan_beta3
  )

  cos_eta <- 1 / sqrt(1 + tan_eta^2)
  sin_eta <- tan_eta * cos_eta
  modal_sum <- sum(
    (-1)^m * nu * sin_eta * (cos_eta - 1i * sin_eta),
    na.rm = TRUE
  )

  A <- k_sw * length_body * cos(theta_body)
  sinc_factor <- if (abs(A) < 1e-10) 1 else sin(A) / A

  -(length_body / pi) * sinc_factor * modal_sum
}

frequency <- seq(12e3, 200e3, by = 2e3)
sound_speed_sw <- 1477.3
density_sw <- 1026.8
density_body <- 2800
sound_speed_longitudinal_body <- 6398
sound_speed_transversal_body <- 3122
length_body <- 0.04
radius_body <- 0.005
theta_body <- pi / 2

k_sw <- 2 * pi * frequency / sound_speed_sw
k_l <- 2 * pi * frequency / sound_speed_longitudinal_body
k_t <- 2 * pi * frequency / sound_speed_transversal_body
m_limit <- ceiling(pmax(k_sw, k_l, k_t) * radius_body * abs(sin(theta_body))) + 10

elapsed <- system.time({
  f_bs <- mapply(
    FUN = .ecms_reference_scalar,
    k_sw = k_sw,
    k_l = k_l,
    k_t = k_t,
    m_limit = m_limit,
    MoreArgs = list(
      length_body = length_body,
      radius_body = radius_body,
      theta_body = theta_body,
      density_sw = density_sw,
      density_body = density_body
    )
  )
})["elapsed"]

sigma_bs <- abs(f_bs)^2
ts_db <- 10 * log10(sigma_bs)

write.csv(
  data.frame(
    frequency_hz = frequency,
    frequency_khz = frequency * 1e-3,
    ka = k_sw * radius_body,
    f_bs_real = Re(f_bs),
    f_bs_imag = Im(f_bs),
    sigma_bs = sigma_bs,
    TS = ts_db
  ),
  file.path(data_dir, "ecms_reference_direct.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(
    implementation = "direct-reference",
    elapsed_s = as.numeric(elapsed),
    frequency_min_khz = min(frequency) * 1e-3,
    frequency_max_khz = max(frequency) * 1e-3,
    frequency_step_khz = 2
  ),
  file.path(data_dir, "ecms_reference_direct_timing.csv"),
  row.names = FALSE
)

acoustic_df <- read.csv(file.path(data_dir, "ecms_reference_acousticts.csv"))
reference_df <- read.csv(file.path(data_dir, "ecms_reference_direct.csv"))
timing_acoustic <- read.csv(file.path(data_dir, "ecms_reference_acousticts_timing.csv"))
timing_reference <- read.csv(file.path(data_dir, "ecms_reference_direct_timing.csv"))

merged <- merge(
  acoustic_df,
  reference_df,
  by = c("frequency_hz", "frequency_khz", "ka"),
  suffixes = c("_acousticts", "_reference")
)

# Reorder rows explicitly because base::merge() does not preserve the intended
# frequency sweep ordering for the ECMS comparison table.
merged <- merged[
  order(merged$frequency_hz),
  ,
  drop = FALSE
]

write.csv(
  merged,
  file.path(data_dir, "ecms_reference_compare.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(
    comparison = "acousticTS vs direct reference",
    max_abs_delta_TS_dB = max(abs(merged$TS_acousticts - merged$TS_reference)),
    mean_abs_delta_TS_dB = mean(abs(merged$TS_acousticts - merged$TS_reference)),
    frequency_at_max_delta_kHz = merged$frequency_khz[
      which.max(abs(merged$TS_acousticts - merged$TS_reference))
    ],
    elapsed_acousticts_s = timing_acoustic$elapsed_s,
    elapsed_reference_s = timing_reference$elapsed_s
  ),
  file.path(data_dir, "ecms_reference_compare_summary.csv"),
  row.names = FALSE
)
