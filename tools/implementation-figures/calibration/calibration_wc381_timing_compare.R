source("tools/implementation-figures/helpers/common.R")
impl_load_all()
data_dir <- impl_data_path("")
output_path <- file.path(data_dir, "calibration_wc381_timing_compare.csv")

calc_old_ts <- function(freq) {
  obj <- cal_generate(material = 'WC', diameter = 38.1e-3)
  model <- calibration_initialize(obj, frequency = freq, sound_speed_sw = 1477.3, density_sw = 1026.8)@model_parameters$calibration
  ka_sw <- model$parameters$acoustics$k_sw * model$body$radius
  ka_l <- model$parameters$acoustics$k_l * model$body$radius
  ka_t <- model$parameters$acoustics$k_t * model$body$radius
  m_limit <- round(ka_sw) + 10
  f_j <- mapply(FUN = function(ka_sw, ka_l, ka_t, theta, density_body, density_sw, ml) {
    m <- 0:ml
    Pl <- Pn(m, cos(theta))[,1]
    js_mat <- js(m, ka_sw)
    js_mat_l <- js(m, ka_l)
    js_mat_t <- js(m, ka_t)
    ys_mat <- ys(m, ka_sw)
    jsd_mat <- jsd(m, ka_sw)
    jsd_mat_l <- jsd(m, ka_l)
    jsd_mat_t <- jsd(m, ka_t)
    ysd_mat <- ysd(m, ka_sw)
    g <- density_body / density_sw
    tan_sw <- -ka_sw * jsd_mat / js_mat
    tan_l <- -ka_l * jsd_mat_l / js_mat_l
    tan_t <- -ka_t * jsd_mat_t / js_mat_t
    tan_beta <- -ka_sw * ysd_mat / ys_mat
    tan_diff <- -js_mat / ys_mat
    along_m <- (m * m + m)
    tan_l_add <- tan_l + 1
    tan_t_div <- along_m - 1 - ka_t * ka_t / 2 + tan_t
    numerator <- (tan_l / tan_l_add) - (along_m / tan_t_div)
    denominator1 <- (along_m - ka_t * ka_t / 2 + 2 * tan_l) / tan_l_add
    denominator2 <- along_m * (tan_t + 1) / tan_t_div
    denominator <- denominator1 - denominator2
    ratio <- -0.5 * (ka_t * ka_t) * numerator / denominator
    phi <- -ratio / g
    eta_tan <- tan_diff * (phi + tan_sw) / (phi + tan_beta)
    cos_eta <- 1 / sqrt(1 + eta_tan * eta_tan)
    sin_eta <- eta_tan * cos_eta
    sum((2 * m + 1) * Pl * (sin_eta * (1i * cos_eta - sin_eta)))
  }, ka_sw, ka_l, ka_t, model$body$theta, model$body$density, model$medium$density, m_limit)
  f_bs <- abs(-2i * f_j / ka_sw) * model$body$radius / 2
  10 * log10(f_bs * f_bs)
}

calc_new_ts <- function(freq) {
  obj <- cal_generate(material = 'WC', diameter = 38.1e-3)
  obj <- target_strength(obj, frequency = freq, model = 'calibration', sound_speed_sw = 1477.3, density_sw = 1026.8)
  obj@model$calibration$TS
}

freq <- seq(1e3, 360e3, 1e3)

if (impl_should_refresh_timings() || !file.exists(output_path)) {
  old_times <- numeric(5)
  new_times <- numeric(5)
  for (i in seq_len(5)) {
    old_times[i] <- system.time(calc_old_ts(freq))[['elapsed']]
    new_times[i] <- system.time(calc_new_ts(freq))[['elapsed']]
  }
  res <- data.frame(
    version = c('old_fixed_cutoff', 'new_adaptive_cutoff'),
    mean_elapsed_s = c(mean(old_times), mean(new_times)),
    median_elapsed_s = c(median(old_times), median(new_times)),
    min_elapsed_s = c(min(old_times), min(new_times)),
    max_elapsed_s = c(max(old_times), max(new_times))
  )
  res <- impl_round_timing_columns(res)

  write.csv(
    res,
    output_path,
    row.names = FALSE
  )
} else {
  res <- read.csv(output_path)
}
print(res)
