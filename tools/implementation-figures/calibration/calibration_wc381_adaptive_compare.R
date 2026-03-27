source("tools/implementation-figures/helpers/common.R")
impl_load_all()
data_dir <- impl_data_path("")
ref <- read.csv(file.path(data_dir, 'calibration_wc381_package_compare.csv'))
obj <- cal_generate(material = 'WC', diameter = 38.1e-3)
model <- calibration_initialize(obj, frequency = ref$frequency, sound_speed_sw = 1477.3, density_sw = 1026.8)@model_parameters$calibration

calc_adaptive <- function(model, tol = 1e-10) {
  ka_sw_v <- model$parameters$acoustics$k_sw * model$body$radius
  ka_l_v <- model$parameters$acoustics$k_l * model$body$radius
  ka_t_v <- model$parameters$acoustics$k_t * model$body$radius
  g <- model$body$density / model$medium$density
  theta <- model$body$theta
  out <- mapply(function(ka_sw, ka_l, ka_t) {
    ml <- as.integer(round(ka_sw) + 10L)
    repeat {
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
      terms <- (2 * m + 1) * Pl * (sin_eta * (1i * cos_eta - sin_eta))
      if (Mod(tail(terms, 1)) <= tol || ml >= 500L) break
      ml <- ml + 1L
    }
    f_j <- sum(terms)
    f_bs <- abs(-2i * f_j / ka_sw) * model$body$radius / 2
    10 * log10(f_bs * f_bs)
  }, ka_sw_v, ka_l_v, ka_t_v)
  as.numeric(out)
}

ts <- calc_adaptive(model)
d1 <- abs(ts - ref$echoSMs)
d2 <- abs(ts - ref$sphereTS)
cat('adaptive vs echoSMs mean/max\n')
print(c(mean=mean(d1), max=max(d1)))
cat('adaptive vs sphereTS mean/max\n')
print(c(mean=mean(d2), max=max(d2)))
