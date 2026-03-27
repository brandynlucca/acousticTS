source("tools/implementation-figures/helpers/common.R")
impl_load_all()
data_dir <- impl_data_path("")

calc_variant <- function(freq, m_extra = 10L, adaptive = FALSE, tol = 1e-10) {
  obj <- cal_generate(material = 'WC', diameter = 38.1e-3)
  model <- calibration_initialize(obj, frequency = freq, sound_speed_sw = 1477.3, density_sw = 1026.8)@model_parameters$calibration
  ka_sw <- model$parameters$acoustics$k_sw * model$body$radius
  ka_l <- model$parameters$acoustics$k_l * model$body$radius
  ka_t <- model$parameters$acoustics$k_t * model$body$radius
  g <- model$body$density / model$medium$density
  theta <- model$body$theta
  one_freq <- function(ka_sw, ka_l, ka_t) {
    sum_to_ml <- function(ml) {
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
      list(sum = sum(terms), last = tail(terms, 1), terms = terms)
    }
    ml <- as.integer(round(ka_sw) + m_extra)
    res <- sum_to_ml(ml)
    if (adaptive) {
      while (Mod(res$last) > tol && ml < 500L) {
        ml <- ml + 1L
        res <- sum_to_ml(ml)
      }
    }
    f_j <- res$sum
    f_bs <- abs(-2i * f_j / ka_sw) * model$body$radius / 2
    sigma_bs <- f_bs * f_bs
    c(TS = 10 * log10(sigma_bs), ml = ml, last_mod = Mod(res$last))
  }
  out <- mapply(one_freq, ka_sw, ka_l, ka_t)
  t(out)
}

freqs <- seq(1e3, 360e3, 1e3)
base <- calc_variant(freqs, m_extra = 10, adaptive = FALSE)
more <- calc_variant(freqs, m_extra = 20, adaptive = FALSE)
adpt <- calc_variant(freqs, m_extra = 10, adaptive = TRUE, tol = 1e-10)

cmp <- data.frame(
  frequency = freqs,
  base = base[, 'TS'],
  more = more[, 'TS'],
  adpt = adpt[, 'TS'],
  ml_base = base[, 'ml'],
  ml_more = more[, 'ml'],
  ml_adpt = adpt[, 'ml'],
  last_base = base[, 'last_mod'],
  last_adpt = adpt[, 'last_mod']
)

cmp$delta_more <- cmp$more - cmp$base
cmp$delta_adpt <- cmp$adpt - cmp$base
write.csv(
  cmp,
  file.path(data_dir, 'calibration_wc381_truncation_diagnostic.csv'),
  row.names = FALSE
)
cat('ABS delta more vs base\n')
print(summary(abs(cmp$delta_more)))
cat('ABS delta adaptive vs base\n')
print(summary(abs(cmp$delta_adpt)))
cat('Worst adaptive row\n')
print(cmp[which.max(abs(cmp$delta_adpt)), ])
