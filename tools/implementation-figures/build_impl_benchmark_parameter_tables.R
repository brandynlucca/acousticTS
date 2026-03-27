library(devtools)
load_all(".", quiet = TRUE)

source("tests/testthat/helper-fixtures.R")

source("tools/implementation-figures/helpers/common.R")
impl_set_repo_root()

benchmark_csv <- "C:/Users/Brandyn/GitHub/echoSMs/src/echosms/resources/BenchMark_Data/Benchmark_Frequency_TS.csv"
jech_dir <- "C:/Users/Brandyn/GitHub/echoSMs/src/echosms/resources/Jechetal_allmodels"

bench_df <- utils::read.csv(benchmark_csv, check.names = FALSE)
density_sw <- 1026.8
sound_speed_sw <- 1477.3
sdwba_benchmark_args <- list(
  n_iterations = 100L,
  n_segments_init = 50L,
  phase_sd_init = sqrt(2) / 32,
  length_init = 38.35e-3,
  frequency_init = 120e3
)

safe_mean <- function(x) if (length(x)) mean(x, na.rm = TRUE) else NA_real_
safe_max <- function(x) if (length(x)) max(x, na.rm = TRUE) else NA_real_

compare_ts <- function(calc_ts, ref_ts) {
  ok <- is.finite(calc_ts) & is.finite(ref_ts)
  delta <- calc_ts[ok] - ref_ts[ok]
  data.frame(
    n_freq = sum(ok),
    max_abs_delta_ts_db = safe_max(abs(delta)),
    mean_abs_delta_ts_db = safe_mean(abs(delta))
  )
}

time_target_strength <- function(expr) {
  elapsed <- system.time(obj <- force(expr))[["elapsed"]]
  list(object = obj, elapsed = elapsed)
}

extract_model_df <- function(object, model_name) {
  acousticTS::extract(object, "model")[[model_name]]
}

run_model <- function(object, model_name, frequency_hz, extra_args = list()) {
  do.call(
    acousticTS::target_strength,
    c(
      list(
        object = object,
        frequency = frequency_hz,
        model = model_name,
        density_sw = density_sw,
        sound_speed_sw = sound_speed_sw
      ),
      extra_args
    )
  )
}

read_jech <- function(file) {
  utils::read.csv(file.path(jech_dir, file), check.names = FALSE)
}

make_weak_scattering <- function(shape, theta_body = pi / 2, n_segments = 201L) {
  if (shape == "sphere") {
    fls_generate(
      shape = "sphere",
      radius_body = 0.01,
      theta_body = theta_body,
      n_segments = n_segments,
      density_body = 1028.9,
      sound_speed_body = 1480.3
    )
  } else if (shape == "prolate_spheroid") {
    fls_generate(
      shape = "prolate_spheroid",
      length_body = 0.14,
      radius_body = 0.01,
      theta_body = theta_body,
      n_segments = n_segments,
      density_body = 1028.9,
      sound_speed_body = 1480.3
    )
  } else if (shape == "cylinder") {
    fls_generate(
      shape = "cylinder",
      length_body = 0.07,
      radius_body = 0.01,
      theta_body = theta_body,
      n_segments = n_segments,
      density_body = 1028.9,
      sound_speed_body = 1480.3
    )
  } else {
    stop("Unknown shape: ", shape)
  }
}

## SPHMS truncation sensitivity ------------------------------------------------

freq_all <- bench_df$Frequency_kHz * 1e3

sphms_liquid_mlimit <- do.call(rbind, lapply(
  list(NA_integer_, 20L, 10L),
  function(m_limit) {
    obj <- fixture_sphere("liquid_filled")
    args <- list(boundary = "liquid_filled")
    label <- "default"
    if (!is.na(m_limit)) {
      args$m_limit <- m_limit
      label <- as.character(m_limit)
    }
    run <- time_target_strength(run_model(obj, "sphms", freq_all, args))
    ts_df <- extract_model_df(run$object, "SPHMS")
    cmp <- compare_ts(ts_df$TS, bench_df$Sphere_WeaklyScattering)
    data.frame(
      boundary = "liquid_filled",
      m_limit = label,
      n_freq = cmp$n_freq,
      max_abs_delta_ts_db = cmp$max_abs_delta_ts_db,
      mean_abs_delta_ts_db = cmp$mean_abs_delta_ts_db,
      elapsed_sec = run$elapsed
    )
  }
))

## FCMS truncation sensitivity -------------------------------------------------

fcms_liquid_mlimit <- do.call(rbind, lapply(
  list(NA_integer_, 20L, 10L),
  function(m_limit) {
    obj <- fixture_cylinder("liquid_filled")
    args <- list(boundary = "liquid_filled")
    label <- "default"
    if (!is.na(m_limit)) {
      args$m_limit <- m_limit
      label <- as.character(m_limit)
    }
    run <- time_target_strength(run_model(obj, "fcms", freq_all, args))
    ts_df <- extract_model_df(run$object, "FCMS")
    cmp <- compare_ts(ts_df$TS, bench_df$Cylinder_WeaklyScattering)
    data.frame(
      boundary = "liquid_filled",
      m_limit = label,
      n_freq = cmp$n_freq,
      max_abs_delta_ts_db = cmp$max_abs_delta_ts_db,
      mean_abs_delta_ts_db = cmp$mean_abs_delta_ts_db,
      elapsed_sec = run$elapsed
    )
  }
))

## ESSMS truncation status -----------------------------------------------------

essms_cases <- list(
  shelled_pressure_release = fixture_sphere("shelled_pressure_release"),
  shelled_gas = fixture_sphere("shelled_gas"),
  shelled_liquid = fixture_sphere("shelled_liquid")
)

essms_mlimit_status <- do.call(rbind, lapply(names(essms_cases), function(case_name) {
  obj <- essms_cases[[case_name]]
  do.call(rbind, lapply(list(NA_integer_, 40L, 20L), function(m_limit) {
    args <- list()
    label <- "default"
    if (!is.na(m_limit)) {
      args$m_limit <- m_limit
      label <- as.character(m_limit)
    }
    run <- try(time_target_strength(run_model(obj, "essms", freq_all, args)), silent = TRUE)
    if (inherits(run, "try-error")) {
      return(data.frame(
        case = case_name,
        m_limit = label,
        benchmark_frequencies = length(freq_all),
        finite_points = 0L,
        max_abs_delta_ts_db = NA_real_,
        mean_abs_delta_ts_db = NA_real_,
        elapsed_sec = NA_real_,
        status = "error"
      ))
    }
    ts_df <- extract_model_df(run$object, "ESSMS")
    finite_idx <- is.finite(ts_df$TS)
    data.frame(
      case = case_name,
      m_limit = label,
      benchmark_frequencies = length(freq_all),
      finite_points = sum(finite_idx),
      max_abs_delta_ts_db = NA_real_,
      mean_abs_delta_ts_db = NA_real_,
      elapsed_sec = run$elapsed,
      status = if (sum(finite_idx)) "partial finite" else "no finite full-grid TS"
    )
  }))
}))

## SDWBA iteration sensitivity -------------------------------------------------

sdwba_refs <- list(
  sphere = list(
    file = "Figure_01_WeakScatt-sphere.csv",
    shape = "sphere"
  ),
  prolate_spheroid = list(
    file = "Figure_03-04_WeakScatt-Pspheroid.csv",
    shape = "prolate_spheroid"
  ),
  cylinder = list(
    file = "Figure_06_WeakScatt-Cylinder.csv",
    shape = "cylinder"
  )
)

sdwba_iterations <- do.call(rbind, lapply(names(sdwba_refs), function(case_name) {
  ref <- read_jech(sdwba_refs[[case_name]]$file)
  do.call(rbind, lapply(c(25L, 100L, 500L), function(n_iterations) {
    obj <- make_weak_scattering(sdwba_refs[[case_name]]$shape)
    set.seed(1)
    run <- time_target_strength(run_model(
      obj,
      "sdwba",
      ref$Frequency_kHz * 1e3,
      extra_args = utils::modifyList(
        sdwba_benchmark_args,
        list(n_iterations = n_iterations)
      )
    ))
    ts_df <- extract_model_df(run$object, "SDWBA")
    exact <- compare_ts(ts_df$TS, ref$Benchmark)
    published <- compare_ts(ts_df$TS, ref$Demer_SDWBA_TS)
    data.frame(
      geometry = case_name,
      n_iterations = n_iterations,
      n_freq = exact$n_freq,
      max_abs_delta_vs_benchmark_db = exact$max_abs_delta_ts_db,
      mean_abs_delta_vs_benchmark_db = exact$mean_abs_delta_ts_db,
      max_abs_delta_vs_published_db = published$max_abs_delta_ts_db,
      mean_abs_delta_vs_published_db = published$mean_abs_delta_ts_db,
      elapsed_sec = run$elapsed
    )
  }))
}))

utils::write.csv(sphms_liquid_mlimit, "scratch/sphms_mlimit_compare.csv", row.names = FALSE)
utils::write.csv(fcms_liquid_mlimit, "scratch/fcms_mlimit_compare.csv", row.names = FALSE)
utils::write.csv(essms_mlimit_status, "scratch/essms_mlimit_status.csv", row.names = FALSE)
utils::write.csv(sdwba_iterations, "scratch/sdwba_iterations_compare.csv", row.names = FALSE)

print("Done")
