source("tools/implementation-figures/helpers/common.R")
impl_load_all()

source("tests/testthat/helper-fixtures.R")
options(error = function() {
  traceback(2)
  quit(save = "no", status = 1)
})

benchmark_csv <- "C:/Users/Brandyn/GitHub/echoSMs/src/echosms/resources/BenchMark_Data/Benchmark_Frequency_TS.csv"
jech_dir <- "C:/Users/Brandyn/GitHub/echoSMs/src/echosms/resources/Jechetal_allmodels"
noaa_dir <- "C:/Users/Brandyn/GitHub/echoSMs/src/echosms/resources"

bench_df <- utils::read.csv(benchmark_csv, check.names = FALSE)

density_sw <- 1026.8
sound_speed_sw <- 1477.3
sdwba_benchmark_args <- list(
  n_iterations = 100,
  n_segments_init = 50,
  phase_sd_init = sqrt(2) / 32,
  length_init = 38.35e-3,
  frequency_init = 120e3
)

freq_all <- bench_df$Frequency_kHz * 1e3

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

read_jech <- function(file) {
  utils::read.csv(file.path(jech_dir, file), check.names = FALSE)
}

read_noaa_krm <- function(file) {
  df <- utils::read.csv(
    file.path(noaa_dir, file),
    check.names = FALSE,
    row.names = NULL
  )
  out <- data.frame(
    frequency_khz = suppressWarnings(as.numeric(df[[5]])),
    ts_db = suppressWarnings(as.numeric(df[[6]]))
  )
  out[stats::complete.cases(out), ]
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

make_dwba_reference_shape <- function(shape, theta_body = pi / 2, spacing = 1e-4) {
  n_segments <- max(2L, round(switch(shape,
    sphere = (2 * 0.01) / spacing,
    prolate_spheroid = (2 * 0.07) / spacing,
    cylinder = 0.07 / spacing,
    stop("Unknown DWBA reference shape: ", shape)
  )))

  shape_obj <- if (shape == "sphere") {
    sphere(radius_body = 0.01, n_segments = n_segments)
  } else if (shape == "prolate_spheroid") {
    prolate_spheroid(
      length_body = 0.14,
      radius_body = 0.01,
      n_segments = n_segments
    )
  } else if (shape == "cylinder") {
    cylinder(
      length_body = 0.07,
      radius_body = 0.01,
      n_segments = n_segments
    )
  } else {
    stop("Unknown DWBA reference shape: ", shape)
  }

  fls_generate(
    shape = shape_obj,
    theta_body = theta_body,
    density_body = 1028.9,
    sound_speed_body = 1480.3
  )
}

make_trcm_bent_reference <- function(
  ka = seq(0.1, 10, length.out = 195L),
  length_body = 10.5e-3,
  radius_body = 1e-3,
  radius_curvature_ratio = 1.5,
  g_body = 1.0357,
  h_body = 1.0279,
  n_segments = 401L
) {
  frequency <- ka * sound_speed_sw / (2 * pi * radius_body)
  body_args <- list(
    density_body = density_sw * g_body,
    sound_speed_body = sound_speed_sw * h_body,
    theta_body = pi / 2
  )

  straight <- do.call(
    fls_generate,
    c(
      list(shape = cylinder(
        length_body = length_body,
        radius_body = radius_body,
        n_segments = n_segments
      )),
      body_args
    )
  )

  bent <- do.call(
    fls_generate,
    c(
      list(shape = cylinder(
        length_body = length_body,
        radius_body = radius_body,
        radius_curvature_ratio = radius_curvature_ratio,
        n_segments = n_segments
      )),
      body_args
    )
  )

  straight_fcms <- target_strength(
    straight,
    frequency,
    "fcms",
    boundary = "liquid_filled",
    density_sw = density_sw,
    sound_speed_sw = sound_speed_sw
  )

  lebc <- .trcm_equivalent_length_fresnel(
    k1 = 2 * pi * frequency / sound_speed_sw,
    l = length_body,
    a = radius_body,
    rho_c = radius_curvature_ratio * length_body
  )

  list(
    ka = ka,
    frequency = frequency,
    reference_ts = 20 * log10(abs(lebc * straight_fcms@model$FCMS$f_bs / length_body)),
    object = bent
  )
}

## -------------------------- Exact modal series ------------------------------

run_modal_case <- function(object, model_name, boundary = NULL, frequency_hz, ref_ts, extra_args = list()) {
  run <- time_target_strength(do.call(
    acousticTS::target_strength,
    c(
      list(
        object = object,
        frequency = frequency_hz,
        model = model_name,
        density_sw = density_sw,
        sound_speed_sw = sound_speed_sw
      ),
      if (!is.null(boundary)) list(boundary = boundary) else list(),
      extra_args
    )
  ))
  mod_name <- toupper(model_name)
  if (model_name == "calibration") mod_name <- "calibration"
  ts_df <- extract_model_df(run$object, mod_name)
  cbind(
    compare_ts(ts_df$TS, ref_ts),
    elapsed_sec = run$elapsed
  )
}

sphms_cases <- list(
  fixed_rigid = list(obj = fixture_sphere("fixed_rigid"), boundary = "fixed_rigid", ref = bench_df$Sphere_Rigid),
  pressure_release = list(obj = fixture_sphere("pressure_release"), boundary = "pressure_release", ref = bench_df$Sphere_PressureRelease),
  gas_filled = list(obj = fixture_sphere("gas_filled"), boundary = "gas_filled", ref = bench_df$Sphere_Gas),
  liquid_filled = list(obj = fixture_sphere("liquid_filled"), boundary = "liquid_filled", ref = bench_df$Sphere_WeaklyScattering)
)

message("Running SPHMS")
sphms_summary <- do.call(rbind, lapply(names(sphms_cases), function(case_name) {
  info <- sphms_cases[[case_name]]
  cbind(
    case = case_name,
    run_modal_case(info$obj, "sphms", info$boundary, freq_all, info$ref)
  )
}))

fcms_cases <- list(
  fixed_rigid = list(obj = fixture_cylinder("fixed_rigid"), boundary = "fixed_rigid", ref = bench_df$Cylinder_Rigid),
  pressure_release = list(obj = fixture_cylinder("pressure_release"), boundary = "pressure_release", ref = bench_df$Cylinder_PressureRelease),
  gas_filled = list(obj = fixture_cylinder("gas_filled"), boundary = "gas_filled", ref = bench_df$Cylinder_Gas),
  liquid_filled = list(obj = fixture_cylinder("liquid_filled"), boundary = "liquid_filled", ref = bench_df$Cylinder_WeaklyScattering)
)

message("Running FCMS")
fcms_summary <- do.call(rbind, lapply(names(fcms_cases), function(case_name) {
  info <- fcms_cases[[case_name]]
  cbind(
    case = case_name,
    run_modal_case(info$obj, "fcms", info$boundary, freq_all, info$ref)
  )
}))

essms_cases <- list(
  shelled_pressure_release = list(obj = fixture_sphere("shelled_pressure_release"), ref = bench_df$ShellSphere_PressureRelease),
  shelled_gas = list(obj = fixture_sphere("shelled_gas"), ref = bench_df$ShellSphere_Gas),
  shelled_liquid = list(obj = fixture_sphere("shelled_liquid"), ref = bench_df$ShellSphere_WeaklyScattering)
)

message("Running ESSMS")
essms_summary <- do.call(rbind, lapply(names(essms_cases), function(case_name) {
  info <- essms_cases[[case_name]]
  res <- try(run_modal_case(info$obj, "essms", NULL, freq_all, info$ref), silent = TRUE)
  if (inherits(res, "try-error")) {
    data.frame(
      case = case_name,
      n_freq = 0L,
      max_abs_delta_ts_db = NA_real_,
      mean_abs_delta_ts_db = NA_real_,
      elapsed_sec = NA_real_,
      status = "error"
    )
  } else {
    cbind(case = case_name, res, status = "ok")
  }
}))

## ------------------------- Approximation families ---------------------------

run_simple_model <- function(object, model_name, frequency_hz, extra_args = list()) {
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

dwba_refs <- list(
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

message("Running DWBA")
dwba_summary <- do.call(rbind, lapply(names(dwba_refs), function(case_name) {
  ref <- read_jech(dwba_refs[[case_name]]$file)
  obj <- make_dwba_reference_shape(dwba_refs[[case_name]]$shape)
  run <- time_target_strength(run_simple_model(obj, "dwba", ref$Frequency_kHz * 1e3))
  ts_df <- extract_model_df(run$object, "DWBA")
  exact <- compare_ts(ts_df$TS, ref$Benchmark)
  data.frame(
    geometry = case_name,
    max_abs_delta_ts_db = exact$max_abs_delta_ts_db,
    mean_abs_delta_ts_db = exact$mean_abs_delta_ts_db,
    elapsed_sec = run$elapsed
  )
}))

message("Running SDWBA")
sdwba_refs <- dwba_refs[c("sphere", "prolate_spheroid", "cylinder")]

sdwba_summary <- do.call(rbind, lapply(names(sdwba_refs), function(case_name) {
  ref <- read_jech(dwba_refs[[case_name]]$file)
  obj <- make_weak_scattering(sdwba_refs[[case_name]]$shape)
  set.seed(1)
  run <- time_target_strength(run_simple_model(
    obj,
    "sdwba",
    ref$Frequency_kHz * 1e3,
    extra_args = sdwba_benchmark_args
  ))
  ts_df <- extract_model_df(run$object, "SDWBA")
  exact <- compare_ts(ts_df$TS, ref$Benchmark)
  data.frame(
    geometry = case_name,
    max_abs_delta_ts_db = exact$max_abs_delta_ts_db,
    mean_abs_delta_ts_db = exact$mean_abs_delta_ts_db,
    elapsed_sec = run$elapsed
  )
}))

hpa_cases <- list(
  sphere_johnson = list(shape = "sphere", method = "johnson", ref = bench_df$Sphere_WeaklyScattering),
  sphere_stanton = list(shape = "sphere", method = "stanton", ref = bench_df$Sphere_WeaklyScattering),
  prolate_spheroid_stanton = list(shape = "prolate_spheroid", method = "stanton", ref = bench_df$ProlateSpheroid_WeaklyScattering),
  cylinder_stanton = list(shape = "cylinder", method = "stanton", ref = bench_df$Cylinder_WeaklyScattering)
)

message("Running HPA")
hpa_summary <- do.call(rbind, lapply(names(hpa_cases), function(case_name) {
  info <- hpa_cases[[case_name]]
  obj <- make_weak_scattering(info$shape)
  run <- time_target_strength(run_simple_model(
    obj,
    "hpa",
    freq_all,
    extra_args = list(method = info$method)
  ))
  ts_df <- extract_model_df(run$object, "HPA")
  exact <- compare_ts(ts_df$TS, info$ref)
  data.frame(
    case = case_name,
    n_freq = exact$n_freq,
    max_abs_delta_ts_db = exact$max_abs_delta_ts_db,
    mean_abs_delta_ts_db = exact$mean_abs_delta_ts_db,
    elapsed_sec = run$elapsed
  )
}))

message("Running TRCM")
trcm_ref_straight <- read_jech("Figure_06_WeakScatt-Cylinder.csv")
trcm_straight_run <- time_target_strength(run_simple_model(
  make_weak_scattering("cylinder"),
  "trcm",
  trcm_ref_straight$Frequency_kHz * 1e3
))
trcm_straight_df <- extract_model_df(trcm_straight_run$object, "TRCM")
trcm_straight_cmp <- compare_ts(trcm_straight_df$TS, trcm_ref_straight$Benchmark)

trcm_bent_ref <- make_trcm_bent_reference()
trcm_bent_summary <- do.call(rbind, lapply(c(FALSE, TRUE), function(stationary_phase) {
  run <- time_target_strength(run_simple_model(
    trcm_bent_ref$object,
    "trcm",
    trcm_bent_ref$frequency,
    extra_args = list(stationary_phase = stationary_phase)
  ))
  ts_df <- extract_model_df(run$object, "TRCM")
  exact <- compare_ts(ts_df$TS, trcm_bent_ref$reference_ts)
  data.frame(
    geometry = "bent_cylinder",
    implementation = if (stationary_phase) "stationary_phase" else "fresnel_integral",
    reference = "modal_series_bent_cylinder",
    n_freq = exact$n_freq,
    max_abs_delta_ts_db = exact$max_abs_delta_ts_db,
    mean_abs_delta_ts_db = exact$mean_abs_delta_ts_db,
    elapsed_sec = run$elapsed
  )
}))

trcm_summary <- rbind(
  data.frame(
    geometry = "straight_cylinder",
    implementation = "standard",
    reference = "weak_cylinder_benchmark",
    n_freq = trcm_straight_cmp$n_freq,
    max_abs_delta_ts_db = trcm_straight_cmp$max_abs_delta_ts_db,
    mean_abs_delta_ts_db = trcm_straight_cmp$mean_abs_delta_ts_db,
    elapsed_sec = trcm_straight_run$elapsed
  ),
  trcm_bent_summary
)

## ----------------------------- KRM ------------------------------------------

message("Running KRM")
make_krm_reference_shape <- function(shape, kind, theta_body = pi / 2,
                                     n_segments = 200L) {
  shape_obj <- switch(shape,
    sphere = sphere(radius_body = 0.01, n_segments = n_segments),
    prolate_spheroid = prolate_spheroid(
      length_body = 0.14,
      radius_body = 0.01,
      n_segments = n_segments
    ),
    cylinder = cylinder(
      length_body = 0.07,
      radius_body = 0.01,
      n_segments = n_segments
    )
  )

  props <- if (kind == "gas") {
    list(density_body = 1.24, sound_speed_body = 345)
  } else {
    list(density_body = 1028.9, sound_speed_body = 1480.3)
  }

  do.call(
    fls_generate,
    c(list(shape = shape_obj, theta_body = theta_body), props)
  )
}

krm_cases <- list(
  sphere_gas = list(shape = "sphere", file = "Figure_01_Gas-sphere.csv"),
  sphere_weak = list(shape = "sphere", file = "Figure_01_WeakScatt-sphere.csv"),
  prolate_gas = list(shape = "prolate_spheroid", file = "Figure_03-04_Gas-Pspheroid.csv"),
  prolate_weak = list(shape = "prolate_spheroid", file = "Figure_03-04_WeakScatt-Pspheroid.csv"),
  cylinder_gas = list(shape = "cylinder", file = "Figure_06_Gas-Cylinder.csv"),
  cylinder_weak = list(shape = "cylinder", file = "Figure_06_WeakScatt-Cylinder.csv")
)

krm_summary <- do.call(rbind, lapply(names(krm_cases), function(case_name) {
  info <- krm_cases[[case_name]]
  kind <- if (grepl("_gas$", case_name)) "gas" else "weak"
  ref <- read_jech(info$file)
  run <- time_target_strength(run_simple_model(
    make_krm_reference_shape(info$shape, kind),
    "krm",
    ref$Frequency_kHz * 1e3
  ))
  ts_df <- extract_model_df(run$object, "KRM")
  exact <- compare_ts(ts_df$TS, ref$Benchmark)
  data.frame(
    case = case_name,
    n_freq = exact$n_freq,
    max_abs_delta_ts_db = exact$max_abs_delta_ts_db,
    mean_abs_delta_ts_db = exact$mean_abs_delta_ts_db,
    elapsed_sec = run$elapsed
  )
}))

## -------------------------- Calibration -------------------------------------

message("Running calibration")
cal_ref <- data.frame(
  frequency = c(38e3, 70e3, 120e3, 200e3),
  ts = c(-42.3296920811525, -41.07033000976835, -39.50264169247166, -39.43821198162192)
)

cal_obj <- cal_generate()
cal_run <- time_target_strength(run_simple_model(cal_obj, "calibration", cal_ref$frequency))
cal_df <- extract_model_df(cal_run$object, "calibration")
cal_cmp <- compare_ts(cal_df$TS, cal_ref$ts)
calibration_summary <- data.frame(
  case = "default_WC_38.1mm",
  n_freq = cal_cmp$n_freq,
  max_abs_delta_ts_db = cal_cmp$max_abs_delta_ts_db,
  mean_abs_delta_ts_db = cal_cmp$mean_abs_delta_ts_db,
  elapsed_sec = cal_run$elapsed
)

utils::write.csv(sphms_summary, "scratch/sphms_benchmark_summary_full.csv", row.names = FALSE)
utils::write.csv(fcms_summary, "scratch/fcms_benchmark_summary_full.csv", row.names = FALSE)
utils::write.csv(essms_summary, "scratch/essms_benchmark_summary_full.csv", row.names = FALSE)
utils::write.csv(dwba_summary, "scratch/dwba_benchmark_summary_full.csv", row.names = FALSE)
utils::write.csv(sdwba_summary, "scratch/sdwba_benchmark_summary_full.csv", row.names = FALSE)
utils::write.csv(hpa_summary, "scratch/hpa_benchmark_summary_full.csv", row.names = FALSE)
utils::write.csv(trcm_summary, "scratch/trcm_benchmark_summary_full.csv", row.names = FALSE)
utils::write.csv(krm_summary, "scratch/krm_benchmark_summary_full.csv", row.names = FALSE)
utils::write.csv(calibration_summary, "scratch/calibration_benchmark_summary_full.csv", row.names = FALSE)

print("Done")
