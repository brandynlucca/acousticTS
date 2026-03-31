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
freq_all <- bench_df$Frequency_kHz * 1e3

read_jech <- function(file) {
  utils::read.csv(file.path(jech_dir, file), check.names = FALSE)
}

extract_model_df <- function(object, model_name) {
  acousticTS::extract(object, "model")[[model_name]]
}

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

time_target_strength <- function(expr) {
  elapsed <- system.time(obj <- force(expr))[["elapsed"]]
  list(object = obj, elapsed = elapsed)
}

detail_compare <- function(case_name, frequency_hz, calc_ts, ref_ts) {
  ok <- is.finite(calc_ts) & is.finite(ref_ts)
  data.frame(
    case = case_name,
    frequency_khz = frequency_hz[ok] * 1e-3,
    model_ts_db = calc_ts[ok],
    benchmark_ts_db = ref_ts[ok],
    delta_ts_db = calc_ts[ok] - ref_ts[ok],
    abs_delta_ts_db = abs(calc_ts[ok] - ref_ts[ok])
  )
}

top_delta <- function(df, n = 5L) {
  do.call(
    rbind,
    lapply(split(df, df$case), function(x) {
      x[order(x$abs_delta_ts_db, decreasing = TRUE), ][seq_len(min(n, nrow(x))), ]
    })
  )
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

## Modal series implementation detail tables ==================================

sphms_cases <- list(
  fixed_rigid = list(obj = fixture_sphere("fixed_rigid"), boundary = "fixed_rigid", ref = bench_df$Sphere_Rigid),
  pressure_release = list(obj = fixture_sphere("pressure_release"), boundary = "pressure_release", ref = bench_df$Sphere_PressureRelease),
  gas_filled = list(obj = fixture_sphere("gas_filled"), boundary = "gas_filled", ref = bench_df$Sphere_Gas),
  liquid_filled = list(obj = fixture_sphere("liquid_filled"), boundary = "liquid_filled", ref = bench_df$Sphere_WeaklyScattering)
)

sphms_detail <- do.call(rbind, lapply(names(sphms_cases), function(case_name) {
  info <- sphms_cases[[case_name]]
  obj <- acousticTS::target_strength(
    object = info$obj,
    frequency = freq_all,
    model = "sphms",
    boundary = info$boundary,
    density_sw = density_sw,
    sound_speed_sw = sound_speed_sw
  )
  ts_df <- extract_model_df(obj, "SPHMS")
  detail_compare(case_name, freq_all, ts_df$TS, info$ref)
}))

fcms_cases <- list(
  fixed_rigid = list(obj = fixture_cylinder("fixed_rigid"), boundary = "fixed_rigid", ref = bench_df$Cylinder_Rigid),
  pressure_release = list(obj = fixture_cylinder("pressure_release"), boundary = "pressure_release", ref = bench_df$Cylinder_PressureRelease),
  gas_filled = list(obj = fixture_cylinder("gas_filled"), boundary = "gas_filled", ref = bench_df$Cylinder_Gas),
  liquid_filled = list(obj = fixture_cylinder("liquid_filled"), boundary = "liquid_filled", ref = bench_df$Cylinder_WeaklyScattering)
)

fcms_detail <- do.call(rbind, lapply(names(fcms_cases), function(case_name) {
  info <- fcms_cases[[case_name]]
  obj <- acousticTS::target_strength(
    object = info$obj,
    frequency = freq_all,
    model = "fcms",
    boundary = info$boundary,
    density_sw = density_sw,
    sound_speed_sw = sound_speed_sw
  )
  ts_df <- extract_model_df(obj, "FCMS")
  detail_compare(case_name, freq_all, ts_df$TS, info$ref)
}))

hpa_cases <- list(
  sphere_johnson = list(shape = "sphere", method = "johnson", ref = bench_df$Sphere_WeaklyScattering),
  sphere_stanton = list(shape = "sphere", method = "stanton", ref = bench_df$Sphere_WeaklyScattering),
  prolate_spheroid_stanton = list(shape = "prolate_spheroid", method = "stanton", ref = bench_df$ProlateSpheroid_WeaklyScattering),
  cylinder_stanton = list(shape = "cylinder", method = "stanton", ref = bench_df$Cylinder_WeaklyScattering)
)

hpa_detail <- do.call(rbind, lapply(names(hpa_cases), function(case_name) {
  info <- hpa_cases[[case_name]]
  obj <- run_simple_model(
    make_weak_scattering(info$shape),
    "hpa",
    freq_all,
    extra_args = list(method = info$method)
  )
  ts_df <- extract_model_df(obj, "HPA")
  detail_compare(case_name, freq_all, ts_df$TS, info$ref)
}))

## Approximation-family detail tables =========================================

dwba_refs <- list(
  sphere = list(file = "Figure_01_WeakScatt-sphere.csv", shape = "sphere"),
  prolate_spheroid = list(file = "Figure_03-04_WeakScatt-Pspheroid.csv", shape = "prolate_spheroid"),
  cylinder = list(file = "Figure_06_WeakScatt-Cylinder.csv", shape = "cylinder")
)

dwba_detail <- do.call(rbind, lapply(names(dwba_refs), function(case_name) {
  ref <- read_jech(dwba_refs[[case_name]]$file)
  obj <- run_simple_model(
    make_dwba_reference_shape(dwba_refs[[case_name]]$shape),
    "dwba",
    ref$Frequency_kHz * 1e3
  )
  ts_df <- extract_model_df(obj, "DWBA")
  detail_compare(case_name, ref$Frequency_kHz * 1e3, ts_df$TS, ref$Benchmark)
}))

sdwba_refs <- dwba_refs[c("sphere", "prolate_spheroid", "cylinder")]

set.seed(1)
sdwba_detail <- do.call(rbind, lapply(names(sdwba_refs), function(case_name) {
  ref <- read_jech(sdwba_refs[[case_name]]$file)
  obj <- run_simple_model(
    make_weak_scattering(sdwba_refs[[case_name]]$shape),
    "sdwba",
    ref$Frequency_kHz * 1e3
  )
  ts_df <- extract_model_df(obj, "SDWBA")
  detail_compare(case_name, ref$Frequency_kHz * 1e3, ts_df$TS, ref$Benchmark)
}))

krm_cases <- list(
  sphere_gas = list(shape = "sphere", file = "Figure_01_Gas-sphere.csv"),
  sphere_weak = list(shape = "sphere", file = "Figure_01_WeakScatt-sphere.csv"),
  prolate_gas = list(shape = "prolate_spheroid", file = "Figure_03-04_Gas-Pspheroid.csv"),
  prolate_weak = list(shape = "prolate_spheroid", file = "Figure_03-04_WeakScatt-Pspheroid.csv"),
  cylinder_gas = list(shape = "cylinder", file = "Figure_06_Gas-Cylinder.csv"),
  cylinder_weak = list(shape = "cylinder", file = "Figure_06_WeakScatt-Cylinder.csv")
)

krm_detail <- do.call(rbind, lapply(names(krm_cases), function(case_name) {
  info <- krm_cases[[case_name]]
  kind <- if (grepl("_gas$", case_name)) "gas" else "weak"
  ref <- read_jech(info$file)
  obj <- run_simple_model(
    make_krm_reference_shape(info$shape, kind),
    "krm",
    ref$Frequency_kHz * 1e3
  )
  ts_df <- extract_model_df(obj, "KRM")
  detail_compare(case_name, ref$Frequency_kHz * 1e3, ts_df$TS, ref$Benchmark)
}))

cal_ref <- data.frame(
  frequency = c(38e3, 70e3, 120e3, 200e3),
  ts = c(-42.3296920811525, -41.07033000976835, -39.50264169247166, -39.43821198162192)
)

cal_obj <- run_simple_model(cal_generate(), "calibration", cal_ref$frequency)
cal_df <- extract_model_df(cal_obj, "calibration")
calibration_detail <- detail_compare("default_WC_38.1mm", cal_ref$frequency, cal_df$TS, cal_ref$ts)

## Write detail and top-delta tables ==========================================

write_pair <- function(name, df) {
  utils::write.csv(df, file.path("scratch", paste0(name, "_benchmark_detail.csv")), row.names = FALSE)
  utils::write.csv(top_delta(df), file.path("scratch", paste0(name, "_benchmark_topdelta.csv")), row.names = FALSE)
}

write_pair("sphms", sphms_detail)
write_pair("fcms", fcms_detail)
write_pair("hpa", hpa_detail)
write_pair("dwba", dwba_detail)
write_pair("sdwba", sdwba_detail)
write_pair("krm", krm_detail)
write_pair("calibration", calibration_detail)
