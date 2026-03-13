# Testing benchmark models for comparison with Jech et al. (2015)

library(acousticTS)

test_that("Test spherical model", {

  # Read in benchmark values
  data(benchmark_ts)

  # Indexing
  frequency <- benchmark_ts$frequency_spectra$index$frequency

  # MEDIUM
  density_sw <- 1026.8
  sound_speed_sw <- 1477.3

  # HELPER FUNCTION
  check_sphere_mss <- function(boundary) {
    scatterer_ess <- fixture_sphere(boundary)
    # ----> COMPUTE AND EXTRACT TS
    scatterer_ess <- target_strength(
      scatterer_ess,
      frequency = frequency,
      model = "sphms",
      boundary = boundary,
      density_sw = density_sw,
      sound_speed_sw = sound_speed_sw
    )
    # ----> RETURN TS
    expect_equal(
      acousticTS::extract(scatterer_ess, "model")$SPHMS$TS,
      benchmark_ts$frequency_spectra$sphere[[boundary]],
      tolerance = 1e-2
    )
  }

  # ESS - SPHERE - FIXED-RIGID
  check_sphere_mss(boundary="fixed_rigid")

  # ESS - SPHERE - PRESSURE RELEASE
  check_sphere_mss(boundary="pressure_release")

  # ESS - SPHERE - GAS-FILLED
  check_sphere_mss(boundary="gas_filled")

  # FLS - SPHERE - LIQUID-FILLED
  check_sphere_mss(boundary="liquid_filled")

  # ESS - SPHERE - SHELLED PRESSURE RELEASE
  check_sphere_mss(boundary = "shelled_pressure_release")

  # ESS - SPHERE - SHELLED GAS
  check_sphere_mss(boundary = "shelled_gas")

  # ESS - SPHERE - SHELLED LIQUID
  check_sphere_mss(boundary = "shelled_liquid")
})

test_that("Test cylindrical model", {

  # Read in benchmark values
  data(benchmark_ts)

  # Indexing
  frequency <- benchmark_ts$frequency_spectra$index$frequency

  # MEDIUM
  density_sw <- 1026.8
  sound_speed_sw <- 1477.3

  # HELPER FUNCTION
  check_cylinder_mss <- function(boundary) {
    scatterer_cyl <- fixture_cylinder(boundary)
    # ----> COMPUTE AND EXTRACT TS
    scatterer_cyl <- target_strength(
      scatterer_cyl,
      frequency = frequency,
      model = "fcms",
      boundary = boundary,
      density_sw = density_sw,
      sound_speed_sw = sound_speed_sw
    )
    # ----> RETURN TS
    expect_equal(
      acousticTS::extract(scatterer_cyl, "model")$FCMS$TS,
      benchmark_ts$frequency_spectra$cylinder[[boundary]],
      tolerance = 1e-2
    )
  }

  # ESS - CYLINDER - FIXED-RIGID
  check_cylinder_mss(boundary="fixed_rigid")

  # ESS - CYLINDER - PRESSURE RELEASE
  check_cylinder_mss(boundary="pressure_release")

  # ESS - CYLINDER - GAS-FILLED
  # check_cylinder_mss(boundary="gas_filled")

  # FLS - CYLINDER - LIQUID-FILLED
  check_cylinder_mss(boundary="liquid_filled")

})

test_that("Test prolate spheroid model", {

  # Read in benchmark values
  data(benchmark_ts)

  # MEDIUM
  density_sw <- 1026.8
  sound_speed_sw <- 1477.3

  # HELPER FUNCTION
  check_ps_mss <- function(boundary) {
    scatterer_ps <- fixture_ps(boundary)
    # Indexing
    frequency <- benchmark_ts$frequency_spectra$index$frequency
    if (boundary == "liquid_filled") {
      idx <- which(
        benchmark_ts$frequency_spectra$index$frequency %in%
          c(12e3, 18e3, 38e3, 100e3, 200e3)
      )
    } else {
      idx <- which(
        !is.na(benchmark_ts$frequency_spectra$prolate_spheroid[boundary])
      )
    }

    # ----> COMPUTE AND EXTRACT TS
    scatterer_ps <- target_strength(
      scatterer_ps,
      frequency = frequency[idx],
      model = "psms",
      boundary = boundary,
      simplify_Amn = FALSE,
      density_sw = density_sw,
      sound_speed_sw = sound_speed_sw,
      precision = "quad"
    )
    # ----> RETURN TS
    expect_equal(
      acousticTS::extract(scatterer_ps, "model")$PSMS$TS,
      benchmark_ts$frequency_spectra$prolate_spheroid[[boundary]][idx],
      tolerance = 1e-2
    )
  }

  # ESS - PROLATE SPHEROID - FIXED-RIGID
  check_ps_mss(boundary="fixed_rigid")

  # ESS - PROLATE SPHEROID - PRESSURE RELEASE
  check_ps_mss(boundary="pressure_release")

  # FLS - PROLATE SPHEROID - LIQUID-FILLED
  check_ps_mss(boundary="liquid_filled")

})
