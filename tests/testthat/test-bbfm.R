library(acousticTS)

test_that("BBF generation creates a composite body-plus-backbone scatterer", {
  body_shape <- arbitrary(
    x_body = c(0, 0.03, 0.06, 0.09),
    zU_body = c(0.0005, 0.0035, 0.0035, 0.0005),
    zL_body = c(-0.0005, -0.0035, -0.0035, -0.0005)
  )
  backbone_shape <- cylinder(
    length_body = 0.05,
    radius_body = 0.0006,
    n_segments = 20
  )

  fish <- bbf_generate(
    body_shape = body_shape,
    backbone_shape = backbone_shape,
    density_body = 1070,
    sound_speed_body = 1570,
    density_backbone = 1900,
    sound_speed_longitudinal_backbone = 3500,
    sound_speed_transversal_backbone = 1700,
    x_offset_backbone = 0.015,
    z_offset_backbone = 0.0002
  )

  expect_s4_class(fish, "BBF")
  expect_s4_class(fish, "CSC")
  expect_equal(names(extract(fish, "components")), "backbone")
  expect_equal(extract(fish, c("components", "backbone", "density")), 1900)
  expect_equal(
    extract(fish, c("shape_parameters", "backbone", "shape")),
    "Cylinder"
  )
})

test_that("BBFM matches the explicit DWBA-plus-ECMS coherent sum", {
  body_shape <- cylinder(
    length_body = 0.08,
    radius_body = 0.003,
    n_segments = 40
  )
  backbone_shape <- cylinder(
    length_body = 0.05,
    radius_body = 0.0006,
    n_segments = 20
  )
  frequency <- seq(38000, 42000, by = 2000)
  density_sw <- 1026.8
  sound_speed_sw <- 1477.3

  fish <- bbf_generate(
    body_shape = body_shape,
    backbone_shape = backbone_shape,
    density_body = 1070,
    sound_speed_body = 1570,
    density_backbone = 1900,
    sound_speed_longitudinal_backbone = 3500,
    sound_speed_transversal_backbone = 1700,
    x_offset_backbone = 0.015
  )

  fish_model <- target_strength(
    fish,
    frequency = frequency,
    model = "bbfm",
    density_sw = density_sw,
    sound_speed_sw = sound_speed_sw
  )
  composite_out <- extract(fish_model, "model")$BBFM

  body_object <- methods::new("FLS",
    metadata = list(ID = "body"),
    model_parameters = list(),
    model = list(),
    body = extract(fish, "body"),
    shape_parameters = extract(fish, c("shape_parameters", "body"))
  )
  body_object <- target_strength(
    body_object,
    frequency = frequency,
    model = "dwba",
    density_sw = density_sw,
    sound_speed_sw = sound_speed_sw
  )

  backbone_shape_params <- extract(fish, c("shape_parameters", "backbone"))
  backbone_object <- methods::new("FLS",
    metadata = list(ID = "backbone"),
    model_parameters = list(),
    model = list(),
    body = extract(fish, "backbone"),
    shape_parameters = backbone_shape_params
  )
  backbone_object <- target_strength(
    backbone_object,
    frequency = frequency,
    model = "ecms",
    density_sw = density_sw,
    sound_speed_sw = sound_speed_sw,
    density_body = extract(fish, c("backbone", "density")),
    sound_speed_longitudinal_body = extract(
      fish,
      c("backbone", "sound_speed_longitudinal")
    ),
    sound_speed_transversal_body = extract(
      fish,
      c("backbone", "sound_speed_transversal")
    )
  )

  body_out <- extract(body_object, "model")$DWBA
  backbone_out <- extract(backbone_object, "model")$ECMS
  backbone_body <- extract(backbone_object, "body")
  x_center <- mean(range(backbone_body$rpos["x", ], na.rm = TRUE))
  z_center <- mean(range(backbone_body$rpos["z", ], na.rm = TRUE))
  phase_shift <- exp(
    2i * acousticTS::wavenumber(frequency, sound_speed_sw) *
      (x_center * cos(backbone_body$theta) +
         z_center * sin(backbone_body$theta))
  )
  expected_fbs <- body_out$f_bs + backbone_out$f_bs * phase_shift

  expect_equal(composite_out$f_bs, expected_fbs, tolerance = 1e-12)
  expect_equal(composite_out$sigma_bs, abs(expected_fbs)^2, tolerance = 1e-12)
})

test_that("BBFM enforces the documented composite-input requirements", {
  density_sw <- 1026.8
  sound_speed_sw <- 1477.3
  frequency <- c(38e3, 70e3)

  valid_fish <- bbf_generate(
    body_shape = cylinder(
      length_body = 0.08,
      radius_body = 0.003,
      n_segments = 40
    ),
    backbone_shape = cylinder(
      length_body = 0.05,
      radius_body = 0.0006,
      n_segments = 20
    ),
    density_body = 1070,
    sound_speed_body = 1570,
    density_backbone = 1900,
    sound_speed_longitudinal_backbone = 3500,
    sound_speed_transversal_backbone = 1700,
    x_offset_backbone = 0.015
  )

  expect_error(
    target_strength(
      fls_generate(
        shape = cylinder(length_body = 0.08, radius_body = 0.003, n_segments = 40),
        density_body = 1070,
        sound_speed_body = 1570
      ),
      frequency = frequency,
      model = "bbfm",
      density_sw = density_sw,
      sound_speed_sw = sound_speed_sw
    ),
    "requires a 'BBF' composite scatterer"
  )

  missing_body_density <- valid_fish
  missing_body_density@body$density <- NULL
  missing_body_density@body$g <- NULL
  expect_error(
    target_strength(
      missing_body_density,
      frequency = frequency,
      model = "bbfm",
      density_sw = density_sw,
      sound_speed_sw = sound_speed_sw
    ),
    "requires the flesh body to carry either absolute density"
  )

  missing_body_speed <- valid_fish
  missing_body_speed@body$sound_speed <- NULL
  missing_body_speed@body$h <- NULL
  expect_error(
    target_strength(
      missing_body_speed,
      frequency = frequency,
      model = "bbfm",
      density_sw = density_sw,
      sound_speed_sw = sound_speed_sw
    ),
    "requires the flesh body to carry either absolute sound speed"
  )

  missing_backbone_props <- valid_fish
  missing_backbone_props@backbone$sound_speed_longitudinal <- NULL
  expect_error(
    target_strength(
      missing_backbone_props,
      frequency = frequency,
      model = "bbfm",
      density_sw = density_sw,
      sound_speed_sw = sound_speed_sw
    ),
    "requires backbone density plus longitudinal and transversal wave speeds"
  )
})

test_that("BBFM returns the documented component bookkeeping columns", {
  object <- target_strength(
    bbf_generate(
      body_shape = cylinder(
        length_body = 0.08,
        radius_body = 0.003,
        n_segments = 40
      ),
      backbone_shape = cylinder(
        length_body = 0.05,
        radius_body = 0.0006,
        n_segments = 20
      ),
      density_body = 1070,
      sound_speed_body = 1570,
      density_backbone = 1900,
      sound_speed_longitudinal_backbone = 3500,
      sound_speed_transversal_backbone = 1700,
      x_offset_backbone = 0.015
    ),
    frequency = c(38e3, 70e3),
    model = "bbfm",
    density_sw = 1026.8,
    sound_speed_sw = 1477.3
  )

  out <- extract(object, "model")$BBFM
  expect_true(all(c(
    "f_body", "f_backbone", "f_backbone_aligned", "f_bs",
    "sigma_body", "sigma_backbone", "sigma_bs",
    "TS_body", "TS_backbone", "TS"
  ) %in% names(out)))
  expect_equal(out$sigma_bs, Mod(out$f_bs)^2, tolerance = 1e-12)
})
