library(acousticTS)

.espsms_rebuild_fixture <- function() {
  ess_generate(
    shape = prolate_spheroid(
      length_body = 0.07,
      radius_body = 0.007,
      n_segments = 100
    ),
    shell_thickness = 0.0005,
    density_shell = 2565,
    density_fluid = 1077.3,
    sound_speed_fluid = 1575,
    E = 7.0e10,
    nu = 0.32
  )
}

test_that("ESPSMS rebuild geometry matches the canonical validation case", {
  obj <- .espsms_rebuild_fixture()
  body <- acousticTS:::.elastic_shell_prolate_body(
    obj,
    sound_speed_sw = 1477.3,
    density_sw = 1026.8
  )
  eta <- acousticTS:::.espsms_uniform_eta_grid(eta_grid_size = 129L)
  geom <- acousticTS:::.espsms_shell_geometry_state(body, eta)

  expect_equal(geom$semimajor_length, 0.035, tolerance = 1e-14)
  expect_equal(geom$semiminor_length, 0.007, tolerance = 1e-14)
  expect_equal(geom$a_shape, 1.0206207261596576, tolerance = 1e-12)
  expect_equal(geom$interfocal_distance, 0.068585712797929, tolerance = 1e-12)
  expect_equal(geom$bending_epsilon, 1.7006802721088435e-05, tolerance = 1e-14)
})

test_that("ESPSMS rebuild shell operator matches the external rebuild report", {
  obj <- .espsms_rebuild_fixture()
  body <- acousticTS:::.elastic_shell_prolate_body(
    obj,
    sound_speed_sw = 1477.3,
    density_sw = 1026.8
  )

  expected <- data.frame(
    frequency_hz = c(12000, 38000, 70000, 120000),
    nondimensional_frequency = c(
      0.47859179131096113,
      1.5155406724847102,
      2.7917854493139398,
      4.7859179131096115
    ),
    dynamic_fro = c(
      20485.119335691877,
      20466.865311343823,
      20418.393333759724,
      20285.572275542752
    ),
    dynamic_max_abs = c(
      5191.596010143484,
      5189.506426888276,
      5183.951350172893,
      5168.681318692528
    ),
    dynamic_relative_skew = c(
      0.06117285956880672,
      0.06122741852793515,
      0.06137276855670121,
      0.06177461061245796
    ),
    stringsAsFactors = FALSE
  )

  for (i in seq_len(nrow(expected))) {
    freq <- expected$frequency_hz[[i]]
    sys <- acousticTS:::.espsms_axisymmetric_shell_system(
      body,
      frequency_hz = freq,
      eta_grid_size = 129L
    )
    dynamic <- sys$dynamic_matrix
    rel_skew <- norm(dynamic - t(dynamic), type = "F") / norm(dynamic, type = "F")

    expect_equal(dim(dynamic), c(387L, 387L))
    expect_equal(
      sys$nondimensional_frequency,
      expected$nondimensional_frequency[[i]],
      tolerance = 1e-12
    )
    expect_equal(
      norm(dynamic, type = "F"),
      expected$dynamic_fro[[i]],
      tolerance = 1e-9
    )
    expect_equal(
      max(abs(dynamic)),
      expected$dynamic_max_abs[[i]],
      tolerance = 1e-10
    )
    expect_equal(
      rel_skew,
      expected$dynamic_relative_skew[[i]],
      tolerance = 1e-12
    )
  }
})

test_that("ESPSMS rebuild keeps the K_uw block aligned with the shell report", {
  obj <- .espsms_rebuild_fixture()
  body <- acousticTS:::.elastic_shell_prolate_body(
    obj,
    sound_speed_sw = 1477.3,
    density_sw = 1026.8
  )
  sys <- acousticTS:::.espsms_axisymmetric_shell_system(
    body,
    frequency_hz = 12000,
    eta_grid_size = 129L
  )
  kuw <- sys$structural_blocks$K_uw

  expect_equal(norm(kuw, type = "F"), 93.41762341313695, tolerance = 1e-11)
  expect_equal(
    norm(kuw - t(kuw), type = "F"),
    155.40040305023714,
    tolerance = 1e-11
  )
})
