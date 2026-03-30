library(acousticTS)

test_that("direct PSMS kernels cover fluid vectorization and exact rigid or soft branches", {
  density_sw <- 1026.8
  sound_speed_sw <- 1477.3
  quad32 <- gauss_legendre(n = 32, a = -1, b = 1)

  fluid_init <- acousticTS:::psms_initialize(
    object = fixture_ps("liquid_filled"),
    frequency = c(12e3, 18e3),
    boundary = "liquid_filled",
    adaptive = FALSE,
    precision = "double",
    n_integration = 32L,
    simplify_Amn = FALSE,
    sound_speed_sw = sound_speed_sw,
    density_sw = density_sw
  )
  fluid_params <- fluid_init@model_parameters[[1]]
  fluid_fixed <- acousticTS:::prolate_spheroid_fbs(
    acoustics = fluid_params$parameters$acoustics,
    body = fluid_params$body,
    medium = fluid_params$medium,
    integration_pts = quad32,
    precision = "double",
    Amn_method = "Amn_fluid",
    adaptive = FALSE,
    vectorized = FALSE
  )
  fluid_vectorized <- acousticTS:::prolate_spheroid_fbs(
    acoustics = fluid_params$parameters$acoustics,
    body = fluid_params$body,
    medium = fluid_params$medium,
    integration_pts = quad32,
    precision = "double",
    Amn_method = "Amn_fluid",
    adaptive = FALSE,
    vectorized = TRUE
  )
  fluid_simplified <- acousticTS:::prolate_spheroid_fbs(
    acoustics = fluid_params$parameters$acoustics,
    body = fluid_params$body,
    medium = fluid_params$medium,
    integration_pts = quad32,
    precision = "double",
    Amn_method = "Amn_fluid_simplify",
    adaptive = FALSE,
    vectorized = FALSE
  )
  fluid_adaptive <- acousticTS:::prolate_spheroid_fbs(
    acoustics = fluid_params$parameters$acoustics,
    body = fluid_params$body,
    medium = fluid_params$medium,
    integration_pts = quad32,
    precision = "double",
    Amn_method = "Amn_fluid",
    adaptive = TRUE,
    vectorized = TRUE
  )

  expect_equal(fluid_vectorized, fluid_fixed, tolerance = 1e-8)
  expect_length(fluid_simplified, 2L)
  expect_true(all(is.finite(fluid_simplified)))
  expect_length(fluid_adaptive, 2L)
  expect_true(all(is.finite(fluid_adaptive)))

  rigid_init <- acousticTS:::psms_initialize(
    object = fixture_ps("fixed_rigid"),
    frequency = 38e3,
    boundary = "fixed_rigid",
    adaptive = TRUE,
    precision = "quad",
    sound_speed_sw = sound_speed_sw,
    density_sw = density_sw
  )
  rigid_params <- rigid_init@model_parameters[[1]]
  rigid_direct <- acousticTS:::prolate_spheroid_fbs(
    acoustics = rigid_params$parameters$acoustics,
    body = rigid_params$body,
    medium = rigid_params$medium,
    integration_pts = quad32,
    precision = "quad",
    Amn_method = "Amn_fixed_rigid",
    adaptive = TRUE,
    vectorized = FALSE
  )
  rigid_public <- acousticTS:::PSMS(rigid_init)@model$PSMS$f_bs
  expect_equal(rigid_direct, rigid_public, tolerance = 1e-10)

  soft_init <- acousticTS:::psms_initialize(
    object = fixture_ps("pressure_release"),
    frequency = 38e3,
    boundary = "pressure_release",
    adaptive = TRUE,
    precision = "double",
    sound_speed_sw = sound_speed_sw,
    density_sw = density_sw
  )
  soft_params <- soft_init@model_parameters[[1]]
  soft_direct <- acousticTS:::prolate_spheroid_fbs(
    acoustics = soft_params$parameters$acoustics,
    body = soft_params$body,
    medium = soft_params$medium,
    integration_pts = quad32,
    precision = "double",
    Amn_method = "Amn_pressure_release",
    adaptive = TRUE,
    vectorized = FALSE
  )
  soft_public <- acousticTS:::PSMS(soft_init)@model$PSMS$f_bs
  expect_equal(soft_direct, soft_public, tolerance = 1e-10)
})

test_that("direct retained prolate TMM kernels match stored double-precision helpers", {
  density_sw <- 1026.8
  sound_speed_sw <- 1480
  theta_body <- pi / 2
  phi_body <- pi / 2
  theta_scatter <- pi / 3
  phi_scatter <- pi
  theta_grid <- c(pi / 3, pi / 2)
  phi_grid <- c(0, pi)
  quad96 <- gauss_legendre(n = 96, a = -1, b = 1)

  object <- target_strength(
    object = fls_generate(
      shape = prolate_spheroid(length_body = 0.07, radius_body = 0.01, n_segments = 80),
      g_body = 1,
      h_body = 1,
      theta_body = theta_body
    ),
    frequency = 38e3,
    model = "tmm",
    boundary = "pressure_release",
    density_sw = density_sw,
    sound_speed_sw = sound_speed_sw,
    store_t_matrix = TRUE
  )
  params <- object@model_parameters[[1]]

  direct_tmm <- acousticTS:::prolate_spheroid_tmatrix_cpp(
    acoustics = params$parameters$acoustics,
    body = params$body,
    medium = params$medium,
    integration_pts = quad96,
    precision = params$parameters$precision,
    Amn_method = "Amn_pressure_release"
  )
  expect_length(
    direct_tmm$t_matrix[[1]],
    params$parameters$acoustics$m_max[1] + 1L
  )
  expect_equal(
    (-2i / params$parameters$acoustics$k_sw) * direct_tmm$f_scat,
    object@model$TMM$f_bs,
    tolerance = 1e-10
  )

  direct_point <- acousticTS:::prolate_spheroid_scattering_from_tmatrix_cpp(
    acoustics = params$parameters$acoustics,
    t_matrix = params$parameters$t_matrix,
    theta_body = theta_body,
    phi_body = phi_body,
    theta_scatter = theta_scatter,
    phi_scatter = phi_scatter,
    precision = params$parameters$precision
  )
  public_point <- tmm_scattering(
    object,
    theta_body = theta_body,
    phi_body = phi_body,
    theta_scatter = theta_scatter,
    phi_scatter = phi_scatter
  )$f_scat
  expect_equal(direct_point, public_point, tolerance = 1e-10)

  direct_grid <- acousticTS:::prolate_spheroid_scattering_grid_from_tmatrix_cpp(
    acoustics = params$parameters$acoustics,
    t_matrix = params$parameters$t_matrix[[1]],
    theta_body = theta_body,
    phi_body = phi_body,
    theta_scatter = theta_grid,
    phi_scatter = phi_grid,
    precision = params$parameters$precision
  )
  public_grid <- tmm_scattering_grid(
    object,
    frequency = 38e3,
    theta_body = theta_body,
    phi_body = phi_body,
    theta_scatter = theta_grid,
    phi_scatter = phi_grid
  )$f_scat
  expect_equal(direct_grid, public_grid, tolerance = 1e-10)
})

test_that("direct retained prolate TMM kernels match stored quad-precision helpers", {
  density_sw <- 1026.8
  sound_speed_sw <- 1477.3
  theta_body <- pi / 2
  phi_body <- pi / 2
  theta_scatter <- pi / 2
  phi_scatter <- 3 * pi / 2
  theta_grid <- c(pi / 3, pi / 2)
  phi_grid <- c(pi / 2, 3 * pi / 2)
  quad96 <- gauss_legendre(n = 96, a = -1, b = 1)

  object <- target_strength(
    object = fls_generate(
      shape = prolate_spheroid(length_body = 0.07, radius_body = 0.01, n_segments = 80),
      density_body = 1028.9,
      sound_speed_body = 1480.3,
      theta_body = theta_body
    ),
    frequency = 38e3,
    model = "tmm",
    boundary = "liquid_filled",
    density_sw = density_sw,
    sound_speed_sw = sound_speed_sw,
    store_t_matrix = TRUE
  )
  params <- object@model_parameters[[1]]

  direct_tmm <- acousticTS:::prolate_spheroid_tmatrix_cpp(
    acoustics = params$parameters$acoustics,
    body = params$body,
    medium = params$medium,
    integration_pts = quad96,
    precision = params$parameters$precision,
    Amn_method = "Amn_fluid"
  )
  expect_equal(
    (-2i / params$parameters$acoustics$k_sw) * direct_tmm$f_scat,
    object@model$TMM$f_bs,
    tolerance = 1e-8
  )

  direct_point <- acousticTS:::prolate_spheroid_scattering_from_tmatrix_cpp(
    acoustics = params$parameters$acoustics,
    t_matrix = params$parameters$t_matrix,
    theta_body = theta_body,
    phi_body = phi_body,
    theta_scatter = theta_scatter,
    phi_scatter = phi_scatter,
    precision = params$parameters$precision
  )
  public_point <- tmm_scattering(
    object,
    theta_body = theta_body,
    phi_body = phi_body,
    theta_scatter = theta_scatter,
    phi_scatter = phi_scatter
  )$f_scat
  expect_equal(direct_point, public_point, tolerance = 1e-8)

  direct_grid <- acousticTS:::prolate_spheroid_scattering_grid_from_tmatrix_cpp(
    acoustics = params$parameters$acoustics,
    t_matrix = params$parameters$t_matrix[[1]],
    theta_body = theta_body,
    phi_body = phi_body,
    theta_scatter = theta_grid,
    phi_scatter = phi_grid,
    precision = params$parameters$precision
  )
  public_grid <- tmm_scattering_grid(
    object,
    frequency = 38e3,
    theta_body = theta_body,
    phi_body = phi_body,
    theta_scatter = theta_grid,
    phi_scatter = phi_grid
  )$f_scat
  expect_equal(direct_grid, public_grid, tolerance = 1e-8)
})

test_that("direct liquid-filled exact spheroidal field matches stored TMM general-angle output", {
  density_sw <- 1026.8
  sound_speed_sw <- 1477.3
  theta_body <- pi / 2
  phi_body <- pi / 2
  scatter_angles <- data.frame(
    theta_scatter = c(pi / 3, pi / 2),
    phi_scatter = c(pi, 3 * pi / 2)
  )
  quad64 <- gauss_legendre(n = 64, a = -1, b = 1)

  object <- target_strength(
    object = fls_generate(
      shape = prolate_spheroid(length_body = 0.07, radius_body = 0.01, n_segments = 80),
      density_body = 1028.9,
      sound_speed_body = 1480.3,
      theta_body = theta_body
    ),
    frequency = 18e3,
    model = "tmm",
    boundary = "liquid_filled",
    density_sw = density_sw,
    sound_speed_sw = sound_speed_sw,
    store_t_matrix = TRUE
  )
  params <- object@model_parameters[[1]]
  acoustics <- params$parameters$acoustics
  body <- params$body
  medium <- params$medium

  for (i in seq_len(nrow(scatter_angles))) {
    body$theta_body <- theta_body
    body$phi_body <- phi_body
    body$theta_scatter <- scatter_angles$theta_scatter[i]
    body$phi_scatter <- scatter_angles$phi_scatter[i]

    exact_raw <- acousticTS:::prolate_spheroid_fbs(
      acoustics = acoustics,
      body = body,
      medium = medium,
      integration_pts = quad64,
      precision = "quad",
      Amn_method = "Amn_fluid",
      adaptive = FALSE,
      vectorized = FALSE
    )
    exact_f <- (-2i / acoustics$k_sw) * exact_raw
    stored_f <- tmm_scattering(
      object,
      theta_body = theta_body,
      phi_body = phi_body,
      theta_scatter = body$theta_scatter,
      phi_scatter = body$phi_scatter
    )$f_scat

    expect_equal(stored_f, exact_f, tolerance = 1e-7)
  }

  simplified_raw <- acousticTS:::prolate_spheroid_fbs(
    acoustics = acoustics,
    body = body,
    medium = medium,
    integration_pts = quad64,
    precision = "double",
    Amn_method = "Amn_fluid_simplify",
    adaptive = FALSE,
    vectorized = FALSE
  )
  expect_true(all(is.finite(simplified_raw)))
})
