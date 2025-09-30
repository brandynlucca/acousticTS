library(acousticTS)

test_that("Simple acoustic utility functions work as intended", {
  # Test linear
  expect_equal(linear(-70), 1e-7)

  # Test db
  expect_equal(db(1e-7), -70)

  # Test sigma_bs
  expected <- 2e-20
  expect_equal(sigma_bs(1e-10 * 1e-10 * 1i), expected)

  # Test transmission_coefficient
  expected <- 0.957670166
  expect_equal(
    transmission_coefficient(
      data.frame(sound_speed = 1500, density = 1026),
      data.frame(sound_speed = 1400, density = 1010)
    ),
    expected,
    tolerance = 1e-10
  )

  # Test kappa
  expected <- 0.1661446757
  expect_equal(
    kappa(
      data.frame(sound_speed = 1500, density = 1026),
      data.frame(sound_speed = 1400, density = 1010)
    ),
    expected,
    tolerance = 1e-10
  )

  # Test pois
  # ---- Case: K & E
  expect_equal(pois(K = 1, E = 1), 1 / 3)
  # ---- Case: K & G
  expect_equal(pois(K = 1, G = 1), 1)
  # ---- Case: E & G
  expect_equal(pois(E = 1, G = 1), -0.5)
  # ---- Error: sparse parameters -> K
  expect_error(
    pois(K = 1),
    paste0(
      "At least two elasticity moduli values are required to calculate ",
      "Poisson's ratio."
    )
  )
  # ---- Error: sparse parameters -> E
  expect_error(
    pois(E = 1),
    paste0(
      "At least two elasticity moduli values are required to calculate ",
      "Poisson's ratio."
    )
  )
  # ---- Error: sparse parameters -> G
  expect_error(
    pois(G = 1),
    paste0(
      "At least two elasticity moduli values are required to calculate ",
      "Poisson's ratio."
    )
  )
  # ---- Overloaded parameters: K & E & G
  expect_equal(pois(K = 1, E = 1, G = 1), -0.5)

  # Test bulk
  # ---- Case: E & G
  expect_equal(bulk(E = 1, G = 1), 1 / 6)
  # ---- Case: E & nu
  expect_equal(bulk(E = 1, nu = 1), -1 / 3)
  # ---- Case: G & nu
  expect_equal(bulk(G = 1, nu = 1), -1 - 1 / 3)
  # ---- Error: sparse parameters -> E
  expect_error(
    bulk(E = 1),
    paste0(
      "At least two elasticity moduli values are required to calculate ",
      "the bulk modulus."
    )
  )
  # ---- Error: sparse parameters -> G
  expect_error(
    bulk(G = 1),
    paste0(
      "At least two elasticity moduli values are required to calculate ",
      "the bulk modulus."
    )
  )
  # ---- Error: sparse parameters -> nu
  expect_error(
    bulk(nu = 1),
    paste0(
      "At least two elasticity moduli values are required to calculate ",
      "the bulk modulus."
    )
  )
  # ---- Overloaded parameters: E & G & nu
  expect_equal(bulk(E = 1, G = 1, nu = 1), 1 / 6)

  # Test young
  # ---- Case: K & G
  expect_equal(young(K = 1, G = 1), 2.25)
  # ---- Case: K & nu
  expect_equal(young(K = 1, nu = 1), -3)
  # ---- Case: G & nu
  expect_equal(young(G = 1, nu = 1), 4)
  # ---- Error: sparse parameters -> E
  expect_error(
    young(K = 1),
    paste0(
      "At least two elasticity moduli values are required to calculate ",
      "Young's modulus."
    )
  )
  # ---- Error: sparse parameters -> G
  expect_error(
    young(G = 1),
    paste0(
      "At least two elasticity moduli values are required to calculate ",
      "Young's modulus."
    )
  )
  # ---- Error: sparse parameters -> nu
  expect_error(
    young(nu = 1),
    paste0(
      "At least two elasticity moduli values are required to calculate ",
      "Young's modulus."
    )
  )
  # ---- Overloaded parameters: K & G & nu
  expect_equal(young(K = 1, G = 1, nu = 1), 2.25)

  # Test shear
  # ---- Case: K & E
  expect_equal(shear(K = 1, E = 1), 0.375)
  # ---- Case: K & nu
  expect_equal(shear(K = 1, nu = 1), -0.75)
  # ---- Case: E & nu
  expect_equal(shear(E = 1, nu = 1), 0.25)
  # ---- Error: sparse parameters -> E
  expect_error(
    shear(K = 1),
    paste0(
      "At least two elasticity moduli values are required to calculate ",
      "the shear modulus."
    )
  )
  # ---- Error: sparse parameters -> E
  expect_error(
    shear(E = 1),
    paste0(
      "At least two elasticity moduli values are required to calculate ",
      "the shear modulus."
    )
  )
  # ---- Error: sparse parameters -> nu
  expect_error(
    shear(nu = 1),
    paste0(
      "At least two elasticity moduli values are required to calculate ",
      "the shear modulus."
    )
  )
  # ---- Overloaded parameters: K & E & nu
  expect_equal(shear(K = 1, E = 1, nu = 1), 0.375)

  # Test lame
  # ---- Case: K & E
  expect_equal(lame(K = 1, E = 1), 0.75)
  # ---- Case: K & G
  expect_equal(lame(K = 1, G = 1), 1 / 3)
  # ---- Case: K & nu
  expect_equal(lame(K = 1, nu = 1), 1.5)
  # ---- Case: E & G
  expect_equal(lame(E = 1, G = 1), -0.5)
  # ---- Case: E & nu
  expect_equal(lame(E = 1, nu = 1), -0.5)
  # ---- Case: G & nu
  expect_equal(lame(G = 1, nu = 1), -2)
  # ---- Error: sparse parameters -> K
  expect_error(
    lame(K = 1),
    paste0(
      "At least two elasticity moduli values are required to calculate ",
      "the Lamé parameter."
    )
  )
  # ---- Error: sparse parameters -> E
  expect_error(
    lame(E = 1),
    paste0(
      "At least two elasticity moduli values are required to calculate ",
      "the Lamé parameter."
    )
  )
  # ---- Error: sparse parameters -> G
  expect_error(
    lame(G = 1),
    paste0(
      "At least two elasticity moduli values are required to calculate ",
      "the Lamé parameter."
    )
  )
  # ---- Error: sparse parameters -> G
  expect_error(
    lame(nu = 1),
    paste0(
      "At least two elasticity moduli values are required to calculate ",
      "the Lamé parameter."
    )
  )
  # ---- Overloaded parameters: K & E & G & nu
  expect_equal(lame(K = 1, E = 1, G = 1, nu = 1), 1 / 3)

  # Test rho
  expected <- -0.015841584158
  expect_equal(
    rho(
      data.frame(sound_speed = 1500, density = 1026),
      data.frame(sound_speed = 1400, density = 1010)
    ),
    expected,
    tolerance = 1e-10
  )
})

test_that("Curvature edge-case", {
  data(krill)

  # "Brake" the krill
  w <- capture_warnings(brake(krill, radius_curvature = 20e-3))
  expect_match(w, ".*Arc angle per segment", all = FALSE)
  expect_match(w, ".*One or more body segments", all = FALSE)
})

test_that("Rotation functions for the KRM work correctly", {
  data(sardine, package = "acousticTS")

  # Extract
  body <- acousticTS::extract(sardine, "body")
  sb <- acousticTS::extract(sardine, "bladder")

  # Compute the along-matrix sums
  body_sums <- acousticTS::along_sum(body$rpos, 5)
  bladder_sums <- acousticTS::along_sum(sb$rpos, 5)

  # Test body rotation
  rotated_body <- acousticTS::body_rotation(
    body_sums, body$rpos, pi / 4, 10
  )

  expect_true(is.list(rotated_body))
  expect_equal(names(rotated_body), c("vbU", "vbL", "delta_u"))
  expect_true(is.matrix(rotated_body$vbL))
  expect_true(is.matrix(rotated_body$vbU))
  expect_true(is.vector(rotated_body$delta_u))
  expect_type(rotated_body$vbL, "double")
  expect_type(rotated_body$vbU, "double")
  expect_type(rotated_body$delta_u, "double")
  expect_equal(
    dim(rotated_body$vbU),
    c(10, 4)
  )
  expect_equal(
    dim(rotated_body$vbL),
    c(10, 4)
  )
  expect_equal(
    length(rotated_body$delta_u),
    length(body$rpos[1, ]) - 1
  )

  # Test bladder rotation
  rotated_bladder <- acousticTS::bladder_rotation(
    bladder_sums, sb$rpos, pi / 4, 10
  )

  expect_true(is.list(rotated_bladder))
  expect_equal(names(rotated_bladder), c("v", "delta_u"))
  expect_true(is.matrix(rotated_bladder$v))
  expect_true(is.vector(rotated_bladder$delta_u))
  expect_type(rotated_bladder$v, "double")
  expect_type(rotated_bladder$delta_u, "double")
  expect_equal(
    dim(rotated_bladder$v),
    c(10, 4)
  )
  expect_equal(
    length(rotated_bladder$delta_u),
    length(sb$rpos[1, ]) - 1
  )
})

test_that("SDWBA resampling works as intended", {
  data(krill, package = "acousticTS")
  data(sardine, package = "acousticTS")

  # Test object type-error
  expect_error(
    sdwba_resample(sardine, 100),
    "Object must be of class FLS"
  )

  # Test base case
  krill_resampled <- sdwba_resample(krill, 50)

  # Extract body
  body <- acousticTS::extract(krill_resampled, "body")
  shape <- acousticTS::extract(krill_resampled, "shape_parameters")

  # Test dimensions
  expect_equal(dim(body$rpos), c(5, 51))
  expect_equal(length(body$radius), 51)
  expect_equal(shape$n_segments, 50)
})
