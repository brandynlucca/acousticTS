test_that("Simple acoustic utility functions work as intended", {
  library(acousticTS)
  
  # Test sigma_bs
  expected <- 2e-20
  expect_equal(sigma_bs(1e-10 * 1e-10*1i), expected)

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
  expected  =  0.1661446757
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
  expect_equal(pois(K = 1, E = 1), 1/3)
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
  expect_equal(bulk(E = 1, G = 1), 1/6)
  # ---- Case: E & nu
  expect_equal(bulk(E = 1, nu = 1), -1/3)
  # ---- Case: G & nu
  expect_equal(bulk(G = 1, nu = 1), -1 - 1/3)
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
  expect_equal(bulk(E = 1, G = 1, nu = 1), 1/6)
  
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
  expect_equal(lame(K = 1, G = 1), 1/3)
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
  expect_equal(lame(K = 1, E = 1, G = 1, nu = 1), 1/3)
  
  # Test rho
  expected <- -0.015841584158
  expect_equal(
    rho(
      data.frame(sound_speed=1500, density=1026), 
      data.frame(sound_speed=1400, density=1010)
    ),
    expected,
    tolerance=1e-10
  )
  
})
