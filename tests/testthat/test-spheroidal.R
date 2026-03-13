test_that("Spheroidal wave functions work correctly", {
  
  # Angular wave function, Smn
  expect_equal(
    Smn(2, 3, 1, 0.5)$value,
    5.650368053851631
  )
  expect_equal(
    Smn(2, 3, 1, 0.5)$derivative,
    3.454326321112444
  )
  expect_equal(
    Smn(2, 3, 1, 0.0)$value,
    0
  )
  expect_equal(
    Smn(2, 3, 1, 0.0)$derivative,
    15.2772229786631266542735101
  )
  expect_equal(
    Smn(0, 3, 1, 0.5)$value,
    -0.4302211279618986
  )
  expect_equal(
    Smn(0, 3, 1, 0.5)$derivative,
    0.4316245260227451
  )
})