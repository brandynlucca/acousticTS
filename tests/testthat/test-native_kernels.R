test_that("native numerical kernels pass their independent assertions", {
  skip_on_os("solaris")
  output <- capture.output(
    passed <- acousticTS:::native_kernel_tests_cpp()
  )
  expect_true(passed, info = paste(output, collapse = "\n"))
})
