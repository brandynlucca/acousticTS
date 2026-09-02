if (requireNamespace("testthat", quietly = TRUE)) {
  testthat::test_check(
    "acousticTS",
    reporter = testthat::SummaryReporter$new()
  )
}
