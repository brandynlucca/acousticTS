library(testthat)
library(acousticTS)

test_check("acousticTS", reporter = testthat::SummaryReporter$new())
