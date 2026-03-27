library(acousticTS)

test_that("pkgdown helper markup reflects validation metadata", {
  badge <- acousticTS:::.model_status_badge("benchmarked")
  expect_match(badge, "Benchmarked")
  expect_match(badge, "model-tag")

  header <- as.character(
    acousticTS:::.model_family_header(
      family = "hpa",
      pages = c(
        Overview = "index.html",
        Theory = "theory.html",
        Empty = ""
      )
    )
  )
  expect_match(header, "Validated")
  expect_match(header, "Overview")
  expect_false(grepl("Empty", header, fixed = TRUE))

  overview <- as.character(
    acousticTS:::.model_family_overview(
      summary = "Summary text",
      family = "hpa",
      core_idea = "Core idea text",
      best_for = c("screening", "compact targets"),
      supports = "sphere",
      assumptions = "high ka",
      reading = c("Overview: index.html", "Theory: theory.html")
    )
  )
  expect_match(overview, "Summary text")
  expect_match(overview, "## Core idea")
  expect_match(overview, "## Validation status")
  expect_match(overview, "sphere")
})

test_that("validation helper pages render expected sections", {
  policy <- as.character(acousticTS:::.validation_status_policy())
  library_page <- as.character(acousticTS:::.validation_model_library())

  expect_match(policy, "Status tags are derived")
  expect_match(policy, "Unvalidated")
  expect_match(library_page, "Modal-series families")
  expect_match(library_page, "Approximation and ray-based families")
  expect_match(library_page, "Composite and emerging families")
})
