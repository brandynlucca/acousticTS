library(acousticTS)

test_that("pkgdown helper markup reflects validation metadata", {
  badge <- acousticTS:::.model_status_badge("benchmarked")
  expect_match(badge, "Benchmarked")
  expect_match(badge, "model-tag")

  partial_badge <- acousticTS:::.model_status_badge("partially_validated")
  expect_match(partial_badge, "Partially validated")

  header <- as.character(
    acousticTS:::.model_family_header(
      family = "hpa",
      pages = c(
        Theory = "theory.html",
        Implementation = "implementation.html",
        Overview = "index.html",
        Empty = ""
      ),
      current_page = "Theory"
    )
  )
  expect_match(header, "Validated")
  expect_true(grepl("href=\"index.html\"", header, fixed = TRUE))
  expect_true(grepl("href=\"implementation.html\"", header, fixed = TRUE))
  expect_false(grepl("href=\"theory.html\"", header, fixed = TRUE))
  expect_false(grepl("Empty", header, fixed = TRUE))
  expect_lt(
    regexpr("href=\"index.html\"", header, fixed = TRUE)[[1]],
    regexpr("href=\"implementation.html\"", header, fixed = TRUE)[[1]]
  )

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

  expect_match(policy, "<ul>", fixed = TRUE)
  expect_match(policy, "Partially validated", fixed = TRUE)
  expect_match(
    policy,
    "These tags are intended to be read in three pieces:",
    fixed = TRUE
  )
  expect_false(grepl(
    ".validation_evidence_registry()",
    policy,
    fixed = TRUE
  ))
  expect_match(library_page, "Modal-series families")
  expect_match(library_page, "Approximation and ray-based families")
  expect_match(library_page, "Composite and emerging families")
})

test_that("vignette figure helpers render png-based markup", {
  static_figure <- as.character(
    acousticTS:::.vignette_static_figure(
      src = "example.png",
      alt = "Example figure"
    )
  )

  expect_match(static_figure, 'src="example.png"', fixed = TRUE)
  expect_match(static_figure, 'alt="Example figure"', fixed = TRUE)
  expect_match(static_figure, 'class="vignette-figure"', fixed = TRUE)

  clickable_figure <- as.character(
    acousticTS:::.vignette_clickable_figure("getting_started_workflow")
  )

  expect_match(
    clickable_figure,
    'getting-started-workflow.png',
    fixed = TRUE
  )
  expect_match(
    clickable_figure,
    'href="../building-shapes/building-shapes.html"',
    fixed = TRUE
  )
  expect_match(
    clickable_figure,
    'aria-label="Build a shape"',
    fixed = TRUE
  )
  expect_match(
    clickable_figure,
    'class="clickable-figure__link"',
    fixed = TRUE
  )

  model_selection <- as.character(
    acousticTS:::.vignette_clickable_figure("model_selection_flowchart")
  )

  expect_match(
    model_selection,
    'href="../tmm/index.html"',
    fixed = TRUE
  )
  expect_false(grepl('href="../ttm/index.html"', model_selection, fixed = TRUE))
})

test_that("pkgdown helpers cover fallback output and ordering branches", {
  testthat::local_mocked_bindings(
    requireNamespace = function(...) FALSE,
    .package = "base"
  )

  expect_null(acousticTS:::.model_family_current_page())
  expect_equal(acousticTS:::.model_family_pages_ordered(character()), character())

  unnamed_pages <- acousticTS:::.model_family_pages_ordered(
    c("implementation.html", "index.html")
  )
  expect_equal(names(unnamed_pages), c("", ""))

  header <- acousticTS:::.model_family_header(
    family = NULL,
    status = character(),
    pages = character(),
    current_page = NULL
  )
  expect_equal(header, "")

  overview <- acousticTS:::.model_family_overview(
    summary = "Summary only",
    family = NULL
  )
  expect_match(overview, "Summary only")

  static_figure <- acousticTS:::.vignette_static_figure("example.png", "Example")
  expect_match(static_figure, 'src="example.png"', fixed = TRUE)

  expect_error(
    acousticTS:::.vignette_clickable_figure("unknown_id"),
    "Unknown clickable figure id"
  )
  clickable <- acousticTS:::.vignette_clickable_figure("getting_started_workflow")
  expect_match(clickable, 'class="clickable-figure"', fixed = TRUE)
})
