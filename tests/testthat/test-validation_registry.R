library(acousticTS)

test_that("validation registry derives family statuses consistently", {
  expect_equal(
    acousticTS:::.validation_family_status("hpa"),
    c("benchmarked", "partially_validated")
  )
  expect_equal(
    acousticTS:::.validation_family_status("dwba"),
    c("benchmarked", "partially_validated")
  )
  expect_equal(
    acousticTS:::.validation_family_status("sdwba"),
    c("benchmarked", "partially_validated")
  )
  expect_equal(
    acousticTS:::.validation_family_status("ecms"),
    c("partially_validated", "experimental")
  )
  expect_equal(
    acousticTS:::.validation_family_status("vesm"),
    c("benchmarked", "partially_validated", "experimental")
  )
  expect_equal(
    acousticTS:::.validation_family_status("bcms"),
    c("partially_validated", "experimental")
  )
  expect_equal(
    acousticTS:::.validation_family_status("bbfm"),
    c("partially_validated", "experimental")
  )
  expect_equal(
    acousticTS:::.validation_family_status("tmm"),
    c("benchmarked", "partially_validated", "experimental")
  )
  expect_equal(
    acousticTS:::.validation_family_status("krm"),
    c("benchmarked", "partially_validated")
  )
  expect_equal(
    acousticTS:::.validation_family_status("fcms"),
    c("benchmarked", "partially_validated")
  )
  expect_equal(
    acousticTS:::.validation_family_status("psms"),
    c("benchmarked", "partially_validated")
  )
})

test_that("validation registry tables cover every documented family", {
  meta <- acousticTS:::.validation_family_registry()
  status_table <- acousticTS:::.validation_family_status_table()
  evidence_table <- acousticTS:::.validation_evidence_table(
    type = c(
      "benchmark", "validated", "partially_validated", "experimental"
    )
  )

  expect_equal(nrow(status_table), nrow(meta))
  expect_setequal(status_table$Family, meta$display)
  expect_true(all(c("Family", "Status") %in% names(status_table)))
  expect_true(all(c("Family", "Evidence type", "Source", "Scope", "Summary") %in% names(evidence_table)))
})

test_that("status-implied benchmark and validation articles exist and are listed", {
  meta <- acousticTS:::.validation_family_registry()
  root_candidates <- c(".", "../..")
  root <- root_candidates[file.exists(file.path(root_candidates, "_pkgdown.yml"))][[1]]
  pkgdown_yml <- readLines(file.path(root, "_pkgdown.yml"), warn = FALSE)

  for (i in seq_len(nrow(meta))) {
    family <- meta$family[[i]]
    stem <- tolower(meta$display[[i]])
    if (identical(family, "calibration")) {
      stem <- "calibration"
    }

    if (isTRUE(meta$benchmarked[[i]])) {
      article <- file.path(root, "vignettes", family, paste0(stem, "-benchmarks.Rmd"))
      yml_entry <- paste0("      - ", family, "/", stem, "-benchmarks")
      expect_true(file.exists(article), info = article)
      expect_true(any(pkgdown_yml == yml_entry), info = yml_entry)
    }

    if (meta$validation_status[[i]] %in% c("validated", "partially_validated")) {
      article <- file.path(root, "vignettes", family, paste0(stem, "-validation.Rmd"))
      yml_entry <- paste0("      - ", family, "/", stem, "-validation")
      expect_true(file.exists(article), info = article)
      expect_true(any(pkgdown_yml == yml_entry), info = yml_entry)
    }
  }
})

test_that("family status metadata stays aligned with evidence rows", {
  meta <- acousticTS:::.validation_family_registry()

  for (family in meta$family) {
    meta_row <- acousticTS:::.validation_family_meta(family)
    evidence <- acousticTS:::.validation_family_evidence(family)
    validation_rows <- evidence$evidence_type %in% c(
      "validated", "partially_validated"
    )

    expect_identical(
      isTRUE(meta_row$benchmarked[[1]]),
      any(evidence$evidence_type == "benchmark")
    )

    if (identical(meta_row$validation_status[[1]], "validated")) {
      expect_true(any(evidence$evidence_type == "validated"))
      expect_false(any(evidence$evidence_type == "partially_validated"))
    } else if (identical(meta_row$validation_status[[1]], "partially_validated")) {
      expect_true(any(evidence$evidence_type == "partially_validated"))
    } else {
      expect_false(any(validation_rows))
    }

    expect_identical(
      isTRUE(meta_row$experimental[[1]]),
      any(evidence$evidence_type == "experimental")
    )
  }
})

test_that("software-distribution comparisons are benchmark evidence", {
  evidence <- acousticTS:::.validation_evidence_registry()
  benchmark_only_rows <- grepl(
    "echoSMs|KRMr|NOAA|ZooScatR|Echopop|software|McGehee|Published|formula",
    evidence$source
  )

  expect_true(all(evidence$evidence_type[benchmark_only_rows] == "benchmark"))
  expect_false(any(
    evidence$evidence_type[benchmark_only_rows] %in%
      c("validated", "partially_validated")
  ))
})

test_that("validation registry helper branches reject invalid inputs cleanly", {
  expect_error(
    acousticTS:::.validation_status_badge("mystery"),
    "Unknown model status"
  )
  expect_equal(
    acousticTS:::.validation_normalize_family("SoEmS"),
    "calibration"
  )
  expect_error(
    acousticTS:::.validation_normalize_family(""),
    "`family` must be a single non-empty string."
  )
  expect_error(
    acousticTS:::.validation_normalize_family("unknown"),
    "Unknown model family"
  )

  meta <- acousticTS:::.validation_family_meta("tmm")
  expect_equal(meta$validation_status[[1]], "partially_validated")
  expect_true(any(grepl("validated", acousticTS:::.validation_family_validation("tmm"))))
  policy <- as.character(acousticTS:::.validation_status_policy())
  expect_match(policy, "<ul>", fixed = TRUE)
  expect_match(policy, "Benchmarked", fixed = TRUE)
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
  expect_false(grepl(
    "Allowed `evidence_type` values are exactly:",
    policy,
    fixed = TRUE
  ))
  expect_match(acousticTS:::.validation_model_library(), "Composite and emerging families")
  expect_equal(acousticTS:::.validation_family_meta("calibration")$display[[1]], "SOEMS")
  expect_match(
    acousticTS:::.validation_family_validation("essms"),
    "spherical-shell validation coverage",
    fixed = TRUE
  )
})

test_that("validation badge tooltips come from registry summaries", {
  tmm_tooltip <- acousticTS:::.validation_status_tooltip(
    "tmm",
    "partially_validated"
  )
  expect_match(tmm_tooltip, "TMM is partially validated", fixed = TRUE)

  essms_tooltip <- acousticTS:::.validation_status_tooltip(
    "essms",
    "unvalidated"
  )
  expect_identical(
    essms_tooltip,
    "The package does not yet claim independent validation across the current public scope."
  )

  badge <- acousticTS:::.validation_status_badge(
    "validated",
    family = "sphms",
    tooltip = TRUE
  )
  expect_match(badge, 'data-tooltip="', fixed = TRUE)
  expect_match(badge, 'title="', fixed = TRUE)
  expect_match(badge, 'tabindex="0"', fixed = TRUE)
  expect_match(
    badge,
    "Validation uses monostatic spherical BEM far-field checks",
    fixed = TRUE
  )

  library_page <- as.character(acousticTS:::.validation_model_library())
  expect_match(library_page, 'data-tooltip="', fixed = TRUE)
  expect_match(
    library_page,
    "TMM is partially validated",
    fixed = TRUE
  )
})
