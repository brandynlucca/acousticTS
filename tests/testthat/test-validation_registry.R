library(acousticTS)

test_that("validation registry derives family statuses consistently", {
  expect_equal(
    acousticTS:::.validation_family_status("hpa"),
    c("benchmarked", "validated")
  )
  expect_equal(
    acousticTS:::.validation_family_status("bcms"),
    c("unvalidated", "experimental")
  )
  expect_equal(
    acousticTS:::.validation_family_status("bbfm"),
    c("unvalidated", "experimental")
  )
  expect_equal(
    acousticTS:::.validation_family_status("tmm"),
    c("benchmarked", "partially_validated", "experimental")
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
  expect_equal(
    acousticTS:::.validation_family_validation("essms"),
    "The package does not yet claim external validation across the current public scope."
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
    "The package does not yet claim external validation across the current public scope."
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
    "Validated against `KRMr` and `echoSMs`",
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
