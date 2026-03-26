library(acousticTS)

test_that("validation registry derives family statuses consistently", {
  expect_equal(
    acousticTS:::.validation_family_status("hpa"),
    c("benchmarked", "validated")
  )
  expect_equal(
    acousticTS:::.validation_family_status("bcms"),
    c("experimental", "unvalidated")
  )
  expect_equal(
    acousticTS:::.validation_family_status("bbfm"),
    c("experimental", "unvalidated")
  )
  expect_equal(
    acousticTS:::.validation_family_status("tmm"),
    c("benchmarked", "validated", "experimental")
  )
})

test_that("validation registry tables cover every documented family", {
  meta <- acousticTS:::.validation_family_registry()
  status_table <- acousticTS:::.validation_family_status_table()
  evidence_table <- acousticTS:::.validation_evidence_table(
    type = c("benchmark", "external"),
    include_notes = TRUE
  )

  expect_equal(nrow(status_table), nrow(meta))
  expect_setequal(status_table$Family, meta$display)
  expect_true(all(c("Family", "Status") %in% names(status_table)))
  expect_true(all(c("Family", "Evidence type", "Source", "Scope", "Summary") %in% names(evidence_table)))
})

test_that("unvalidated status is only assigned when no evidence tags exist", {
  meta <- acousticTS:::.validation_family_registry()

  for (family in meta$family) {
    status <- acousticTS:::.validation_family_status(family)
    evidence <- acousticTS:::.validation_family_evidence(family)
    has_evidence_tag <- any(evidence$evidence_type %in% c("benchmark", "external"))

    if ("unvalidated" %in% status) {
      expect_false(has_evidence_tag)
    } else {
      expect_true(has_evidence_tag || "experimental" %in% status || "benchmarked" %in% status || "validated" %in% status)
    }
  }
})
