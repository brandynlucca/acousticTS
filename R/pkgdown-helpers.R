################################################################################
# PKGDOWN / VIGNETTE PRESENTATION HELPERS
################################################################################

#' @noRd
.model_status_badge <- function(status) {
  .validation_status_badge(status)
}

#' @noRd
.model_family_header <- function(family = NULL,
                                 status = character(),
                                 pages = c()) {
  if (!is.null(family) && length(status) == 0L) {
    status <- .validation_family_status(family)
  }

  status_markup <- if (length(status) > 0) {
    paste(vapply(status, .model_status_badge, character(1)), collapse = " ")
  } else {
    ""
  }

  pages <- pages[!is.na(pages) & nzchar(pages)]
  nav_markup <- if (length(pages) > 0) {
    paste0(
      "*Model-family pages:* ",
      paste0("[", names(pages), "](", unname(pages), ")", collapse = " ")
    )
  } else {
    ""
  }

  out <- paste(Filter(nzchar, c(status_markup, nav_markup)), collapse = "\n\n")

  if (requireNamespace("knitr", quietly = TRUE)) {
    knitr::asis_output(paste0(out, "\n"))
  } else {
    out
  }
}

#' @noRd
.model_family_overview <- function(summary,
                                   family = NULL,
                                   core_idea = NULL,
                                   best_for = character(),
                                   supports = character(),
                                   assumptions = character(),
                                   validation = character(),
                                   reading = character()) {
  if (!is.null(family) && length(validation) == 0L) {
    validation <- .validation_family_validation(family)
  }

  sections <- c(summary)

  if (!is.null(core_idea) && nzchar(core_idea)) {
    sections <- c(sections, paste0("## Core idea\n\n", core_idea))
  }

  if (length(best_for) > 0) {
    sections <- c(
      sections,
      paste0("## Best for\n\n", paste0("- ", best_for, collapse = "\n"))
    )
  }

  if (length(supports) > 0) {
    sections <- c(
      sections,
      paste0("## Supports\n\n", paste0("- ", supports, collapse = "\n"))
    )
  }

  if (length(assumptions) > 0) {
    sections <- c(
      sections,
      paste0("## Main assumptions\n\n", paste0("- ", assumptions, collapse = "\n"))
    )
  }

  if (length(validation) > 0) {
    sections <- c(
      sections,
      paste0("## Validation status\n\n", paste0("- ", validation, collapse = "\n"))
    )
  }

  if (length(reading) > 0) {
    sections <- c(
      sections,
      paste0("## Family pages\n\n", paste0("- ", reading, collapse = "\n"))
    )
  }

  out <- paste(sections, collapse = "\n\n")

  if (requireNamespace("knitr", quietly = TRUE)) {
    knitr::asis_output(paste0(out, "\n"))
  } else {
    out
  }
}
