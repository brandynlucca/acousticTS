################################################################################
# PKGDOWN / VIGNETTE PRESENTATION HELPERS
################################################################################

#' @noRd
.model_status_badge <- function(status) {
  status_map <- list(
    benchmarked = list(css = "established", label = "Benchmarked"),
    validated = list(css = "validated", label = "Validated"),
    experimental = list(css = "experimental", label = "Experimental"),
    unvalidated = list(css = "unvalidated", label = "Unvalidated")
  )

  spec <- status_map[[status]]
  if (is.null(spec)) {
    stop("Unknown model status: ", status, call. = FALSE)
  }

  paste0(
    '<span class="model-tag ', spec$css, '">',
    spec$label,
    "</span>"
  )
}

#' @noRd
.model_family_header <- function(status = character(), pages = c()) {
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
