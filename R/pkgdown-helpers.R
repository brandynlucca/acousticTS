################################################################################
# PKGDOWN / VIGNETTE PRESENTATION HELPERS
################################################################################
# Generate modal status badge for pkgdown site
#' @noRd
.model_status_badge <- function(status) {
  .validation_status_badge(status)
}

# Model subdirectory linking for pkgdown site
#' @noRd
.model_family_page_order <- function() {
  c("Overview", "Theory", "Implementation")
}

# Detect current page to hide from list to avoid a duplicated listing
#' @noRd
.model_family_current_page <- function() {
  if (!requireNamespace("knitr", quietly = TRUE)) {
    return(NULL)
  }

  input <- tryCatch(knitr::current_input(), error = function(...) "")
  if (length(input) != 1L || is.na(input) || !nzchar(input)) {
    return(NULL)
  }

  input <- tolower(basename(input))

  if (grepl("^index\\.[rq]md$", input)) {
    return("Overview")
  }
  if (grepl("-theory\\.[rq]md$", input)) {
    return("Theory")
  }
  if (grepl("-implementation\\.[rq]md$", input)) {
    return("Implementation")
  }

  NULL
}

# Order the landing pages
#' @noRd
.model_family_pages_ordered <- function(pages) {
  pages <- pages[!is.na(pages) & nzchar(pages)]
  if (length(pages) == 0L) {
    return(pages)
  }

  page_names <- names(pages)
  if (is.null(page_names)) {
    page_names <- rep("", length(pages))
    names(pages) <- page_names
  }

  canonical <- .model_family_page_order()
  is_canonical <- nzchar(page_names) & page_names %in% canonical
  canonical_pages <- pages[is_canonical]
  canonical_pages <- canonical_pages[
    order(match(names(canonical_pages), canonical))
  ]

  other_pages <- pages[!is_canonical]
  c(canonical_pages, other_pages)
}

# Format the badges
#' @noRd
.model_family_nav_badge <- function(label, href) {
  paste0(
    '<a class="model-nav-badge" href="', href, '">',
    label,
    "</a>"
  )
}

# Format and organize model family header landing pages
#' @noRd
.model_family_header <- function(family = NULL,
                                 status = character(),
                                 pages = c(),
                                 current_page = NULL) {
  if (!is.null(family) && length(status) == 0L) {
    status <- .validation_family_status(family)
  }

  status_markup <- if (length(status) > 0) {
    paste0(
      '<div class="model-status-line">',
      paste(vapply(status, .model_status_badge, character(1)), collapse = " "),
      "</div>"
    )
  } else {
    ""
  }

  pages <- .model_family_pages_ordered(pages)
  if (is.null(current_page)) {
    current_page <- .model_family_current_page()
  }
  page_labels <- names(pages)
  if (is.null(page_labels)) {
    page_labels <- rep("", length(pages))
  }
  blank_labels <- is.na(page_labels) | !nzchar(page_labels)
  page_labels[blank_labels] <- unname(pages)[blank_labels]

  if (!is.null(current_page) && length(current_page) == 1L &&
    !is.na(current_page) &&
    nzchar(current_page)) {
    keep <- page_labels != current_page
    pages <- pages[keep]
    page_labels <- page_labels[keep]
  }

  nav_markup <- if (length(pages) > 0) {
    paste0(
      '<div class="model-nav"><div class="model-nav-links">',
      paste(
        vapply(seq_along(pages), function(i) {
          .model_family_nav_badge(page_labels[[i]], unname(pages)[[i]])
        }, character(1)),
        collapse = " "
      ),
      "</div></div>"
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

# Model family overiew landing page
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
      paste0("## Main assumptions\n\n", paste0("- ", assumptions,
        collapse = "\n"
      ))
    )
  }

  if (length(validation) > 0) {
    sections <- c(
      sections,
      paste0("## Validation status\n\n", paste0("- ", validation,
        collapse = "\n"
      ))
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

# Utility function for handling HTML formatting
#' @noRd
.html_escape <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("\"", "&quot;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  gsub(">", "&gt;", x, fixed = TRUE)
}

# Static, non-interactive vignette figure
#' @noRd
.vignette_static_figure <- function(src, alt) {
  out <- paste0(
    '<img class="vignette-figure" src="', src,
    '" alt="', .html_escape(alt), '" />'
  )

  if (requireNamespace("knitr", quietly = TRUE)) {
    knitr::asis_output(paste0(out, "\n"))
  } else {
    out
  }
}

# Interactive vignette figure configurations that routes to specified API/page
#' @noRd
.vignette_clickable_figure_specs <- function() {
  list(
    getting_started_workflow = list(
      image = "getting-started-workflow.png",
      alt = "Getting started workflow",
      areas = list(
        list(
          href = "../building-shapes/building-shapes.html",
          label = "Build a shape",
          left = 1.81, top = 26.82, width = 21.12, height = 32.11
        ),
        list(
          href = "../building-scatterers/building-scatterers.html",
          label = "Wrap as scatterer",
          left = 34.68, top = 27.24, width = 24.97, height = 32.11
        ),
        list(
          href = "../running-models/running-models.html",
          label = "Run model(s)",
          left = 71.42, top = 27.43, width = 26.77, height = 32.11
        )
      )
    ),
    building_scatterers_map = list(
      image = "building-scatterers-map.png",
      alt = "Updated scatterer-generation and class-hierarchy map",
      areas = list(
        list(
          href = "../../reference/Scatterer-class.html",
          label = "Scatterer class",
          left = 40.94, top = 34.76, width = 20.38, height = 11.23
        ),
        list(
          href = "../../reference/FLS-class.html",
          label = "FLS class",
          left = 5.16, top = 59.45, width = 19.14, height = 10.43
        ),
        list(
          href = "../../reference/GAS-class.html",
          label = "GAS class",
          left = 5.16, top = 71.63, width = 19.16, height = 10.41
        ),
        list(
          href = "../../reference/CAL-class.html",
          label = "CAL class",
          left = 28.58, top = 72.45, width = 19.16, height = 10.41
        ),
        list(
          href = "../../reference/ELA-class.html",
          label = "ELA class",
          left = 5.16, top = 83.83, width = 19.16, height = 10.41
        ),
        list(
          href = "../../reference/CSC-class.html",
          label = "CSC class",
          left = 53.16, top = 65.53, width = 19.22, height = 10.32
        ),
        list(
          href = "../../reference/SBF-class.html",
          label = "SBF class",
          left = 76.79, top = 59.51, width = 19.22, height = 10.32
        ),
        list(
          href = "../../reference/BBF-class.html",
          label = "BBF class",
          left = 76.79, top = 71.67, width = 19.22, height = 10.32
        ),
        list(
          href = "../../reference/ESS-class.html",
          label = "ESS class",
          left = 76.79, top = 83.88, width = 19.22, height = 10.32
        )
      )
    ),
    package_concepts_architecture = list(
      image = "package-concepts-architecture.png",
      alt = "Conceptual layers",
      areas = list(
        list(
          href = "../building-shapes/building-shapes.html",
          label = "Shape concepts",
          left = 0.53, top = 26.49, width = 38.36, height = 43.69
        ),
        list(
          href = "../building-scatterers/building-scatterers.html",
          label = "Scatterer concepts",
          left = 46.03, top = 16.80, width = 23.76, height = 61.72
        ),
        list(
          href = "../running-models/running-models.html",
          label = "Model execution",
          left = 76.70, top = 35.98, width = 22.73, height = 24.41
        )
      )
    ),
    shape_manipulation_schematic = list(
      image = "shape-manipulation-schematic.png",
      alt = "Shape manipulation",
      areas = list(
        list(
          href = "../../reference/brake.html",
          label = "brake()",
          left = 22.29, top = 29.03, width = 15.00, height = 21.34
        ),
        list(
          href = "../../reference/reforge.html",
          label = "reforge()",
          left = 22.29, top = 63.23, width = 15.00, height = 21.34
        )
      )
    ),
    model_selection_flowchart = list(
      image = "model-selection-flowchart.png",
      alt = "Revised general model-selection flowchart",
      areas = list(
        list(
          href = "../krm/index.html",
          label = "KRM overview",
          left = 75.75, top = 24.09, width = 5.47, height = 2.86
        ),
        list(
          href = "../bbfm/index.html",
          label = "BBFM overview",
          left = 55.68, top = 36.55, width = 6.08, height = 2.86
        ),
        list(
          href = "../hpa/index.html",
          label = "HPA overview",
          left = 53.59, top = 56.20, width = 5.34, height = 2.86
        ),
        list(
          href = "../sphms/index.html",
          label = "SPHMS overview",
          left = 77.31, top = 63.58, width = 6.58, height = 2.59
        ),
        list(
          href = "../calibration/index.html",
          label = "Calibration overview",
          left = 77.68, top = 66.65, width = 5.85, height = 2.59
        ),
        list(
          href = "../fcms/index.html",
          label = "FCMS overview",
          left = 90.04, top = 66.00, width = 5.64, height = 2.59
        ),
        list(
          href = "../essms/index.html",
          label = "ESSMS overview",
          left = 77.86, top = 69.72, width = 5.49, height = 2.59
        ),
        list(
          href = "../trcm/index.html",
          label = "TRCM overview",
          left = 89.70, top = 69.09, width = 6.32, height = 2.56
        ),
        list(
          href = "../bcms/index.html",
          label = "BCMS overview",
          left = 89.70, top = 72.15, width = 6.32, height = 2.56
        ),
        list(
          href = "../ecms/index.html",
          label = "ECMS overview",
          left = 89.70, top = 75.21, width = 6.32, height = 2.56
        ),
        list(
          href = "../vesm/index.html",
          label = "VESM overview",
          left = 77.54, top = 72.79, width = 6.12, height = 2.59
        ),
        list(
          href = "../sdwba/index.html",
          label = "SDWBA overview",
          left = 47.03, top = 73.83, width = 7.38, height = 2.83
        ),
        list(
          href = "../dwba/index.html",
          label = "DWBA overview",
          left = 7.12, top = 75.46, width = 6.79, height = 2.84
        ),
        list(
          href = "../psms/index.html",
          label = "PSMS overview",
          left = 69.80, top = 78.14, width = 5.64, height = 2.59
        ),
        list(
          href = "../pcdwba/index.html",
          label = "PCDWBA overview",
          left = 47.03, top = 77.08, width = 7.38, height = 2.83
        ),
        list(
          href = "../tmm/index.html",
          label = "TMM overview",
          left = 90.19, top = 78.17, width = 5.34, height = 2.59
        ),
        list(
          href = "../tmm/index.html",
          label = "TMM overview",
          left = 59.98, top = 78.13, width = 5.34, height = 2.59
        ),
        list(
          href = "../tmm/index.html",
          label = "TMM overview",
          left = 77.93, top = 75.86, width = 5.34, height = 2.59
        ),
        list(
          href = "../tmm/index.html",
          label = "TMM overview",
          left = 69.98, top = 81.22, width = 5.34, height = 2.59
        )
      )
    )
  )
}

# Interactive vignette figure that routes to specified API/page
#' @noRd
.vignette_clickable_figure <- function(id) {
  specs <- .vignette_clickable_figure_specs()
  spec <- specs[[id]]

  if (is.null(spec)) {
    stop("Unknown clickable figure id: ", id, call. = FALSE)
  }

  links <- vapply(spec$areas, function(area) {
    paste0(
      '<a class="clickable-figure__link" href="', area$href,
      '" aria-label="', .html_escape(area$label),
      '" title="', .html_escape(area$label),
      '" style="left:', sprintf("%.2f%%", area$left),
      ";top:", sprintf("%.2f%%", area$top),
      ";width:", sprintf("%.2f%%", area$width),
      ";height:", sprintf("%.2f%%", area$height),
      ';"></a>'
    )
  }, character(1))

  out <- paste0(
    '<div class="clickable-figure">',
    '<img class="clickable-figure__image" src="', spec$image,
    '" alt="', .html_escape(spec$alt), '" />',
    paste(links, collapse = ""),
    "</div>"
  )

  if (requireNamespace("knitr", quietly = TRUE)) {
    knitr::asis_output(paste0(out, "\n"))
  } else {
    out
  }
}
