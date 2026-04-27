################################################################################
# INTERNAL VALIDATION / BENCHMARK REGISTRY
################################################################################

#' @noRd
.validation_status_badge <- function(status,
                                     family = NULL,
                                     tooltip = FALSE) {
  status_map <- list(
    benchmarked = list(css = "established", label = "Benchmarked"),
    validated = list(css = "validated", label = "Validated"),
    partially_validated = list(
      css = "partially-validated",
      label = "Partially validated"
    ),
    experimental = list(css = "experimental", label = "Experimental"),
    unvalidated = list(css = "unvalidated", label = "Unvalidated")
  )

  spec <- status_map[[status]]
  if (is.null(spec)) {
    stop("Unknown model status: ", status, call. = FALSE)
  }

  tooltip_attr <- ""
  if (isTRUE(tooltip) && !is.null(family)) {
    tooltip_text <- .validation_status_tooltip(family, status)
    if (length(tooltip_text) == 1L && !is.na(tooltip_text) && nzchar(tooltip_text)) {
      escaped_tooltip <- .validation_html_escape(tooltip_text)
      tooltip_attr <- paste0(
        ' title="', escaped_tooltip,
        '" data-tooltip="', escaped_tooltip,
        '" tabindex="0"'
      )
    }
  }

  paste0(
    '<span class="model-tag ', spec$css, '"', tooltip_attr, ">",
    spec$label,
    "</span>"
  )
}

#' @noRd
.validation_html_escape <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("\"", "&quot;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  gsub(">", "&gt;", x, fixed = TRUE)
}

#' @noRd
.validation_normalize_family <- function(family) {
  if (length(family) != 1L || is.na(family) || !nzchar(family)) {
    stop("`family` must be a single non-empty string.", call. = FALSE)
  }

  family <- tolower(family)

  aliases <- c(
    soems = "calibration"
  )

  if (family %in% names(aliases)) {
    family <- aliases[[family]]
  }

  valid <- .validation_family_registry()$family
  if (!family %in% valid) {
    stop("Unknown model family: ", family, call. = FALSE)
  }

  family
}

#' @noRd
.validation_evidence_types <- function() {
  c("benchmark", "validated", "partially_validated", "experimental")
}

#' @noRd
.validation_validation_status_types <- function() {
  c("validated", "partially_validated", "unvalidated")
}

#' @noRd
.validation_family_registry <- function() {
  data.frame(
    family = c(
      "sphms", "fcms", "psms", "calibration", "essms", "bcms", "ecms",
      "dwba", "sdwba", "krm", "hpa", "trcm", "pcdwba",
      "bbfm", "vesm", "tmm"
    ),
    display = c(
      "SPHMS", "FCMS", "PSMS", "SOEMS", "ESSMS", "BCMS", "ECMS",
      "DWBA", "SDWBA", "KRM", "HPA", "TRCM", "PCDWBA",
      "BBFM", "VESM", "TMM"
    ),
    section = c(
      rep("Modal-series families", 7L),
      rep("Approximation and ray-based families", 6L),
      rep("Composite and emerging families", 3L)
    ),
    section_order = c(rep(1L, 7L), rep(2L, 6L), rep(3L, 3L)),
    family_order = c(
      1L, 2L, 3L, 4L, 5L, 6L, 7L,
      1L, 2L, 3L, 4L, 5L, 6L,
      1L, 2L, 3L
    ),
    url = c(
      "../sphms/index.html",
      "../fcms/index.html",
      "../psms/index.html",
      "../calibration/index.html",
      "../essms/index.html",
      "../bcms/index.html",
      "../ecms/index.html",
      "../dwba/index.html",
      "../sdwba/index.html",
      "../krm/index.html",
      "../hpa/index.html",
      "../trcm/index.html",
      "../pcdwba/index.html",
      "../bbfm/index.html",
      "../vesm/index.html",
      "../tmm/index.html"
    ),
    description = c(
      "Spherical modal-series solution for canonical spherical targets.",
      "Finite-cylinder modal-series solution for straight cylindrical targets.",
      paste(
        "Prolate-spheroidal modal-series solution for smooth",
        "elongated canonical bodies."
      ),
      "Solid elastic spherical model used mainly for calibration spheres.",
      "Elastic-shelled spherical family for layered shell targets.",
      paste(
        "Bent-cylinder modal-series family for straight and uniformly bent",
        "cylinders."
      ),
      "Elastic-cylinder modal-series family for fully elastic solid cylinders.",
      "Weak-scattering elongated-body approximation for fluid-like targets.",
      "Stochastic DWBA family for unresolved phase variability.",
      paste(
        "Kirchhoff-ray mode model for segmented fish-like ",
        "body-plus-inclusion targets."
      ),
      "High-pass approximation for compact asymptotic screening.",
      paste(
        "Two-ray cylindrical family for high-frequency locally",
        "cylindrical targets."
      ),
      "Phase-compensated DWBA for bent weakly scattering targets.",
      paste(
        "Composite flesh-plus-backbone family for swimbladder-less ",
        "fish-like targets."
      ),
      paste(
        "Viscous-elastic layered-sphere family for gas-core, shell,",
        "and viscous-layer targets."
      ),
      paste(
        "Single-target transition-matrix family for retained monostatic",
        "and angle-dependent scattering products across supported",
        "canonical shapes."
      )
    ),
    # Registry maintenance notes:
    # - `benchmarked` controls the public `Benchmarked` badge.
    # - `validation_status` must be exactly one of `validated`,
    #   `partially_validated`, or `unvalidated`.
    # - `experimental` controls the separate public `Experimental` badge.
    # - Update these public badge fields first, then keep the supporting rows
    #   in `.validation_evidence_registry()` aligned with them.
    benchmarked = c(
      TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE,
      TRUE, TRUE, TRUE, TRUE, TRUE, FALSE,
      FALSE, FALSE, TRUE
    ),
    validation_status = c(
      "validated", "validated", "validated", "validated",
      "unvalidated", "unvalidated", "unvalidated",
      "validated", "validated", "validated", "validated", "unvalidated",
      "validated", "unvalidated", "validated", "partially_validated"
    ),
    experimental = c(
      FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE,
      FALSE, FALSE, FALSE, FALSE, FALSE, TRUE,
      TRUE, TRUE, TRUE
    ),
    stringsAsFactors = FALSE
  )
}

#' @noRd
.validation_evidence_registry <- function() {
  # Registry maintenance rules:
  # - `evidence_type` must be exactly one of `benchmark`, `validated`,
  #   `partially_validated`, or `experimental`.
  # - `benchmark` means a documented canonical benchmark ladder or stored
  #   benchmark values.
  # - `validated` means the current public scope has an external software or
  #   independent-comparison check.
  # - `partially_validated` means only part of the public scope has those
  #   external checks.
  # - `experimental` means the workflow is public but still provisional.
  # - `unvalidated` is tracked only in `.validation_family_registry()`.
  # - Do not add legacy umbrella tags such as `note` or `external`.
  data.frame(
    family = c(
      "sphms", "sphms",
      "fcms", "fcms",
      "psms", "psms",
      "calibration", "calibration",
      "bcms",
      "ecms",
      "dwba", "dwba",
      "sdwba", "sdwba",
      "krm", "krm",
      "hpa", "hpa",
      "trcm",
      "pcdwba", "pcdwba",
      "bbfm",
      "vesm", "vesm",
      "tmm", "tmm", "tmm", "tmm", "tmm"
    ),
    evidence_type = c(
      "benchmark", "validated",
      "benchmark", "validated",
      "benchmark", "validated",
      "benchmark", "validated",
      "experimental",
      "experimental",
      "benchmark", "validated",
      "benchmark", "validated",
      "benchmark", "validated",
      "benchmark", "validated",
      "benchmark",
      "validated", "experimental",
      "experimental",
      "validated", "experimental",
      "benchmark", "validated", "validated",
      "partially_validated", "experimental"
    ),
    source = c(
      # SPHMS
      "benchmark_ts / Jech et al. (2015)", "KRMr and echoSMs",
      # FCMS
      "benchmark_ts / Jech et al. (2015)", "echoSMs",
      # PSMS
      "benchmark_ts / Jech et al. (2015)", "Prol_Spheroid and echoSMs",
      # SOEMS
      "Published calibration spheres", "echoSMs, SphereTS, NOAA applet",
      # BCMS
      "Internal FCMS-based reference reconstruction",
      # ECMS
      "Independent algebra transcription",
      # DWBA
      "benchmark_ts / Jech et al. (2015)",
      "McGehee et al. (1998) and echoSMs",
      # SDWBA
      "benchmark_ts / Jech et al. (2015)",
      "CCAMLR MATLAB, NOAA applet, echoSMs",
      # KRM
      "benchmark_ts / Jech et al. (2015)", "KRMr, echoSMs, NOAA applet",
      # HPA
      "benchmark_ts / Jech et al. (2015)", "echoSMs",
      # TRCM
      "benchmark_ts / Jech et al. (2015)",
      # PCDWBA
      "benchmark_ts / Jech et al. (2015)", "ZooScatR",
      # BBFM
      "Internal DWBA + ECMS reconstruction",
      # VESM
      "Reference Python VESM workflow",
      "Current layered-sphere workflow surface",
      # TMM
      "SPHMS / PSMS / FCMS benchmark ladder",
      "BEMPP far-field checks",
      "Exact general-angle spheroidal solution",
      "Cylinder retained-angle scope",
      "Current retained-state branch matrix"
    ),
    summary = c(
      # SPHMS
      paste(
        "Benchmarked against the canonical spherical spectra stored",
        "in benchmark_ts."
      ),
      paste(
        "Validated against `KRMr` and `echoSMs` on shared",
        "penetrable-sphere cases."
      ),
      # FCMS
      paste(
        "Benchmarked against the canonical finite-cylinder spectra",
        "stored in benchmark_ts."
      ),
      "Validated against the echoSMs finite-cylinder implementation.",
      # PSMS
      paste(
        "Benchmarked against the canonical prolate-spheroid spectra",
        "stored in benchmark_ts."
      ),
      paste(
        "Validated against the external Prol_Spheroid and echoSMs ",
        "implementations on shared prolate cases."
      ),
      # SOEMS / CALIBRATION
      paste(
        "Benchmarked against published calibration-sphere targets used",
        "throughout the package documentation."
      ),
      paste(
        "Validated against echoSMs, SphereTS, and the NOAA calibration",
        "applet."
      ),
      # BCMS
      paste(
        "BCMS is currently marked experimental because the documented checks",
        "are internal coherence reconstructions rather than an external",
        "benchmark or software-comparison ladder."
      ),
      # ECMS
      paste(
        "ECMS is currently marked experimental because the documented checks",
        "are independent algebra reconstructions rather than an external",
        "benchmark or software-comparison ladder."
      ),
      # DWBA
      paste(
        "Benchmarked against the canonical spectra stored",
        "in benchmark_ts."
      ),
      paste(
        "Validated against the published McGehee et al (1998) and echoSMs ",
        "workflows."
      ),
      #SDWBA
      paste(
        "Benchmarked against the canonical spectra stored",
        "in benchmark_ts."
      ),
      paste(
        "Validated against the CCAMLR, NOAA applet, and echoSMs ",
        "implementations."
      ),
      # KRM
      paste(
        "Benchmarked against the canonical spectra stored",
        "in benchmark_ts."
      ),
      paste(
        "Validated against`KRMr, echoSMs, and the NOAA KRM applet on",
        "bundled fish objects and shared workflows."
      ),
      # HPA
      paste(
        "Benchmarked against the canonical spherical spectra stored",
        "in benchmark_ts."
      ),
      paste(
        "Validated against the spherical echoSMs implementation."
      ),
      # TRCM
      paste(
        "Benchmarked within the package validation workflow against ",
        "the canonical spectra stored in benchmark_ts. Further ",
        "compared to the straight-cylinder and FCMS-derived ",
        "bent-cylinder reference constructions."
      ),
      # PCDWBA
      paste(
        "Validated against source-level ZooScatR and Echopop PCDWBA ",
        "workflows."
      ),
      paste(
        "PCDWBA is currently marked experimental because the public package",
        "workflow is still being tightened even though the current source-",
        "level comparison cases are documented."
      ),
      # BBFM
      paste(
        "BBFM is currently marked experimental because it has documented",
        "internal reconstruction checks but no external benchmark ladder or",
        "independent public implementation comparison."
      ),
      # VESM
      paste(
        "Validated against the reference Python VESM implementation on",
        "the documented layered-sphere case."
      ),
      paste(
        "VESM is currently marked experimental because the documented public",
        "workflow is still limited to the current layered-sphere scope."
      ),
      # TMM
      paste(
        "Benchmarked against `SPHMS`, `PSMS`, and `FCMS` on the",
        "currently supported canonical shape branches."
      ),
      paste(
        "Validated against external BEMPP far-field checks for sphere,",
        "oblate, and prolate pressure-release cases."
      ),
      paste(
        "Retained prolate angular products are also checked against the",
        "exact general-angle spheroidal solution."
      ),
      paste(
        "TMM is partially validated because the sphere, oblate, and prolate",
        "branches have external checks, but retained general-angle cylinder",
        "products remain outside the validated public scope."
      ),
      paste(
        "TMM is currently marked experimental because the retained-state",
        "workflow and branch matrix are still guarded while shape-specific",
        "support continues to be tightened."
      )
    ),
    scope = c(
      "Sphere spectra across rigid, soft, liquid-filled, and gas-filled cases.",
      "Penetrable sphere spectra on shared software definitions.",
      paste(
        "Finite-cylinder spectra across the canonical cylindrical",
        "benchmark grid."
      ),
      "Rigid, soft, liquid-filled, and gas-filled finite-cylinder spectra.",
      "Prolate-spheroid spectra across the canonical benchmark grid.",
      "Liquid-filled and gas-filled prolate-spheroid software comparisons.",
      "Tungsten-carbide and copper calibration spheres.",
      "Shared calibration-sphere material sets and frequency sweeps.",
      "Uniform-curvature cylinder coherence extension of FCMS.",
      paste(
        "Elastic-cylinder component family and near-broadside",
        "canonical cases."
      ),
      "Weakly scattering sphere, prolate spheroid, and cylinder targets.",
      "Bundled krill geometry and published DWBA reference workflows.",
      paste(
        "Weakly scattering sphere, prolate spheroid, and cylinder",
        "stochastic targets."
      ),
      "Bundled krill stochastic workflow comparisons.",
      paste(
        "Canonical isolated targets used for the package KRM",
        "benchmark ladder."
      ),
      "Bundled sardine and cod software-to-software comparisons.",
      paste(
        "Sphere, prolate spheroid, and cylinder asymptotic",
        "benchmark targets."
      ),
      "Spherical HPModel branch and published asymptotic formulas.",
      paste(
        "Straight and bent cylindrical validation cases documented",
        "in the package."
      ),
      paste(
        "Curved weak-scattering reference workflows on shared",
        "bent-body cases."
      ),
      "Current package-facing PCDWBA workflow and argument surface.",
      "Internal composite-component consistency checks only.",
      paste(
        "Documented spherical layered case used by the original",
        "VESM implementation."
      ),
      "Current documented layered-sphere workflow surface.",
      "Sphere, oblate, prolate, and guarded cylinder monostatic branches.",
      "Pressure-release angular slices for sphere, oblate, and prolate cases.",
      "General-angle prolate retained-state validation.",
      "Cylinder retained-angle scope limitation and guardrails.",
      "Current retained-state branch matrix across supported shapes."
    ),
    stringsAsFactors = FALSE
  )
}

#' @noRd
.validation_registry <- function() {
  list(
    families = .validation_family_registry(),
    evidence = .validation_evidence_registry()
  )
}

#' @noRd
.validation_family_meta <- function(family) {
  family <- .validation_normalize_family(family)
  meta <- .validation_family_registry()
  meta[meta$family == family, , drop = FALSE]
}

#' @noRd
.validation_family_evidence <- function(family, type = NULL) {
  family <- .validation_normalize_family(family)
  evidence <- .validation_evidence_registry()
  rows <- evidence[evidence$family == family, , drop = FALSE]

  if (!is.null(type)) {
    rows <- rows[rows$evidence_type %in% type, , drop = FALSE]
  }

  rows
}

#' @noRd
.validation_family_status <- function(family) {
  family <- .validation_normalize_family(family)
  meta <- .validation_family_meta(family)
  validation_status <- meta$validation_status[[1]]
  valid_validation_status <- .validation_validation_status_types()

  if (!validation_status %in% valid_validation_status) {
    stop(
      "Unknown validation status for family '", family, "': ",
      validation_status,
      call. = FALSE
    )
  }

  status <- character()

  if (isTRUE(meta$benchmarked[[1]])) {
    status <- c(status, "benchmarked")
  }

  status <- c(status, validation_status)

  if (isTRUE(meta$experimental[[1]])) {
    status <- c(status, "experimental")
  }

  unique(status)
}

#' @noRd
.validation_family_validation <- function(family) {
  meta <- .validation_family_meta(family)
  evidence <- .validation_family_evidence(
    family,
    type = .validation_evidence_types()
  )

  if (nrow(evidence) == 0) {
    if (identical(meta$validation_status[[1]], "unvalidated")) {
      return(
        "The package does not yet claim external validation across the current public scope."
      )
    }

    character()
  } else {
    evidence$summary
  }
}

#' @noRd
.validation_status_tooltip <- function(family, status) {
  family <- .validation_normalize_family(family)

  status_type <- c(
    benchmarked = "benchmark",
    validated = "validated",
    partially_validated = "partially_validated",
    experimental = "experimental",
    unvalidated = "unvalidated"
  )[[status]]

  if (is.null(status_type)) {
    stop("Unknown model status: ", status, call. = FALSE)
  }

  if (identical(status_type, "unvalidated")) {
    return(
      "The package does not yet claim external validation across the current public scope."
    )
  }

  evidence <- .validation_family_evidence(family, type = status_type)
  if (nrow(evidence) == 0) {
    return("")
  }

  paste(unique(evidence$summary), collapse = " ")
}

#' @noRd
.validation_status_policy <- function() {
  badge <- function(x) .validation_status_badge(x)

  out <- paste(
    c(
      '<div class="validation-status-policy">',
      "<ul>",
      paste0(
        "<li>", badge("benchmarked"),
        " means the family has a documented comparison against a canonical",
        " benchmark ladder or stored benchmark values.</li>"
      ),
      paste0(
        "<li>", badge("validated"),
        " means the currently supported public scope has a documented",
        " external software or independent-comparison check.</li>"
      ),
      paste0(
        "<li>", badge("partially_validated"),
        " means some supported branches are externally checked, but the",
        " full public scope is not yet closed.</li>"
      ),
      paste0(
        "<li>", badge("unvalidated"),
        " means the package does not yet claim external validation across",
        " the current public scope.</li>"
      ),
      paste0(
        "<li>", badge("experimental"),
        " means the family is available to use, but its interface or",
        " supported workflow should still be treated as provisional.</li>"
      ),
      "</ul>",
      "<p>These tags are intended to be read in three pieces:</p>",
      "<ul>",
      paste(
        "<li><code>Benchmarked</code> is independent of the validation badge",
        "and can appear alongside <code>Validated</code>,",
        "<code>Partially validated</code>, or <code>Unvalidated</code>.</li>"
      ),
      paste(
        "<li>The validation badge is always exactly one of",
        "<code>Validated</code>, <code>Partially validated</code>, or",
        "<code>Unvalidated</code>.</li>"
      ),
      paste(
        "<li><code>Experimental</code> is a separate lifecycle tag and can",
        "coexist with any benchmark or validation badge.</li>"
      ),
      "</ul>",
      "</div>"
    ),
    collapse = "\n"
  )

  if (requireNamespace("knitr", quietly = TRUE)) {
    knitr::asis_output(paste0(out, "\n"))
  } else {
    out
  }
}

#' @noRd
.validation_family_status_table <- function() {
  registry <- .validation_registry()
  families <- registry$families[order(
    registry$families$section_order,
    registry$families$family_order
  ), , drop = FALSE]

  status_text <- vapply(families$family, function(family) {
    paste(
      vapply(.validation_family_status(family), function(x) {
        .validation_family_label(x)
      }, character(1)),
      collapse = ", "
    )
  }, character(1))

  data.frame(
    Family = families$display,
    Section = families$section,
    Status = status_text,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

#' @noRd
.validation_evidence_table <- function(
  type = c("benchmark", "validated", "partially_validated", "experimental")
) {
  evidence <- .validation_evidence_registry()
  type <- match.arg(type, several.ok = TRUE)
  evidence <- evidence[evidence$evidence_type %in% type, , drop = FALSE]

  families <- .validation_family_registry()[, c("family", "display"),
    drop = FALSE
  ]
  evidence <- merge(evidence, families,
    by = "family", all.x = TRUE,
    sort = FALSE
  )

  evidence_type_order <- c(
    benchmark = 1L,
    validated = 2L,
    partially_validated = 3L,
    experimental = 4L
  )
  evidence <- evidence[order(
    match(evidence$family, .validation_family_registry()$family),
    evidence_type_order[evidence$evidence_type]
  ), , drop = FALSE]

  data.frame(
    Family = evidence$display,
    `Evidence type` = vapply(
      evidence$evidence_type,
      .validation_family_label, character(1)
    ),
    Source = evidence$source,
    Scope = evidence$scope,
    Summary = evidence$summary,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    row.names = NULL
  )
}

#' @noRd
.validation_model_library <- function(heading_level = 1L) {
  registry <- .validation_registry()
  families <- registry$families[order(
    registry$families$section_order,
    registry$families$family_order
  ), , drop = FALSE]

  sections <- unique(families$section)
  section_prefix <- paste(rep("#", heading_level), collapse = "")
  family_prefix <- paste(rep("#", heading_level + 1L), collapse = "")

  out <- character()

  for (section in sections) {
    out <- c(out, paste(section_prefix, section))

    rows <- families[families$section == section, , drop = FALSE]
    for (i in seq_len(nrow(rows))) {
      family <- rows$family[[i]]
      display <- rows$display[[i]]
      status <- .validation_family_status(family)
      status_markup <- paste(
        vapply(status, function(x) {
          .validation_status_badge(x, family = family, tooltip = TRUE)
        }, character(1)),
        collapse = " "
      )

      out <- c(
        out,
        "",
        paste0(family_prefix, " [", display, "](", rows$url[[i]], ")"),
        "",
        status_markup,
        "",
        rows$description[[i]]
      )
    }

    out <- c(out, "")
  }

  out <- paste(out, collapse = "\n")

  if (requireNamespace("knitr", quietly = TRUE)) {
    knitr::asis_output(paste0(out, "\n"))
  } else {
    out
  }
}

#' @noRd
.validation_family_label <- function(x) {
  labels <- c(
    benchmark = "Benchmarked",
    benchmarked = "Benchmarked",
    validated = "Validated",
    partially_validated = "Partially validated",
    experimental = "Experimental",
    unvalidated = "Unvalidated"
  )

  if (x %in% names(labels)) labels[[x]] else x
}
