################################################################################
# INTERNAL VALIDATION / BENCHMARK REGISTRY
################################################################################

#' @noRd
.validation_status_badge <- function(status) {
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
      paste("Prolate-spheroidal modal-series solution for smooth",
            "elongated canonical bodies."),
      "Solid elastic spherical model used mainly for calibration spheres.",
      "Elastic-shelled spherical family for layered shell targets.",
      paste("Bent-cylinder modal-series family for straight and uniformly bent",
      "cylinders."),
      "Elastic-cylinder modal-series family for fully elastic solid cylinders.",
      "Weak-scattering elongated-body approximation for fluid-like targets.",
      "Stochastic DWBA family for unresolved phase variability.",
      paste("Kirchhoff-ray mode model for segmented fish-like ",
      "body-plus-inclusion targets."),
      "High-pass approximation for compact asymptotic screening.",
      paste("Two-ray cylindrical family for high-frequency locally",
            "cylindrical targets."),
      "Phase-compensated DWBA for bent weakly scattering targets.",
      paste("Composite flesh-plus-backbone family for swimbladder-less ",
      "fish-like targets."),
      paste("Viscous-elastic layered-sphere family for gas-core, shell,",
            "and viscous-layer targets."),
      paste("Single-target transition-matrix family for retained monostatic",
            "and angle-dependent scattering products across supported",
            "canonical shapes.")
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
  data.frame(
    family = c(
      "sphms", "sphms",
      "fcms", "fcms",
      "psms", "psms",
      "calibration", "calibration",
      "essms",
      "bcms",
      "ecms",
      "dwba", "dwba",
      "sdwba", "sdwba",
      "krm", "krm",
      "hpa", "hpa",
      "trcm",
      "pcdwba",
      "bbfm",
      "vesm",
      "tmm", "tmm", "tmm", "tmm"
    ),
    evidence_type = c(
      "benchmark", "external",
      "benchmark", "external",
      "benchmark", "external",
      "benchmark", "external",
      "note",
      "note",
      "note",
      "benchmark", "external",
      "benchmark", "external",
      "benchmark", "external",
      "benchmark", "external",
      "benchmark",
      "external",
      "note",
      "external",
      "benchmark", "external", "external", "note"
    ),
    source = c(
      # SPHMS
      "benchmark_ts / Jech et al. (2015)", "KRMr and echoSMs",
      # FCMS
      "benchmark_ts / Jech et al. (2015)", "echoSMs",
      # PSMS
      "benchmark_ts / Jech et al. (2015)", "Prol_Spheroid",
      # CAL
      "Published calibration spheres", "echoSMs, sphereTS, NOAA applet",
      # ESSMS
      "Jech shell-sphere benchmark family",
      # BCMS
      "Internal FCMS-based reference reconstruction",
      # ECMS
      "Independent algebra transcription",
      # DWBA
      "Jech weakly scattering benchmark ladder",
      "Published and independent DWBA implementations",
      # SDWBA
      "Published SDWBA weak-scattering ladder",
      "CCAMLR MATLAB and NOAA HTML implementations",
      # KRM
      "Canonical modal families", "KRMr, echoSMs, NOAA applet",
      # HPA
      "Jech canonical asymptotic targets",
      "echoSMs HPModel and published algebra",
      # TRCM
      "Package validation workflow",
      # PCDWBA
      "ZooScatR and echopop source workflows",
      # BBFM
      "Internal DWBA + ECMS reconstruction",
      # VESM
      "Reference Python VESM workflow",
      # TTM
      "SPHMS / PSMS / FCMS benchmark ladder", "BEMPP far-field checks",
      "Exact general-angle spheroidal solution",
      "Cylinder branch remains guarded"
    ),
    summary = c(
      paste("Benchmarked against the canonical spherical spectra stored",
            "in `benchmark_ts`."),
      paste("Validated against `KRMr` and `echoSMs` on shared",
            "penetrable-sphere cases."),
      paste("Benchmarked against the canonical finite-cylinder spectra",
            "stored in `benchmark_ts`."),
      "Validated against the `echoSMs` finite-cylinder implementation.",
      paste("Benchmarked against the canonical prolate-spheroid spectra",
            "stored in `benchmark_ts`."),
      paste("Validated against the external `Prol_Spheroid` implementation",
            "on shared prolate cases."),
      paste("Benchmarked against published calibration-sphere targets used",
            "throughout the package documentation."),
      paste("Validated against `echoSMs`, `sphereTS`, and the NOAA calibration",
            "applet."),
      paste("A direct shell-sphere benchmark family exists, but the current",
            "ESSMS implementation still does not return finite full-grid",
            "benchmark spectra."),
      paste("Current BCMS checks are internal coherence reconstructions; the",
            "family does not yet have an external benchmark or software ladder."
            ),
      paste("Current ECMS checks are independent algebra reconstructions",
            "rather than a documented external benchmark ladder."),
      paste("Benchmarked against the canonical weakly scattering targets",
            "summarized by Jech et al. (2015)."),
      paste("Validated against the published McGehee MATLAB workflow and an",
            "independent DWBA implementation."),
      "Benchmarked against published SDWBA weak-scattering comparison cases.",
      paste(
        "Validated against the CCAMLR MATLAB and NOAA HTML SDWBA ",
        "implementations."
      ),
      paste("Benchmarked against canonical modal-family targets used for",
            "isolated gas-filled and weakly scattering cases."),
      paste("Validated against `KRMr`, `echoSMs`, and the NOAA KRM applet on",
            "bundled fish objects and shared workflows."),
      paste("Benchmarked against canonical asymptotic target families rather",
            "than as an exact modal solver."),
      paste("Validated against the spherical `echoSMs::HPModel` branch and",
            "the published Johnson/Stanton algebra."),
      paste("Benchmarked within the package validation workflow against the",
            "straight-cylinder and FCMS-derived bent-cylinder reference",
            "constructions."),
      paste("Validated against source-level `ZooScatR` and `echopop` PCDWBA ",
            "workflows."),
      paste("BBFM currently has documented internal reconstruction checks",
            "but no external benchmark ladder or independent public",
            "implementation comparison."),
      paste("Validated against the reference Python VESM implementation on",
            "the documented layered-sphere case."),
      paste("Benchmarked against `SPHMS`, `PSMS`, and `FCMS` on the",
            "currently supported canonical shape branches."),
      paste("Validated against external BEMPP far-field checks for sphere,",
            "oblate, and prolate pressure-release cases."),
      paste("Retained prolate angular products are also checked against the",
            "exact general-angle spheroidal solution."),
      paste("The cylinder branch is benchmark-matched only for the exact",
            "monostatic workflow; retained general-angle cylinder products",
            "remain outside the validated public scope.")
    ),
    scope = c(
      "Sphere spectra across rigid, soft, liquid-filled, and gas-filled cases.",
      "Penetrable sphere spectra on shared software definitions.",
      paste("Finite-cylinder spectra across the canonical cylindrical",
            "benchmark grid."),
      "Rigid, soft, liquid-filled, and gas-filled finite-cylinder spectra.",
      "Prolate-spheroid spectra across the canonical benchmark grid.",
      "Liquid-filled and gas-filled prolate-spheroid software comparisons.",
      "Tungsten-carbide and copper calibration spheres.",
      "Shared calibration-sphere material sets and frequency sweeps.",
      paste("Layered shell-sphere benchmark status only; not yet",
            "benchmark-grade agreement."),
      "Uniform-curvature cylinder coherence extension of FCMS.",
      paste("Elastic-cylinder component family and near-broadside",
            "canonical cases."),
      "Weakly scattering sphere, prolate spheroid, and cylinder targets.",
      "Bundled krill geometry and published DWBA reference workflows.",
      paste("Weakly scattering sphere, prolate spheroid, and cylinder",
            "stochastic targets."),
      "Bundled krill stochastic workflow comparisons.",
      paste("Canonical isolated targets used for the package KRM",
            "benchmark ladder."),
      "Bundled sardine and cod software-to-software comparisons.",
      paste("Sphere, prolate spheroid, and cylinder asymptotic",
            "benchmark targets."),
      "Spherical HPModel branch and published asymptotic formulas.",
      paste("Straight and bent cylindrical validation cases documented",
            "in the package."),
      paste("Curved weak-scattering reference workflows on shared",
            "bent-body cases."),
      "Internal composite-component consistency checks only.",
      paste("Documented spherical layered case used by the original",
            "VESM implementation."),
      "Sphere, oblate, prolate, and guarded cylinder monostatic branches.",
      "Pressure-release angular slices for sphere, oblate, and prolate cases.",
      "General-angle prolate retained-state validation.",
      "Cylinder retained-angle scope limitation and guardrails."
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
  evidence <- .validation_family_evidence(family)

  status <- character()

  if (any(evidence$evidence_type == "benchmark")) {
    status <- c(status, "benchmarked")
  }
  if (any(evidence$evidence_type == "external")) {
    status <- c(status, "validated")
  }
  if (isTRUE(meta$experimental)) {
    status <- c(status, "experimental")
  }
  if (!any(evidence$evidence_type %in% c("benchmark", "external"))) {
    status <- c(status, "unvalidated")
  }

  status
}

#' @noRd
.validation_family_validation <- function(family) {
  evidence <- .validation_family_evidence(
    family,
    type = c("benchmark", "external", "note")
  )

  if (nrow(evidence) == 0) {
    character()
  } else {
    evidence$summary
  }
}

#' @noRd
.validation_status_policy <- function() {
  badge <- function(x) .validation_status_badge(x)

  sections <- c(
    paste("Status tags are derived from the package's internal validation",
          "registry and used conservatively:"),
    paste0(
      "- ", badge("benchmarked"),
      " means the family has a documented comparison against a canonical ",
      "benchmark ladder or stored benchmark values."
    ),
    paste0(
      "- ", badge("validated"),
      " means the family has a documented comparison against at least one ",
      "external implementation or software package."
    ),
    paste0(
      "- ", badge("experimental"),
      " means the family is available to use, but its interface or ",
      "validation scope should still be treated as provisional."
    ),
    paste0(
      "- ", badge("unvalidated"),
      " means the package site does not yet document either benchmark ",
      "evidence or an external implementation comparison for that family."
    ),
    "",
    "These statuses are not all mutually exclusive:",
    paste("- `Benchmarked` and `Validated` are evidence tags and can",
          "appear together."),
    paste("- `Experimental` is a lifecycle tag and can coexist with",
          "either evidence tag."),
    paste("- `Unvalidated` is reserved for families lacking both evidence",
          "tags, so it should not be combined with other evidence tags.")
  )

  out <- paste(sections, collapse = "\n")

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
.validation_evidence_table <- function(type = c("benchmark", "external"),
                                       include_notes = FALSE) {
  evidence <- .validation_evidence_registry()
  type <- match.arg(type, several.ok = TRUE)

  keep_types <- unique(c(type, if (include_notes) "note"))
  evidence <- evidence[evidence$evidence_type %in% keep_types, ,
                       drop = FALSE]

  families <- .validation_family_registry()[, c("family", "display"),
                                            drop = FALSE]
  evidence <- merge(evidence, families, by = "family", all.x = TRUE,
                    sort = FALSE)

  evidence_type_order <- c(benchmark = 1L, external = 2L, note = 3L)
  evidence <- evidence[order(
    match(evidence$family, .validation_family_registry()$family),
    evidence_type_order[evidence$evidence_type]
  ), , drop = FALSE]

  data.frame(
    Family = evidence$display,
    `Evidence type` = vapply(evidence$evidence_type,
                             .validation_family_label, character(1)),
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
        vapply(status, .validation_status_badge, character(1)),
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
    external = "Validated",
    validated = "Validated",
    experimental = "Experimental",
    note = "Note",
    unvalidated = "Unvalidated"
  )

  if (x %in% names(labels)) labels[[x]] else x
}
