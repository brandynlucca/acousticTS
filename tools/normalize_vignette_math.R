normalize_preceding_line <- function(lines, start_index) {
  preceding <- seq_len(start_index - 1L)
  nonblank <- which(trimws(lines[preceding]) != "")
  if (length(nonblank) == 0L) {
    return(lines)
  }

  idx <- tail(nonblank, 1L)
  while (idx > 1L && is_anchor_line(lines[idx])) {
    lines[idx] <- sub(":$", "", trimws(lines[idx]))
    prior_nonblank <- which(trimws(lines[seq_len(idx - 1L)]) != "")
    if (length(prior_nonblank) == 0L) {
      return(lines)
    }
    idx <- tail(prior_nonblank, 1L)
  }

  while (idx > 1L && is_math_fence(lines[idx])) {
    lines[idx] <- "$$"
    prior_nonblank <- which(trimws(lines[seq_len(idx - 1L)]) != "")
    if (length(prior_nonblank) == 0L) {
      return(lines)
    }
    idx <- tail(prior_nonblank, 1L)
  }

  if ((idx + 1L) <= (start_index - 1L)) {
    between <- trimws(lines[seq.int(idx + 1L, start_index - 1L)])
    if (any(between %in% c("$$", "$$:"))) {
      return(lines)
    }
  }

  if (trimws(lines[idx]) == ":" && idx > 1L) {
    prior_nonblank <- which(trimws(lines[seq_len(idx - 1L)]) != "")
    if (length(prior_nonblank) == 0L) {
      return(lines)
    }
    prior_idx <- tail(prior_nonblank, 1L)
    lines[prior_idx] <- force_trailing_colon(lines[prior_idx])
    lines[idx] <- ""
    return(lines)
  }

  if (should_force_colon(lines[idx])) {
    lines[idx] <- force_trailing_colon(lines[idx])
  }

  lines
}

should_force_colon <- function(line) {
  trimmed <- trimws(line)
  if (trimmed == "" || grepl(":$", trimmed)) {
    return(FALSE)
  }

  !grepl("^(#|\\* |\\- |\\| |```|~~~|:::)", trimmed) &&
    !is_anchor_line(trimmed) &&
    !is_math_fence(trimmed)
}

force_trailing_colon <- function(line) {
  trimmed <- sub("[[:space:]]+$", "", line)
  if (grepl(":$", trimmed)) {
    return(trimmed)
  }
  if (grepl("[.;,]$", trimmed)) {
    return(sub("[.;,]$", ":", trimmed))
  }
  paste0(trimmed, ":")
}

is_anchor_line <- function(line) {
  grepl("^<div id=\"[^\"]+\"></div>:?$", trimws(line))
}

is_math_fence <- function(line) {
  trimws(line) %in% c("$$", "$$:")
}

normalize_math_block <- function(lines) {
  compact <- character()

  for (line in lines) {
    stripped <- trimws(sub("[[:space:]]+$", "", line))
    stripped <- sub(":$", "", stripped)
    if (stripped == "") {
      compact <- c(compact, "")
      next
    }

    if (stripped %in% c("=", "+", "-")) {
      if (length(compact) > 0L) {
        compact[length(compact)] <- paste0(
          sub("[[:space:]]+$", "", compact[length(compact)]),
          " ",
          stripped
        )
      } else {
        compact <- c(compact, stripped)
      }
      next
    }

    if (grepl("^[=+\\-][[:space:]]+", stripped)) {
      operator <- substr(stripped, 1L, 1L)
      remainder <- trimws(sub("^[=+\\-][[:space:]]*", "", stripped))
      if (length(compact) > 0L) {
        compact[length(compact)] <- paste0(
          sub("[[:space:]]+$", "", compact[length(compact)]),
          " ",
          operator
        )
        compact <- c(compact, remainder)
      } else {
        compact <- c(compact, stripped)
      }
      next
    }

    compact <- c(compact, stripped)
  }

  compact <- strip_recycled_prefix(compact)

  normalized <- character()
  continuation <- FALSE

  for (line in compact) {
    if (line == "") {
      normalized <- c(normalized, "")
      next
    }

    indent <- if (continuation) "    " else "  "
    normalized <- c(normalized, paste0(indent, line))
    continuation <- grepl("[=+\\-]$", line)
  }

  normalized
}

strip_recycled_prefix <- function(lines) {
  if (length(lines) < 2L) {
    return(lines)
  }

  max_prefix <- floor((length(lines) - 1L) / 2L)
  for (k in seq.int(max_prefix, 1L)) {
    if (identical(utils::tail(lines, k), utils::head(lines, k))) {
      return(lines[seq_len(length(lines) - k)])
    }
  }

  lines
}

normalize_file <- function(path) {
  lines <- readLines(path, warn = FALSE)
  if (!length(lines)) {
    return(invisible(FALSE))
  }

  in_yaml <- length(lines) >= 1L && identical(trimws(lines[1L]), "---")
  in_code <- FALSE
  fence_marker <- NULL
  changed <- FALSE
  i <- 1L

  while (i <= length(lines)) {
    trimmed <- trimws(lines[i])

    if (in_yaml) {
      if (i > 1L && identical(trimmed, "---")) {
        in_yaml <- FALSE
      }
      i <- i + 1L
      next
    }

    if (!in_code && grepl("^(?:```|~~~)", trimmed, perl = TRUE)) {
      in_code <- TRUE
      fence_marker <- substr(trimmed, 1L, 3L)
      i <- i + 1L
      next
    }

    if (in_code) {
      if (startsWith(trimmed, fence_marker)) {
        in_code <- FALSE
        fence_marker <- NULL
      }
      i <- i + 1L
      next
    }

    if (trimmed %in% c("$$", "$$:")) {
      lines[i] <- "$$"
      end <- i + 1L
      while (end <= length(lines) &&
             !(trimws(lines[end]) %in% c("$$", "$$:"))) {
        end <- end + 1L
      }

      if (end <= length(lines)) {
        lines[end] <- "$$"
        updated <- normalize_preceding_line(lines, i)
        if (!identical(updated, lines)) {
          lines <- updated
          changed <- TRUE
        }

        block <- lines[(i + 1L):(end - 1L)]
        normalized <- normalize_math_block(block)
        if (!identical(block, normalized)) {
          prefix <- if (i >= 1L) {
            lines[seq_len(i)]
          } else {
            character()
          }
          suffix <- if (end <= length(lines)) {
            lines[end:length(lines)]
          } else {
            character()
          }
          lines <- c(prefix, normalized, suffix)
          changed <- TRUE
          end <- i + length(normalized) + 1L
        }
        i <- end + 1L
        next
      }
    }

    if (should_dedent_prose(lines, i)) {
      lines[i] <- sub("^  ", "", lines[i])
      changed <- TRUE
    }

    i <- i + 1L
  }

  if (changed) {
    writeLines(lines, path, useBytes = TRUE)
  }

  invisible(changed)
}

should_dedent_prose <- function(lines, index) {
  line <- lines[index]
  if (!grepl("^  \\S", line)) {
    return(FALSE)
  }

  trimmed <- trimws(line)
  if (grepl("^([-*+] |\\d+\\. |\\| |<|:::)", trimmed)) {
    return(FALSE)
  }

  previous_nonblank <- which(trimws(lines[seq_len(index - 1L)]) != "")
  if (length(previous_nonblank) > 0L) {
    previous <- trimws(lines[tail(previous_nonblank, 1L)])
    if (grepl("^([-*+] |\\d+\\. )", previous)) {
      return(FALSE)
    }
  }

  TRUE
}

paths <- list.files("vignettes", recursive = TRUE, pattern = "\\.Rmd$",
                    full.names = TRUE)
changes <- vapply(paths, normalize_file, logical(1))
cat(sum(changes), "vignette files updated\n")
