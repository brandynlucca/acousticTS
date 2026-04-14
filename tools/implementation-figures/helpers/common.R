impl_set_repo_root <- function() {
  has_repo_root <- function(path) {
    file.exists(file.path(path, "DESCRIPTION")) &&
      file.exists(file.path(
        path,
        "tools",
        "implementation-figures",
        "manifest.csv"
      ))
  }

  parent_dir <- function(path) {
    normalizePath(file.path(path, ".."), winslash = "/", mustWork = FALSE)
  }

  search_upward <- function(start_path) {
    current <- normalizePath(start_path, winslash = "/", mustWork = FALSE)
    repeat {
      if (has_repo_root(current)) {
        return(current)
      }

      next_path <- parent_dir(current)
      if (identical(next_path, current)) {
        return(NULL)
      }
      current <- next_path
    }
  }

  search_roots <- character()

  if (nzchar(Sys.getenv("GITHUB_WORKSPACE", unset = ""))) {
    search_roots <- c(search_roots, Sys.getenv("GITHUB_WORKSPACE"))
  }

  search_roots <- c(search_roots, getwd())

  args <- commandArgs(FALSE)
  file_arg <- args[grep("^--file=", args)]
  if (length(file_arg) == 1L) {
    search_roots <- c(search_roots, dirname(sub("^--file=", "", file_arg)))
  }

  search_roots <- unique(search_roots[nzchar(search_roots)])
  for (root in search_roots) {
    repo_root <- search_upward(root)
    if (!is.null(repo_root)) {
      setwd(repo_root)
      return(invisible(repo_root))
    }
  }

  stop(
    "Could not locate the repository root containing ",
    "'tools/implementation-figures/manifest.csv'.",
    call. = FALSE
  )
}

impl_load_all <- function() {
  impl_set_repo_root()
  devtools::load_all(".", quiet = TRUE)
  invisible(TRUE)
}

impl_data_path <- function(name) {
  file.path("tools", "implementation-figures", "data", name)
}

impl_output_path <- function(family, name) {
  file.path("vignettes", family, name)
}

impl_with_png <- function(path,
                          code,
                          width = 1600,
                          height = 1200,
                          res = 200) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  grDevices::png(path, width = width, height = height, res = res)
  on.exit(grDevices::dev.off(), add = TRUE)
  force(code)
  invisible(path)
}

impl_ts_label <- function() {
  expression(Target ~ strength ~ (dB ~ re. ~ 1 ~ m^2))
}

impl_delta_label <- function() {
  expression(Delta ~ TS ~ (dB))
}

impl_heatmap_palette <- function(n = 101) {
  grDevices::colorRampPalette(c("#2166ac", "#f7f7f7", "#b2182b"))(n)
}

impl_draw_colorbar <- function(zlim,
                               palette_vals = impl_heatmap_palette(),
                               label = impl_delta_label(),
                               n_ticks = 5) {
  zlim <- range(zlim, finite = TRUE)
  if (length(zlim) != 2L || any(!is.finite(zlim))) {
    zlim <- c(-1, 1)
  }
  if (diff(zlim) <= 0) {
    zlim <- c(-1, 1)
  }

  usr <- graphics::par("usr")
  x_range <- usr[2] - usr[1]
  y_range <- usr[4] - usr[3]
  if (!is.finite(x_range) || x_range == 0) {
    x_range <- 1
  }
  if (!is.finite(y_range) || y_range == 0) {
    y_range <- 1
  }

  x_left <- usr[2] + 0.05 * x_range
  x_right <- usr[2] + 0.11 * x_range
  y_bottom <- usr[3] + 0.08 * y_range
  y_top <- usr[4] - 0.08 * y_range

  old_xpd <- graphics::par("xpd")[[1]]
  graphics::par(xpd = NA)
  on.exit(graphics::par(xpd = old_xpd), add = TRUE)

  y_edges <- seq(y_bottom, y_top, length.out = length(palette_vals) + 1L)
  for (i in seq_along(palette_vals)) {
    graphics::rect(
      xleft = x_left,
      ybottom = y_edges[[i]],
      xright = x_right,
      ytop = y_edges[[i + 1L]],
      col = palette_vals[[i]],
      border = NA
    )
  }
  graphics::rect(
    xleft = x_left,
    ybottom = y_bottom,
    xright = x_right,
    ytop = y_top,
    border = "gray30",
    lwd = 0.8
  )

  ticks <- pretty(zlim, n = n_ticks)
  ticks <- ticks[ticks >= zlim[[1]] & ticks <= zlim[[2]]]
  if (!length(ticks)) {
    ticks <- zlim
  }
  y_ticks <- y_bottom + (ticks - zlim[[1]]) / diff(zlim) * (y_top - y_bottom)
  tick_extent <- 0.015 * x_range
  graphics::segments(
    x0 = x_right,
    y0 = y_ticks,
    x1 = x_right + tick_extent,
    y1 = y_ticks
  )
  graphics::text(
    x = x_right + 0.025 * x_range,
    y = y_ticks,
    labels = format(signif(ticks, 3), trim = TRUE),
    adj = c(0, 0.5),
    cex = 0.68
  )
  graphics::text(
    x = (x_left + x_right) / 2,
    y = y_top + 0.06 * y_range,
    labels = label,
    cex = 0.72
  )
}

impl_sort_compare_df <- function(df) {
  sort_cols <- intersect(
    c(
      "stage",
      "validation_case",
      "boundary_condition",
      "case",
      "size_scale",
      "frequency_hz",
      "frequency",
      "frequency_khz",
      "angle_deg",
      "theta_body_deg",
      "theta_deg",
      "phi_deg",
      "ka"
    ),
    names(df)
  )

  if (!length(sort_cols)) {
    return(df)
  }

  order_args <- lapply(sort_cols, function(col) df[[col]])
  row_idx <- do.call(order, c(order_args, list(na.last = TRUE)))
  df[row_idx, , drop = FALSE]
}

impl_should_refresh_timings <- function() {
  refresh_override <- tolower(
    Sys.getenv("ACOUSTICTS_IMPL_REFRESH_TIMINGS", unset = "")
  )

  identical(refresh_override, "true") ||
    !nzchar(Sys.getenv("CI", unset = ""))
}

impl_round_timing_columns <- function(df, digits = 6) {
  timing_cols <- grep("elapsed", names(df), value = TRUE)

  for (col in timing_cols) {
    if (is.numeric(df[[col]])) {
      df[[col]] <- round(df[[col]], digits = digits)
    }
  }

  df
}
