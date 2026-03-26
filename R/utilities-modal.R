################################################################################
################################################################################
# MODAL-SERIES UTILITIES
################################################################################
################################################################################
#' Pad modal coefficient vector to a target length
#' @param x Modal coefficient vector.
#' @param target_length Desired output length.
#' @param fill Fill value used for padding.
#' @return Vector padded or truncated to `target_length`.
#' @keywords internal
#' @noRd
.pad_modal_terms <- function(x, target_length, fill = NA) {
  # Truncate vectors that already exceed the requested modal cutoff ============
  if (length(x) >= target_length) {
    return(x[seq_len(target_length)])
  }

  # Promote the fill value to complex storage when needed ======================
  fill_value <- if (is.complex(x) || is.complex(fill)) {
    as.complex(fill)
  } else {
    fill
  }

  # Pad the coefficient vector out to the common target length =================
  c(x, rep(fill_value, target_length - length(x)))
}

#' Apply a modal-coefficient builder over frequencies and pad to a matrix
#'
#' @param m_limit Integer vector of per-frequency modal truncation limits.
#' @param FUN Builder returning an unpadded modal vector for one frequency.
#' @param ... Additional vectors recycled by `mapply`.
#' @param fill Fill value used to pad short modal vectors.
#' @return Matrix with one row per retained mode and one column per frequency.
#' @keywords internal
#' @noRd
.modal_series_apply <- function(m_limit, FUN, ..., fill = NA) {
  # Resolve the common retained modal length across all frequencies ============
  target_length <- max(m_limit) + 1L
  out <- mapply(
    FUN = function(ml, ...) {
      .pad_modal_terms(FUN(..., ml), target_length, fill = fill)
    },
    ml = m_limit,
    ...,
    SIMPLIFY = TRUE
  )

  # Restore matrix structure when only one frequency is evaluated ==============
  if (is.null(dim(out))) {
    out <- matrix(out, nrow = target_length)
  }

  # Return the padded modal matrix =============================================
  out
}

#' Compute weighted modal sums column-wise
#' @param coeff_matrix Modal coefficient matrix with modes along rows.
#' @param weights Vector of modal weights.
#' @param na.rm Whether to remove `NA` values.
#' @return Numeric or complex vector of column-wise weighted sums.
#' @keywords internal
#' @noRd
.modal_weighted_sum <- function(coeff_matrix, weights, na.rm = TRUE) {
  # Apply the mode weights row-wise before summing by frequency ================
  colSums(sweep(coeff_matrix, 1, weights, FUN = "*"), na.rm = na.rm)
}

################################################################################
#' Format data for the modal series solution model into the appropriate matrix
#' @param v Vector input.
#' @param limit Modal series limit.
#' @keywords internal
#' @noRd
modal_matrix <- function(v, limit) {
  # Replicate one frequency vector into the common modal matrix layout =========
  matrix(
    data = rep(v, each = limit + 1),
    ncol = length(v),
    nrow = limit + 1
  )
}
