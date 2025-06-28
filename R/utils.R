# rdpartial internal utilities --------------------------------------------------
#
# This file centralises small helper functions that are used across the package
# but should remain **internal** (i.e. not exported).  All helpers are
# prepended with a leading dot so that they are only reachable via
# `rdpartial:::`, which makes accidental user calls unlikely.
#
# Each helper contains an explicit, minimal roxygen block with
# `@keywords internal` so that they are picked up by `R CMD check` without
# creating documentation clutter for end–users.
# -----------------------------------------------------------------------------

#' Assert that a scalar is numeric
#'
#' Thin wrapper around `stopifnot()` providing a clearer error message when an
#' argument is expected to be a single numeric value.
#'
#' @param x   Object to test.
#' @param nm  Character name used in the error message.
#' @keywords internal
.assert_scalar_numeric <- function(x, nm = deparse(substitute(x))) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x)) {
    stop(sprintf("`%s` must be a single, non‑NA numeric value.", nm), call. = FALSE)
  }
  invisible(TRUE)
}

#' Quick column–presence assertion
#'
#' Fails early (and verbosely) if required columns are missing from a data
#' frame.  Used throughout the package wherever we accept arbitrary user data.
#'
#' @param df   A `data.frame`.
#' @param cols Character vector of required column names.
#' @keywords internal
.check_columns <- function(df, cols) {
  missing <- setdiff(cols, names(df))
  if (length(missing)) {
    stop("Missing required column(s): ", paste(missing, collapse = ", "),
         call. = FALSE)
  }
  invisible(TRUE)
}

#' Column–intercept polynomial design matrix
#'
#' Creates a simple polynomial basis `1, x, x^2, …, x^order` that we use for the
#' local polynomial regressions on either side of the RDD cutoff.  The function
#' purposefully avoids `stats::poly()` because we do *not* want orthogonal
#' polynomials here—keeping the design untransformed eases interpretation and
#' reproducibility across languages.
#'
#' @param x     Numeric vector.
#' @param order Non‑negative integer; the highest power of `x` to include.
#' @return A numeric matrix with `length(x)` rows and `order + 1` columns.
#' @keywords internal
.poly_design <- function(x, order = 1L) {
  .assert_scalar_numeric(order, "order")
  if (order < 0 || order != round(order)) {
    stop("`order` must be a non‑negative integer.", call. = FALSE)
  }
  if (!is.numeric(x)) {
    stop("`x` must be numeric.", call. = FALSE)
  }

  # Fast path for order == 0 (intercept only) ------------------------------
  if (order == 0L) {
    return(matrix(1, nrow = length(x), ncol = 1L))
  }

  cbind(1, vapply(seq_len(order), function(p) x^p, numeric(length(x))))
}

#' Truncated tricube kernel weights
#'
#' Implements the weight function `(1 - (d / max(d))^3)^3` for non‑negative
#' distances `d`, returning `1` when `d == 0` and tapering smoothly to `0` at
#' `max(d)`.  This mirrors the weight choice in the original research code.
#'
#' @param dist Numeric vector of *non‑negative* distances.
#' @return Numeric vector of the same length, each element in `[0, 1]`.
#' @keywords internal
.tricube <- function(dist) {
  if (!is.numeric(dist) || any(is.na(dist)) || any(dist < 0)) {
    stop("`dist` must be a non‑negative numeric vector without NAs.", call. = FALSE)
  }
  if (length(dist) == 0L) return(numeric(0L))
  if (all(dist == 0))    return(rep(1, length(dist)))

  (1 - (dist / max(dist))^3)^3
}
