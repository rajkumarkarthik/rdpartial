#' Parametric Bootstrap for Partial-Identification Bounds
#'
#' Draw repeated resamples of the original data **keeping the running-variable
#' distribution fixed** and re-estimate the sharp or fuzzy RDD bounds for each
#' user-specified manipulation region.
#'
#' This function is a high-level wrapper around three internal work-horses:
#' * `.density_estimation()` - estimates non-manipulated counts (`true_counts`).
#' * [bounds_sharp()] / [bounds_fuzzy()] - compute the Manski bounds.
#' * `.tricube()` - constructs local-kernel weights used inside the bounds
#'   estimators.
#'
#' @param data A `data.frame` that contains the running variable, outcome and-
#'   for fuzzy designs-the treatment column.
#' @param running_var,outcome Character names of the running variable and
#'   outcome columns.
#' @param treatment Character name of the treatment take-up column (0/1).
#'   Required when `estimator = "fuzzy"`.
#' @param cutoff Numeric RDD threshold.
#' @param manip_regions List of numeric length-2 vectors `(lower, upper)` giving
#'   suspected manipulation intervals.
#' @param estimator Either "fuzzy" (default) or "sharp".
#' @param treat_direction Either "increase" (treatment prob. jumps *up* at the
#'   cutoff – standard case) or "decrease" (jumps *down*, e.g. donor *deferral*
#'   example).  Determines the sign of the optimisation objective.
#' @param n_boot Integer number of bootstrap replications (default `200`).
#' @param poly_order Local polynomial order (default `1`).
#' @param weight_var Optional character column in `data` holding non-negative
#'   observation weights.  If `NULL`, each row receives weight 1.
#' @param density_args Optional named list forwarded to `.density_estimation()`.
#' @param ci_level Percentile coverage of the bootstrap interval (default
#'   `0.95`).
#' @param parallel Logical.  If `TRUE`, uses `parallel::mclapply()` on
#'   non-Windows systems.
#' @param n_cores Number of cores for parallel execution (defaults to
#'   `parallel::detectCores() - 1`).
#' @param progress Logical - print a progress bar (default `TRUE`).  Progress is
#'   suppressed automatically when running in parallel to avoid garbled output.
#' @param seed Optional integer for reproducibility of resampling.
#'
#' @return An object of class `rdpartial_boot` - a list with elements
#' * `boot_array` - 3-D array `(n_boot, R, 2)` storing lower/upper bounds.
#' * `ci`         - `R × 2` matrix of percentile intervals (`lwr`, `upr`).
#' * `manip_regions`, `estimator`, `call` - metadata for downstream methods.
#'
#' @export
#' @examples
#' set.seed(101)
#' n <- 3000; cutoff <- 16
#' x <- rpois(n, 15)
#' z <- rbinom(n, 1, prob = ifelse(x >= cutoff, 0.7, 0.2))
#' y <- 0.2 * x + 3 * z + rnorm(n)
#' dat <- data.frame(x = x, y = y, z = z)
#'
#' res <- bootstrap_bounds(dat, running_var = "x", outcome = "y",
#'                         treatment = "z", cutoff = cutoff,
#'                         manip_regions = list(c(cutoff - 0.3, cutoff)),
#'                         n_boot = 50, parallel = FALSE)
#' res$ci
#'
#' @importFrom stats quantile
bootstrap_bounds <- function(data, running_var, outcome, treatment = NULL,
                             cutoff, manip_regions,
                             estimator = c("fuzzy", "sharp"),
                             treat_direction = c("increase", "decrease"),
                             n_boot = 200L, poly_order = 1L, weight_var = NULL,
                             density_args = list(), ci_level = 0.95,
                             parallel = TRUE, n_cores = NULL,
                             progress = TRUE, seed = NULL) {

  # Sanity checks -------------------------------------------------------
  if (!is.data.frame(data)) stop("`data` must be a data.frame.")
  .check_columns(data, c(running_var, outcome))

  estimator <- match.arg(estimator)
  if (estimator == "fuzzy") {
    if (is.null(treatment)) stop("`treatment` required for fuzzy estimator.")
    .check_columns(data, treatment)
  }

  .assert_scalar_numeric(cutoff, "cutoff")

  if (!is.list(manip_regions) ||
      !all(vapply(manip_regions, length, 1L) == 2)) {
    stop("`manip_regions` must be a list of length-2 numeric vectors.")
  }

  if (!is.null(seed)) set.seed(seed)

  # User-supplied observation weights ---------------------------------------
  user_wts <- if (is.null(weight_var)) rep(1, nrow(data)) else {
    .check_columns(data, weight_var)
    w <- data[[weight_var]]
    if (any(w < 0, na.rm = TRUE)) stop("`weight_var` contains negative values.")
    w[is.na(w)] <- 0
    w
  }
  names(user_wts) <- rownames(data)

  rv_full <- data[[running_var]]
  # Density estimation on original sample ------------------------------
  hist_df <- as.data.frame(table(rv_full))
  names(hist_df) <- c("x", "freq")
  hist_df$x <- as.numeric(as.character(hist_df$x))

  density_list <- lapply(manip_regions, function(rg) {
    do.call(.density_estimation,
            c(list(hist_df = hist_df, manip_region = rg, cutoff = cutoff),
              density_args))
  })

  # Bootstrap infrastructure -------------------------------------------
  strata <- split(seq_len(nrow(data)), rv_full)        # keep RV dist fixed

  resample_idx_once <- function() {
    unlist(lapply(strata, function(ix) {
      if (length(ix) == 1L) ix else sample(ix, replace = TRUE)
    }), use.names = FALSE)
  }

  boot_indices <- lapply(seq_len(n_boot), function(i) resample_idx_once())

  ## Helper now receives an index vector instead of a data.frame
  eval_bounds <- function(idx_vec) {
    df <- data[idx_vec, , drop = FALSE]
    rv <- df[[running_var]]

    # Kernel weights (tricube) mirroring original script -------------------
    left  <- rv <  cutoff
    right <- !left

    wt_kernel <- numeric(length(rv))
    if (any(left))  wt_kernel[left]  <- .tricube(cutoff - rv[left])
    if (any(right)) wt_kernel[right] <- .tricube(rv[right] - cutoff)

    row_idx <- as.integer(sub("\\..*$", "", rownames(df)))  # strip ".1", ".2", ...
    weights_vec <- wt_kernel * user_wts[row_idx]
    weights_vec <- .sanitize_weights(weights_vec)

    sapply(seq_along(manip_regions), function(j) {
      tc <- density_list[[j]]$counts

      if (estimator == "sharp") {
        bounds_sharp(x          = rv,
                     y          = df[[outcome]],
                     cutoff     = cutoff,
                     weights    = weights_vec,
                     poly_order = poly_order,
                     true_counts = tc,
                     bounds      = "both",
                     treat_direction = treat_direction)
      } else {
        bounds_fuzzy(x          = rv,
                     y          = df[[outcome]],
                     z          = df[[treatment]],
                     cutoff     = cutoff,
                     weights    = weights_vec,
                     poly_order = poly_order,
                     true_counts = tc,
                     bounds      = "both",
                     treat_direction = treat_direction)
      }
    })  # returns 2 × R matrix
  }

  # Run bootstrap -------------------------------------------------------
  if (parallel && .Platform$OS.type != "windows") {
    if (is.null(n_cores)) n_cores <- min(4L, parallel::detectCores() - 1L)  # Cap at 4 cores for memory safety
    runner   <- function(X, FUN) parallel::mclapply(X, FUN, mc.cores = n_cores)
    progress <- FALSE  # suppress progress bar under fork
  } else {
    runner <- lapply
  }

  if (progress) {
    pb <- utils::txtProgressBar(min = 0, max = n_boot, style = 3)
    on.exit(close(pb), add = TRUE)
  }

  boot_list <- runner(seq_len(n_boot), function(i) {
    res <- eval_bounds(boot_indices[[i]])
    if (progress) utils::setTxtProgressBar(pb, i)
    res
  })

  # Convert to array: (n_boot, R, 2) ---------------------------------------
  R <- length(manip_regions)
  boot_array <- array(NA_real_, dim = c(n_boot, R, 2),
                      dimnames = list(NULL,
                                      paste0("region", seq_len(R)),
                                      c("lower", "upper")))
  for (i in seq_len(n_boot)) boot_array[i, , ] <- t(boot_list[[i]])

  # Percentile CIs ------------------------------------------------------
  pr <- c((1 - ci_level) / 2, 1 - (1 - ci_level) / 2)
  ci_mat <- matrix(NA_real_, nrow = R, ncol = 2,
                   dimnames = list(paste0("region", seq_len(R)), c("lwr", "upr")))
  for (j in seq_len(R)) {
    ci_mat[j, 1] <- stats::quantile(boot_array[, j, 1], probs = pr[1], na.rm = TRUE)
    ci_mat[j, 2] <- stats::quantile(boot_array[, j, 2], probs = pr[2], na.rm = TRUE)
  }
  attr(ci_mat, "level") <- ci_level

  structure(list(boot_array    = boot_array,
                 ci            = ci_mat,
                 manip_regions = manip_regions,
                 estimator     = estimator,
                 call          = match.call()),
            class = "rdpartial_boot")
}

# S3 print method -----------------------------------------------------------
#' @export
print.rdpartial_boot <- function(x, ...) {
  lvl <- format(attr(x$ci, "level") * 100, digits = 3)
  cat("Bootstrap percentile CI (", lvl, " %):\n", sep = "")
  print(x$ci, ...)
  invisible(x)
}
