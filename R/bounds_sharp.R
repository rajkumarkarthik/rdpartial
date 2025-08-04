#' Partial‑Identification Bounds for **Sharp** RDDs with Manipulation
#'
#' Compute Manski–type lower/upper bounds on the causal effect at a Regression
#' Discontinuity Design (*sharp* treatment assignment) when the running
#' variable can be strategically manipulated to the right of the cutoff.  The
#' implementation mirrors Algorithm 1 in Rosenman *et al.* (2025) and is the
#' foundation for the fuzzy variant [bounds_fuzzy()].
#'
#' @section Required Inputs:
#' * **`x`** – numeric running variable.
#' * **`y`** – outcome measured under *realised* treatment status (because the
#'   design is sharp, no separate `z` argument is needed).
#' * **`cutoff`** – numeric threshold separating control (< cutoff) from
#'   treatment (≥ cutoff).
#' * **`true_counts`** – `data.frame` with columns `x` (support points ≥ cutoff)
#'   and `n_true` (estimated number of **non‑manipulated** observations at each
#'   support point).  Typically produced by [.density_estimation()].
#'
#' @param x Numeric running variable.
#' @param y Numeric outcome variable.
#' @param cutoff Numeric threshold separating control (< cutoff) from treatment (>= cutoff).
#' @param true_counts Data.frame with columns `x` (support points >= cutoff) and `n_true` (estimated number of non-manipulated observations at each support point).
#' @param weights Numeric vector of observation weights (optional).
#' @param poly_order Integer polynomial order for local regression (default 1L).
#' @param bounds Character – which bound(s) to return; one of
#'   "both" (default), "lower", "upper".
#' @param solver Character – CVXR solver (default taken from
#'   `getOption("rdpartial.solver")`; package default is "ECOS").
#' @param ...    Further arguments passed on to [CVXR::solve()].  Use this to
#'   override CVXR tolerances or provide a different solver control list.
#'
#' @return A numeric vector of length 2 (`lower`, `upper`) when `bounds =
#'   "both"`; otherwise a single numeric.  The full `CVXR` solution objects are
#'   attached as attributes `opt_lwr` / `opt_upr` for diagnostic purposes.
#'
#' @export
#' @examples
#' set.seed(123)
#' n  <- 4000
#' x  <- rpois(n, lambda = 15)
#' y0 <- 0.1 * x + 0.002 * x^2 + rnorm(n)
#' y1 <- y0 + 1                      # treatment adds 1
#' cutoff <- 16
#' y <- ifelse(x >= cutoff, y1, y0)  # sharp assignment
#'
#' # Assume 90 % of post‑cutoff mass is genuine
#' tc <- data.frame(
#'   x      = cutoff:max(x),
#'   n_true = round(tabulate(x[x >= cutoff] + 1L) * 0.9)
#' )
#'
#' bounds_sharp(x, y, cutoff, true_counts = tc)
#'
#' @importFrom stats lm.wfit
#' @importFrom CVXR Variable Maximize Minimize Problem sum_entries solve
# -----------------------------------------------------------------------------
bounds_sharp <- function(x, y, cutoff,
                         true_counts,
                         weights    = NULL,
                         poly_order = 1L,
                         bounds     = c("both", "lower", "upper"),
                         solver     = getOption("rdpartial.solver", "ECOS"),
                         ...) {
  # ---- basic sanity checks --------------------------------------------------
  stopifnot(is.numeric(x), is.numeric(y), length(x) == length(y))
  .assert_scalar_numeric(cutoff, "cutoff")
  if (is.null(weights)) {
    weights <- rep(1, length(x))
  } else {
    if (any(is.na(weights))) weights[is.na(weights)] <- 0
    tiny_neg <- weights < 0 & weights > -sqrt(.Machine$double.eps)
    weights[tiny_neg] <- 0
    stopifnot(is.numeric(weights), length(weights) == length(x),
              all(weights >= 0))
  }
  bounds <- match.arg(bounds)

  .check_columns(true_counts, c("x", "n_true"))
  stopifnot(all(true_counts$x >= cutoff))

  # ---- split sample ---------------------------------------------------------
  left  <- x <  cutoff
  right <- x >= cutoff
  if (!any(left))  stop("No observations to the left of the cutoff.",  call. = FALSE)
  if (!any(right)) stop("No observations to the right of the cutoff.", call. = FALSE)

  # ---- design matrices & left‑hand fit -------------------------------------
  Xl <- .poly_design(x[left],  poly_order)
  Xr <- .poly_design(x[right], poly_order)
  c_star <- .poly_design(cutoff, poly_order)[1, ]

  Wl <- weights[left]
  Wr <- weights[right]

  fit_left <- stats::lm.wfit(x = Xl, y = y[left], w = Wl)
  left_at_cut <- drop(c_star %*% fit_left$coefficients)

  # ---- prepare manipulation constraints ------------------------------------
  xr_vals <- x[right]
  uniq_r  <- sort(unique(xr_vals))
  n_true  <- true_counts$n_true[match(uniq_r, true_counts$x)]
  if (anyNA(n_true)) {
    stop("`true_counts` missing some support points present in data.", call. = FALSE)
  }

  # H matrix: rows = support points, cols = obs.
  H <- matrix(0, nrow = length(uniq_r), ncol = length(xr_vals))
  for (j in seq_along(uniq_r)) H[j, xr_vals == uniq_r[j]] <- 1

  # Pre-compute phi and b as in paper ----------------------------------------
  phi_inv_cstar <- Xr %*% solve(t(Xr) %*% (Xr * Wr)) %*% c_star
  b_vec <- as.vector(y[right] * Wr * phi_inv_cstar)

  # CVXR variables -----------------------------------------------------------
  v <- CVXR::Variable(length(xr_vals))
  constr <- list(0 <= v, v <= 1, H %*% v == n_true)

  out <- c(lower = NA_real_, upper = NA_real_)

  cvx_opts <- modifyList(getOption("rdpartial.cvx_opts", list()), list(...))

  if (bounds != "upper") {
    obj_low <- CVXR::Minimize(sum(b_vec * v) - left_at_cut)
    sol_low <- CVXR::solve(CVXR::Problem(obj_low, constr), solver = solver, .opts = cvx_opts)
    out["lower"] <- sol_low$value
    attr(out, "opt_lwr") <- sol_low
  }
  if (bounds != "lower") {
    obj_up  <- CVXR::Maximize(sum(b_vec * v) - left_at_cut)
    sol_up  <- CVXR::solve(CVXR::Problem(obj_up, constr), solver = solver, .opts = cvx_opts)
    out["upper"] <- sol_up$value
    attr(out, "opt_upr") <- sol_up
  }

  if (bounds != "both") return(out[[bounds]])
  out
}
