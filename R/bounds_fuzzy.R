#' Partial‑Identification Bounds for **Fuzzy** RDDs with Manipulation
#'
#' Extension of [bounds_sharp()] to *fuzzy* designs where treatment take‑up
#' jumps (up **or** down) at the cutoff but is not deterministic.  The estimator
#' solves a linear‑fractional programme as described in Rosenman *et al.* (2025)
#' using `CVXR`.  Manipulation to the right of the cutoff is bounded by
#' `true_counts` obtained from [.density_estimation()].
#'
#' @section Required Inputs:
#' * **`x`** – running variable.
#' * **`y`** – outcome.
#' * **`z`** – realised treatment indicator (0/1).
#' * **`cutoff`** – threshold.
#' * **`true_counts`** – `data.frame(x, n_true)` of non‑manipulated counts.
#'
#' @param x Numeric running variable.
#' @param y Numeric outcome variable.
#' @param z Numeric vector (0/1) – treatment take‑up.
#' @param cutoff Numeric threshold separating control (< cutoff) from treatment (>= cutoff).
#' @param true_counts Data.frame with columns `x` (support points >= cutoff) and `n_true` (estimated number of non-manipulated observations at each support point).
#' @param weights Numeric vector of observation weights (optional).
#' @param poly_order Integer polynomial order for local regression (default 1L).
#' @param bounds Character – which bound(s) to return; one of "both" (default), "lower", "upper".
#' @param treat_direction Either "increase" (treatment prob. jumps *up* at the
#'   cutoff – standard case) or "decrease" (jumps *down*, e.g. donor *deferral*
#'   example).  Determines the sign of the optimisation objective.
#' @param solver Character – CVXR solver (default taken from `getOption("rdpartial.solver")`; package default is "ECOS").
#' @param runVarPlots Logical. If TRUE, generate running variable plots showing outcomes and treatment probabilities.
#' @param ylab Character. Label for the y-axis in plots.
#' @param xlab Character. Label for the x-axis in plots.
#' @param ... Further arguments passed on to [CVXR::solve()]. Use this to override CVXR tolerances or provide a different solver control list.
#'
#' @return Numeric vector of length 2 (`lower`, `upper`) when `bounds = "both"`;
#'   otherwise a single numeric.  CVXR solutions attached as attributes
#'   `opt_lwr`/`opt_upr`.
#'
#' @export
#' @examples
#' set.seed(321)
#' n      <- 6000
#' x      <- rpois(n, 15)
#' cutoff <- 16
#' z      <- rbinom(n, 1, prob = ifelse(x >= cutoff, 0.8, 0.2))
#' y0     <- 0.2 * x + rnorm(n)
#' y1     <- y0 + 2
#' y      <- ifelse(z == 1, y1, y0)
#'
#' freq <- table(factor(x[x >= cutoff], levels = cutoff:max(x)))
#' tc   <- data.frame(
#'   x      = as.integer(names(freq)),
#'   n_true = round(as.vector(freq) * 0.85)
#' )
#'
#' bounds_fuzzy(x, y, z, cutoff, true_counts = tc)
#'
#' @importFrom stats lm.wfit
#' @importFrom CVXR Variable Maximize Minimize Problem solve sum_entries
# -----------------------------------------------------------------------------
bounds_fuzzy <- function(x,
                         y,
                         z,
                         cutoff,
                         true_counts,
                         weights = NULL,
                         poly_order = 1L,
                         bounds = c("both", "lower", "upper"),
                         treat_direction = c("increase", "decrease"),
                         solver = getOption("rdpartial.solver", "ECOS"),
                         runVarPlots = FALSE,
                         ylab = NULL,
                         xlab = NULL,
                         ...) {
  # ---- checks ---------------------------------------------------------------
  stopifnot(
    is.numeric(x),
    is.numeric(y),
    is.numeric(z),
    length(x) == length(y),
    length(x) == length(z)
  )
  .assert_scalar_numeric(cutoff, "cutoff")

  if (is.null(weights)) {
    weights <- rep(1, length(x))
  } else {
    if (any(is.na(weights)))
      weights[is.na(weights)] <- 0
    tiny_neg <- weights < 0 & weights > -sqrt(.Machine$double.eps)
    weights[tiny_neg] <- 0
    stopifnot(is.numeric(weights),
              length(weights) == length(x),
              all(weights >= 0))
  }

  bounds          <- match.arg(bounds)
  treat_direction <- match.arg(treat_direction)

  .check_columns(true_counts, c("x", "n_true"))
  stopifnot(all(true_counts$x >= cutoff))

  if(treat_direction == 'decrease')
    z <- 1 - z

  # ---- split sample ---------------------------------------------------------
  left  <- x <  cutoff
  right <- x >= cutoff
  if (!any(left))
    stop("No observations to the left of the cutoff.", call. = FALSE)
  if (!any(right))
    stop("No observations to the right of the cutoff.", call. = FALSE)

  # ---- design matrices ------------------------------------------------------
  Xl <- .poly_design(x[left], poly_order)
  Xr <- .poly_design(x[right], poly_order)
  c_star <- .poly_design(cutoff, poly_order)[1, ]

  Wl <- weights[left]
  Wr <- weights[right]

  # ---- left-side fits (numerator & denominator) -----------------------------
  fit_left_y <- stats::lm.wfit(x = Xl, y = y[left], w = Wl)
  fit_left_z <- stats::lm.wfit(x = Xl, y = z[left], w = Wl)
  left_num   <- drop(c_star %*% fit_left_y$coefficients)
  left_den   <- drop(c_star %*% fit_left_z$coefficients)

  # ---- manipulation constraints --------------------------------------------
  xr_vals <- x[right]
  uniq_r  <- sort(unique(xr_vals))
  n_true  <- true_counts$n_true[match(uniq_r, true_counts$x)]
  if (anyNA(n_true))
    stop("`true_counts` missing some support points present in data.",
         call. = FALSE)

  wStar <- rep(1, length(xr_vals))
  H <- matrix(0, nrow = length(uniq_r), ncol = length(xr_vals))
  for (j in seq_along(uniq_r)) {
    H[j, xr_vals == uniq_r[j]] <- 1
    wStar[xr_vals == uniq_r[j]] <- c(rep(1, n_true[j]),
                                     rep(0, sum(xr_vals == uniq_r[j]) - n_true[j]))
  }


  # ---- phi terms ------------------------------------------------------------
  phi_inv_c <- Xr %*% solve(t(Xr) %*% (Xr * (Wr * wStar))) %*% c_star
  b_num   <- as.vector(y[right] * Wr * phi_inv_c)
  b_den   <- as.vector(z[right] * Wr * phi_inv_c)

  # ---- CVXR variables -------------------------------------------------------
  y_var <- CVXR::Variable(length(xr_vals))  # weights per obs
  z_var <- CVXR::Variable()                 # scaling factor

  if (treat_direction == "increase") {
    obj_up   <- CVXR::Maximize(sum(b_num * y_var) - left_num * z_var)
    obj_low  <- CVXR::Minimize(sum(b_num * y_var) - left_num * z_var)
    denom_eq <- sum(b_den * y_var) - left_den * z_var == 1
  } else {
    # "decrease" – sign flips
    obj_up   <- CVXR::Maximize(-sum(b_num * y_var) + left_num * z_var)
    obj_low  <- CVXR::Minimize(-sum(b_num * y_var) + left_num * z_var)
    denom_eq <- -sum(b_den * y_var) + left_den * z_var == 1
  }

  constr <- list(z_var >= 0,
                 H %*% y_var == n_true * z_var,
                 denom_eq,
                 0 <= y_var,
                 y_var <= z_var)

  out <- c(lower = NA_real_, upper = NA_real_)
  cvx_opts <- modifyList(getOption("rdpartial.cvx_opts", list()), list(...))

  # ---- optimisation ---------------------------------------------------------
  if (bounds != "upper") {
    sol_low <- CVXR::solve(CVXR::Problem(obj_low, constr),
                           solver = solver,
                           .opts = cvx_opts)
    out["lower"] <- sol_low$value
    attr(out, "opt_lwr") <- sol_low
  }
  if (bounds != "lower") {
    sol_up  <- CVXR::solve(CVXR::Problem(obj_up, constr),
                           solver = solver,
                           .opts = cvx_opts)
    out["upper"] <- sol_up$value
    attr(out, "opt_upr") <- sol_up
  }

  # make the running variable plot
  if(runVarPlots) {

    # get the weights
    upperWeights <- sol_up$getValue(y_var) / sol_up$getValue(z_var)
    lowerWeights <- sol_low$getValue(y_var) / sol_low$getValue(z_var)

    # get the proportion at each value of Xr
    v_y <- tapply(c(y[left], y[right]), c(x[left], x[right]), mean)
    prop.y <- data.frame(
      x_value  = as.numeric(names(v_y)),
      avg.prop = as.numeric(v_y)
    )

    v_z <- tapply(c(z[left], z[right]), c(x[left], x[right]), mean)
    prop.z <- data.frame(
      x_value  = as.numeric(names(v_z)),
      avg.prop = as.numeric(v_z)
    )

    # set appropriate names
    myHist <- data.frame(table(x))
    names(myHist)[1] <- 'x_value'
    names(true_counts) <- c("x_value", "n_true")

    # generate the outcomes plot
    yPlot <- outcomes_plot(
      prop          = prop.y,
      xl            = x[left],
      Yl            = y[left],
      xr            = x[right],
      Yr            = y[right],
      upperWeights  = upperWeights, # weights for right‑side upper bound LOESS
      lowerWeights  = lowerWeights, # weights for right‑side lower bound LOESS
      ylab          = ylab,
      xlab          = xlab,
      title         = "Outcome vs. Running Variable",
      order         = poly_order,   # match the polynomial order used earlier
      hist          = myHist,         # histogram of running variable values
      trueCounts    = true_counts,  # data.frame(x, n_true) from the inputs
      cutoff        = cutoff        # the actual cutoff threshold
    )

    # generate the treatments plot
    zPlot <- outcomes_plot(
      prop          = prop.z,
      xl            = x[left],
      Yl            = z[left],
      xr            = x[right],
      Yr            = z[right],
      upperWeights  = upperWeights, # weights for right‑side upper bound LOESS
      lowerWeights  = lowerWeights, # weights for right‑side lower bound LOESS
      ylab          = ylab,
      xlab          = "Treatment Indicator",
      title         = "Treatment vs. Running Variable",
      order         = poly_order,   # match the polynomial order used earlier
      hist          = myHist,         # histogram of running variable values
      trueCounts    = true_counts,  # data.frame(x, n_true) from the inputs
      cutoff        = cutoff        # the actual cutoff threshold
    )
  }

  # return bounds
  if(!runVarPlots) {
    if (bounds != "both")
      return(out[[bounds]])
    out
  } else {
    return(list(bounds = out,
                yPlot = yPlot,
                zPlot = zPlot))
  }
}
