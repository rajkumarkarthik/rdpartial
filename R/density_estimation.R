#' Estimate the number of *non-manipulated* observations to the right of the cutoff
#'
#' This helper replicates the constrained Poisson-spline procedure from the
#' research code (Rosenman *et al.*, 2025) but in a package-friendly, testable
#' format.  It is **not exported** - downstream functions access it internally
#' via `.density_estimation()`.
#'
#' @details
#' The algorithm treats the running-variable histogram as the realisation of a
#' Poisson process.  A flexible B-spline (possibly augmented with indicator
#' spikes at integer or half-integer values) is fit to the *donut* region -
#' bins outside the user-supplied manipulation interval.  Mass preservation and
#' continuity constraints are imposed through `CVXR`, and the optimal number of
#' knots / spike structure is chosen by $K$-fold cross-validation.  The final
#' step predicts the expected counts in the manipulation interval, which serve
#' as the "true" (non-manipulated) totals fed into the partial-bounds
#' estimators.  Bins with zero counts receive a tiny offset before logging so
#' that the inequality constraints remain finite.
#'
#' @param hist_df `data.frame` with columns `x` (numeric support points)
#'   and `freq` (integer counts).
#' @param manip_region Numeric length-2 vector giving the suspected manipulation
#'   interval *(lower, upper)* **inclusive**.  Must contain the cutoff on its
#'   **upper** bound.
#' @param cutoff Numeric scalar - RDD threshold.
#' @param knot_options Integer vector of candidate spline knot counts;
#'   default `3:15`.
#' @param mod_types Character vector specifying which model structures to
#'   cross-validate over.  Allowed values are:
#'   * "smooth" - plain spline;
#'   * "spike_integer" - spline + indicator for integer spikes;
#'   * "spike_half" - spline + integer + half-integer spikes.
#' @param lambda Penalty weight on total predicted mass inside the manipulation
#'   region (higher ⇒ stricter mass match).  Default `100` (paper default).
#' @param eps Numeric slack for the CVXR continuity constraints (default `1e-6`).
#' @param folds Optional integer vector assigning each *donut* bin to a CV fold.
#'   Length must equal `nrow(hist_df)`.  When `NULL`, a random `num_folds`-way
#'   split is generated.
#' @param num_folds Integer number of folds if `folds` is `NULL` (default `5`).
#' @param make_plot Logical; if `TRUE` return a `ggplot2` object visualising the
#'   final fit.
#' @param seed Optional integer for reproducible fold assignment.
#' @param solver CVXR solver to use.  Defaults to
#'   `.get_rdpartial_opt("solver")`.
#' @param cvx_opts List of solver options passed to `CVXR::solve()`.  Defaults to
#'   `.get_rdpartial_opt("cvx_opts")`.
#'
#' @return A list with components:
#' * `counts`       - `data.frame(x, n_true)` within `manip_region`.
#' * `avg_cv_sse`   - cross-validated SSE (per bin) for the chosen model.
#' * `avg_sse`      - in-sample SSE outside `manip_region`.
#' * `knots`        - optimal knot count.
#' * `model`        - optimal model type (one of `mod_types`).
#' * `plot`         - `ggplot` object when `make_plot = TRUE`.
#'
#' @keywords internal
#' @importFrom stats glm poisson constrOptim quantile
#' @importFrom splines bs
#' @importFrom CVXR Variable Maximize Problem solve
#' @importFrom ggplot2 ggplot aes geom_bar geom_line scale_fill_manual theme
#'   theme_classic geom_vline labs xlim ylim position_nudge
.density_estimation <- function(hist_df, manip_region, cutoff,
                                knot_options = 3:15,
                                mod_types    = c("smooth", "spike_integer"),
                                lambda       = 100, eps = 1e-6,
                                folds = NULL, num_folds = 5L,
                                make_plot = FALSE, seed = NULL,
                                solver = .get_rdpartial_opt("solver"),
                                cvx_opts = .get_rdpartial_opt("cvx_opts")) {
  # ---- sanity checks -------------------------------------------------------
  .check_columns(hist_df, c("x", "freq"))
  stopifnot(is.numeric(hist_df$x), is.numeric(hist_df$freq))
  if (any(hist_df$freq < 0)) stop("`freq` must be non-negative.", call. = FALSE)

  if (!is.null(seed)) set.seed(seed)
  .assert_scalar_numeric(cutoff, "cutoff")
  stopifnot(is.numeric(manip_region), length(manip_region) == 2, manip_region[1] < manip_region[2])

  # Ensure the manipulation interval indeed covers the cutoff ----------------
  if (cutoff < manip_region[1] || cutoff > manip_region[2]) {
    stop("`manip_region` must contain the `cutoff`.", call. = FALSE)
  }

  # Basic indices ------------------------------------------------------------
  in_manip   <- with(hist_df, x >= manip_region[1] & x <= manip_region[2])
  left_sub   <- with(hist_df, x >= manip_region[1] & x <  cutoff)
  right_sub  <- with(hist_df, x >= cutoff          & x <= manip_region[2])
  hist_donut <- hist_df[!in_manip, , drop = FALSE]
  mid_vol    <- sum(hist_df$freq[in_manip])
  # Small offset prevents log(0) when bins have zero counts
  offset     <- .Machine$double.eps
  ci_vec     <- c(log(hist_df$freq[left_sub]  + offset),
                  -log(hist_df$freq[right_sub] + offset))

  # CV folds -----------------------------------------------------------------
  if (is.null(folds)) {
    folds <- sample(rep(seq_len(num_folds), length.out = nrow(hist_donut)))
  } else {
    if (length(folds) != nrow(hist_df)) {
      stop("`folds` must have length equal to `nrow(hist_df)`.", call. = FALSE)
    }
    folds <- folds[!in_manip]
  }

  # Helper: Poisson + mass penalty objective ---------------------------------
  obj_poisson <- function(beta, yy, mm, mm_donut, mid_vol, lambda_val) {
    -(sum(yy * (mm %*% beta) - exp(mm %*% beta))) +
      lambda_val * (sum(exp(mm_donut %*% beta)) - mid_vol)^2
  }

  # Grid search over knot × model -------------------------------------------
  sse_mat <- matrix(Inf, nrow = length(knot_options), ncol = length(mod_types))

  for (k_idx in seq_along(knot_options)) {
    k <- knot_options[k_idx]

    for (m_idx in seq_along(mod_types)) {
      m_type <- tolower(mod_types[m_idx])

      cv_sse <- numeric(num_folds)
      for (fold in seq_len(num_folds)) {
        # Build full model matrix ----------------------------------------
        mm_full <- switch(m_type,
                          smooth        = model.matrix(~ splines::bs(x, k),              data = hist_df),
                          spike_integer = model.matrix(~ splines::bs(x, k) + I(round(x) == x), data = hist_df),
                          spike_half    = model.matrix(~ splines::bs(x, k) + I(round(x) == x) +
                                                         I((10 * x) %% 10 == 5),         data = hist_df),
                          stop("Unknown model type: ", m_type, call. = FALSE)
        )

        mm_donut <- mm_full[in_manip, , drop = FALSE]
        mm_clean <- mm_full[!in_manip, , drop = FALSE]
        ui_mat   <- rbind(mm_full[left_sub,  , drop = FALSE],
                          -mm_full[right_sub, , drop = FALSE])

        # Fit Poisson GLM on donut bins -----------------------------------
        glm_fit <- stats::glm(freq ~ . - 1, family = poisson(),
                              data = data.frame(freq = hist_donut$freq, mm_clean),
                              subset = folds != fold)

        beta_var <- CVXR::Variable(ncol(mm_full))
        obj_cvx  <- CVXR::Maximize(sum(glm_fit$y * (model.matrix(glm_fit) %*% beta_var) -
                                         exp(model.matrix(glm_fit) %*% beta_var)))
        prob_cvx <- CVXR::Problem(obj_cvx, list(ui_mat %*% beta_var - ci_vec >= eps))
        sol      <- CVXR::solve(prob_cvx, solver = solver, .opts = cvx_opts)
        if (sol$status == "solver_error") {
          cv_sse[fold] <- Inf
          next
        }

        theta0 <- as.numeric(sol$getValue(beta_var))
        # Compute per-constraint slack and build an interior RHS
        slack  <- as.vector(ui_mat %*% theta0 - ci_vec)
        min_sl <- min(slack)
        # If the smallest slack ≤ 0, shift *all* constraints inward so that
        # the CVXR point is strictly interior; otherwise just give a tiny buffer.
        if (min_sl <= 0) {
          ci_int <- ci_vec + min_sl - 1e-8
        } else {
          ci_int <- ci_vec - 1e-8
        }

        co_res <- stats::constrOptim(theta = as.numeric(sol$getValue(beta_var)),
                                     f     = obj_poisson,
                                     grad  = NULL,
                                     ui    = ui_mat, ci = ci_int,
                                     yy    = glm_fit$y, mm = model.matrix(glm_fit),
                                     mm_donut = mm_donut, mid_vol = mid_vol,
                                     lambda_val = lambda)

        preds <- exp(mm_clean[folds == fold, , drop = FALSE] %*% co_res$par)
        cv_sse[fold] <- sum((preds - hist_donut$freq[folds == fold])^2)
      }
      sse_mat[k_idx, m_idx] <- sum(cv_sse)
    }
  }

  best_idx  <- arrayInd(which.min(sse_mat), dim(sse_mat))
  best_k    <- knot_options[best_idx[1]]
  best_type <- tolower(mod_types[best_idx[2]])

  # ---- Refit on full dataset using best spec -------------------------------
  mm_full <- switch(best_type,
                    smooth        = model.matrix(~ splines::bs(x, best_k), data = hist_df),
                    spike_integer = model.matrix(~ splines::bs(x, best_k) + I(round(x) == x), data = hist_df),
                    spike_half    = model.matrix(~ splines::bs(x, best_k) + I(round(x) == x) +
                                                   I((10 * x) %% 10 == 5), data = hist_df)
  )

  mm_donut <- mm_full[in_manip, , drop = FALSE]
  ui_mat   <- rbind(mm_full[left_sub,  , drop = FALSE],
                    -mm_full[right_sub, , drop = FALSE])

  glm_fit <- stats::glm(freq ~ . - 1, family = poisson(),
                        data = data.frame(freq = hist_donut$freq, mm_full[!in_manip, , drop = FALSE]))

  beta_var <- CVXR::Variable(ncol(mm_full))
  obj_cvx  <- CVXR::Maximize(sum(glm_fit$y * (model.matrix(glm_fit) %*% beta_var) -
                                   exp(model.matrix(glm_fit) %*% beta_var)))
  prob_cvx <- CVXR::Problem(obj_cvx, list(ui_mat %*% beta_var - ci_vec >= eps))
  sol      <- CVXR::solve(prob_cvx, solver = solver, .opts = cvx_opts)

  theta0 <- as.numeric(sol$getValue(beta_var))
  # Compute per-constraint slack and build an interior RHS
  slack  <- as.vector(ui_mat %*% theta0 - ci_vec)
  min_sl <- min(slack)
  # If the smallest slack ≤ 0, shift *all* constraints inward so that
  # the CVXR point is strictly interior; otherwise just give a tiny buffer.
  if (min_sl <= 0) {
    ci_int <- ci_vec + min_sl - 1e-8
  } else {
    ci_int <- ci_vec - 1e-8
  }

  co_res <- stats::constrOptim(theta = as.numeric(sol$getValue(beta_var)),
                               f     = obj_poisson, grad = NULL,
                               ui = ui_mat, ci = ci_int,
                               yy    = glm_fit$y, mm = model.matrix(glm_fit),
                               mm_donut = mm_donut, mid_vol = mid_vol,
                               lambda_val = lambda)
  preds <- exp(mm_full %*% co_res$par)

  # Output -------------------------------------------------------------------
  right_side <- hist_df$x >= cutoff
  n_true_vec <- hist_df$freq                       # start with observed counts
  n_true_vec[in_manip & right_side] <- preds[in_manip & right_side]

  counts_df <- data.frame(
    x      = hist_df$x[right_side],
    n_true = round(n_true_vec[right_side])
  )

  plot_out <- NULL
  if (make_plot) {
    tmp <- hist_df
    tmp$manip <- ifelse(in_manip, "Yes", "No")
    tmp$pred  <- preds

    plot_out <- ggplot2::ggplot(tmp, ggplot2::aes(x = x, y = freq, fill = manip)) +
      ggplot2::geom_bar(stat = "identity", width = 0.1, colour = "darkgray",
                        position = ggplot2::position_nudge(x = +0.05)) +
      ggplot2::geom_line(ggplot2::aes(y = pred, fill = NULL),
                         colour = "red", linewidth = 0.8,
                         position = ggplot2::position_nudge(x = +0.05)) +
      ggplot2::scale_fill_manual("Manipulation\nRegion", values = c(Yes = "gray", No = "black")) +
      ggplot2::geom_vline(xintercept = cutoff, linetype = 2) +
      ggplot2::theme_classic() +
      ggplot2::labs(x = "Running variable", y = "Count") +
      ggplot2::xlim(range(hist_df$x))
  }

  list(
    counts     = counts_df,
    avg_cv_sse = min(sse_mat) / nrow(hist_donut),
    avg_sse    = mean((preds[!in_manip] - hist_df$freq[!in_manip])^2),
    knots      = best_k,
    model      = best_type,
    plot       = plot_out
  )
}
