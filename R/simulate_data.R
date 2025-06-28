#' Simulate Synthetic Data for Sharp or Fuzzy RDD Examples
#'
#' A convenience generator that mirrors the data‐generating processes used in
#' the **rdpartial** documentation and unit tests.  The function produces a
#' data frame ready for `bounds_sharp()`, `bounds_fuzzy()`, and
#' `bootstrap_bounds()`, while also returning the *true* treatment effect at the
#' cutoff for quick sanity checks.
#'
#' @section Outline of the DGP:
#' 1. Draw an **integer** running variable `x` from either a Poisson or discrete
#'    uniform distribution.
#' 2. (Optional) Manipulate units with `x` in `(cutoff − manip_width, cutoff)`
#'    above the threshold with probability `manip_prob`.
#' 3. Generate potential outcomes `Y(0)` and `Y(1)` as quadratic functions of
#'    `x` plus Normal noise.
#' 4. Assign treatment either *sharply* (`z = 1(x ≥ cutoff)`) or *fuzzily* via
#'    Bernoulli probabilities `p0`, `p1` below / above the cutoff.
#' 5. Observe realised outcome `y = Y(z)`.
#'
#' @param n Integer sample size (≥ 1).
#' @param cutoff Numeric threshold.  Must be integer‐aligned with the support of
#'   `x`.
#' @param dist Either "poisson" (default) or "uniform".
#' @param lambda Poisson mean when `dist = "poisson"`.
#' @param range Integer length‑2 support when `dist = "uniform"`.
#' @param beta Numeric length‑3 vector for the quadratic mean of `Y(0)`.
#' @param tau Constant treatment effect added to `Y(1)` (default `1`).
#' @param sigma Error SD for both potential outcomes (default `1`).
#' @param design "fuzzy" (default) or "sharp".
#' @param p0,p1 Treatment probabilities below / above cutoff when
#'   `design = "fuzzy"`.
#' @param manip_width Numeric ≥ 0 — interval width subject to manipulation just
#'   below the cutoff.
#' @param manip_prob Probability of manipulation for units in that interval
#'   (default `0`).
#' @param seed Optional integer for reproducibility.
#'
#' @return A `data.frame` with columns:
#' * `x` – running variable **after** manipulation.
#' * `y` – realised outcome.
#' * `z` – treatment status.
#' * `x_orig` – running variable before manipulation.
#' * `manipulated` – logical flag.
#'
#' The object carries attributes:
#' * `ate_cutoff` – true treatment effect at the cutoff (`tau`).
#' * `call` – original function call.
#'
#' @export
#' @examples
#' set.seed(42)
#' sim <- simulate_rdd_data(n = 2000, cutoff = 16, design = "fuzzy",
#'                          manip_width = 0.4, manip_prob = 0.25)
#' head(sim)
#' attr(sim, "ate_cutoff")
#'
simulate_rdd_data <- function(n = 4000, cutoff = 16,
                              dist = c("poisson", "uniform"),
                              lambda = 15, range = c(8L, 25L),
                              beta = c(0, 0.1, 0.002), tau = 1, sigma = 1,
                              design = c("fuzzy", "sharp"), p0 = 0.2, p1 = 0.8,
                              manip_width = 0, manip_prob = 0, seed = NULL) {
  # ---- argument checks ------------------------------------------------------
  .assert_scalar_numeric(n, "n"); n <- as.integer(n)
  .assert_scalar_numeric(cutoff, "cutoff")
  if (!is.null(seed)) set.seed(seed)

  dist   <- match.arg(dist)
  design <- match.arg(design)
  if (length(beta) != 3) stop("`beta` must be length 3.", call. = FALSE)
  if (manip_width < 0 || manip_prob < 0 || manip_prob > 1)
    stop("`manip_width` must be ≥0 and `manip_prob` in [0,1].", call. = FALSE)

  # ---- running variable -----------------------------------------------------
  x_raw <- switch(dist,
                  poisson = stats::rpois(n, lambda = lambda),
                  uniform = sample(seq.int(range[1], range[2]), n, replace = TRUE)
  )

  # ---- optional manipulation ------------------------------------------------
  to_manip <- x_raw >= (cutoff - manip_width) & x_raw < cutoff &
    stats::runif(n) < manip_prob
  x_manip  <- ifelse(to_manip, cutoff + (x_raw - (cutoff - manip_width)), x_raw)

  # ---- treatment assignment -------------------------------------------------
  z <- if (design == "sharp") {
    as.integer(x_manip >= cutoff)
  } else {
    prob <- ifelse(x_manip >= cutoff, p1, p0)
    stats::rbinom(n, 1, prob)
  }

  # ---- outcomes -------------------------------------------------------------
  mu0 <- beta[1] + beta[2] * x_manip + beta[3] * x_manip^2
  mu1 <- mu0 + tau
  y   <- stats::rnorm(n, mean = ifelse(z == 1, mu1, mu0), sd = sigma)

  out <- data.frame(x = x_manip,
                    y = y,
                    z = z,
                    x_orig = x_raw,
                    manipulated = to_manip)

  attr(out, "ate_cutoff") <- tau
  attr(out, "call")       <- match.call()
  class(out) <- c("rdpartial_sim", "data.frame")
  out
}

# ---- S3 print method --------------------------------------------------------
#' @export
print.rdpartial_sim <- function(x, ...) {
  NextMethod()
  cat("\nTrue ATE at cutoff:", attr(x, "ate_cutoff"), "\n")
  invisible(x)
}
