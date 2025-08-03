#' Plot Running Variable Outcomes with LOESS Fits
#'
#' Generates a ggplot visualization of outcomes across hemoglobin levels,
#' including standard LOESS fits and upper/lower bounds.  Intended for use with
#' partial‑identification analyses to illustrate treatment and outcome
#' probabilities around the cutoff.
#'
#' @section Required Inputs:
#' * **`prop`** – data.frame with columns `h.level` and `avg.prop`.
#' * **`hist`** – data.frame with histogram counts (`hlevel`, `Freq`).
#' * **`true_counts`** – data.frame with columns `hlevel`, `numTrueSubjects`.
#' * **`xl`, `Yl`** – numeric vectors for left‑side LOESS regression.
#' * **`xr`, `Yr`** – numeric vectors for right‑side LOESS regression.
#'
#' @param prop Data.frame with columns `h.level`, `avg.prop`.
#' @param xl Numeric vector of x-values for left-side LOESS.
#' @param Yl Numeric vector of y-values for left-side LOESS.
#' @param xr Numeric vector of x-values for right-side LOESS.
#' @param Yr Numeric vector of y-values for right-side LOESS.
#' @param upperWeights Numeric vector of weights for the upper bound LOESS.
#' @param lowerWeights Numeric vector of weights for the lower bound LOESS.
#' @param ylab Character. Label for the y-axis.
#' @param title Character. Plot title.
#' @param order Integer. Degree of LOESS polynomial (1 = linear, 2 = quadratic).
#' @param hist Data.frame with histogram counts (`hlevel`, `Freq`).
#' @param trueCounts Data.frame with non-manipulated subject counts (`hlevel`,
#'   `numTrueSubjects`).
#' @param cutoff Numeric. Hemoglobin cutoff value.
#'
#' @return A \code{ggplot} object representing the outcomes plot.
#'
#' @examples
#' \dontrun{
#' outcomes_plot(prop_df, xl, Yl, xr, Yr,
#'               upperWeights, lowerWeights,
#'               ylab = "Treatment Probability",
#'               title = "Treatment by Hemoglobin Level",
#'               order = 1, hist_df, counts_df, cutoff = 12)
#' }
#'
#' @import ggplot2
#' @importFrom dplyr coalesce
#' @export
# -----------------------------------------------------------------------------
outcomes_plot <- function(prop,
                          xl,
                          Yl,
                          xr,
                          Yr,
                          upperWeights,
                          lowerWeights,
                          ylab,
                          xlab,
                          title,
                          order,
                          hist,
                          trueCounts,
                          cutoff) {

  # ---- checks ---------------------------------------------------------------
  .check_columns(prop,       c("h.level", "avg.prop"))
  .check_columns(hist,       c("hlevel", "Freq"))
  .check_columns(trueCounts, c("hlevel", "numTrueSubjects"))
  .assert_scalar_numeric(cutoff, "cutoff")

  # ---- merge and prepare data -----------------------------------------------
  h <- merge(hist, trueCounts, 'hlevel', all.x = TRUE)
  prop <- merge(prop, h, by.x = 'h.level', by.y = 'hlevel')
  prop$numTrueSubjects <- dplyr::coalesce(prop$numTrueSubjects, prop$Freq)
  prop$lr <- ifelse(prop$h.level < cutoff, 'l', 'r')

  # ---- LOESS fits -----------------------------------------------------------
  prop$leftLine <- c(
    predict(stats::loess(Yl ~ xl, span = 1, degree = order,
                         control = stats::loess.control(surface = "direct")),
            newdata = c(sort(unique(xl)), cutoff)),
    rep(NA, sum(prop$lr == 'r') - 1)
  )

  prop$rightLine <- c(
    rep(NA, sum(prop$lr == 'l')),
    predict(stats::loess(Yr ~ xr, span = 1, degree = order),
            newdata = sort(unique(xr)))
  )

  prop$right.upr <- c(
    rep(NA, sum(prop$lr == 'l')),
    predict(stats::loess(Yr ~ xr, span = 1,
                         weights = round(upperWeights, 3),
                         degree = order),
            newdata = sort(unique(xr)))
  )

  prop$right.lwr <- c(
    rep(NA, sum(prop$lr == 'l')),
    predict(stats::loess(Yr ~ xr, span = 1,
                         weights = lowerWeights,
                         degree = order),
            newdata = sort(unique(xr)))
  )

  # ---- build plot -----------------------------------------------------------
  browser()
  m <- ggplot2::ggplot(prop, ggplot2::aes(x = h.level, y = avg.prop)) +
    ggplot2::geom_point(ggplot2::aes(size = numTrueSubjects),
                        colour = 'black', fill = 'gray', shape = 21) +
    ggplot2::theme_classic() +
    ggplot2::theme(text = ggplot2::element_text(size = 16)) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
    ggplot2::labs(size = '# of Non-Manipulators') +
    ggplot2::labs(y = ylab, x = xlab) +
    ggplot2::geom_vline(xintercept = cutoff, lty = 2, col = 'gray') +
    ggplot2::geom_line(ggplot2::aes(x = h.level, y = leftLine,
                                    color = 'Standard Estimate'), lwd = 1.5) +
    ggplot2::geom_line(ggplot2::aes(x = h.level, y = rightLine,
                                    color = 'Standard Estimate'), lwd = 1.5) +
    ggplot2::geom_line(ggplot2::aes(x = h.level, y = right.upr,
                                    color = 'Upper Bound'), lwd = 1.5, lty = 2) +
    ggplot2::geom_line(ggplot2::aes(x = h.level, y = right.lwr,
                                    color = 'Lower Bound'), lwd = 1.5, lty = 2) +
    ggplot2::scale_colour_manual("LOESS Fit",
                                 breaks = c("Standard Estimate", "Upper Bound", "Lower Bound"),
                                 values = c("blue", "green", "red")) +
    ggplot2::ggtitle(title)

  return(m)
}
