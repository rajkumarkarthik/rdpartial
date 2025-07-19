library(testthat)

test_that(".density_estimation() basic invariants", {
  hist_df <- as.data.frame(table(test_data$x))
  names(hist_df) <- c("x", "freq")
  hist_df$x <- as.numeric(as.character(hist_df$x))

  res <- rdpartial:::.density_estimation(
    hist_df      = hist_df,
    manip_region = c(15.7, 16),
    cutoff       = 16,
    make_plot    = FALSE,
    num_folds    = 3
  )

  expect_type(res, "list")
  expect_named(res, c("counts", "avg_cv_sse", "avg_sse", "knots", "model", "plot"),
               ignore.order = TRUE)

  # Counts are non-negative integers ---------------------------------------
  expect_true(all(res$counts$n_true >= 0))
  expect_true(all(res$counts$n_true == round(res$counts$n_true)))

  # Predicted non‑manipulated counts should never exceed observed counts ----
  joined <- merge(res$counts, hist_df, by = "x", all.x = TRUE)
  expect_true(all(joined$n_true <= joined$freq))
})
