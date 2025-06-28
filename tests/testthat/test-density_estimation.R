test_that("density_estimation handles zero counts", {
  hist_df <- data.frame(
    hlevel = 0:4,
    freq   = c(10, 0, 0, 5, 3)
  )

  expect_error(
    res <- rdpartial:::.density_estimation(
      hist_df,
      manip_region = c(1, 3),
      cutoff       = 2,
      knot_options = 3,
      mod_types    = "smooth",
      make_plot    = FALSE,
      seed         = 1
    ),
    regexp = NA
  )

  expect_true(is.list(res))
  expect_true(all(c("counts", "avg_cv_sse", "avg_sse", "knots", "model", "plot") %in% names(res)))
})
