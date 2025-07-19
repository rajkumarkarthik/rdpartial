post_counts_df <- function(x_vec, cutoff, keep_frac = 1) {
  tbl <- as.data.frame(table(x_vec[x_vec >= cutoff]), stringsAsFactors = FALSE)
  names(tbl) <- c("x", "n_true")
  tbl$x <- as.numeric(tbl$x)
  tbl$n_true <- round(tbl$n_true * keep_frac)
  tbl
}

cutoff <- 16

# 1) Finite bounds & ordering ----------------------------------------------

test_that("bounds_fuzzy returns finite, ordered bounds (increase)", {
  dat <- test_data
  tc  <- post_counts_df(dat$x, cutoff, keep_frac = 0.9)

  b_inc <- bounds_fuzzy(dat$x, dat$y, dat$z, cutoff,
                        true_counts     = tc,
                        poly_order      = 1,
                        treat_direction = "increase")
  b_inc_num <- as.numeric(b_inc)

  expect_true(all(is.finite(b_inc_num)))
  expect_length(b_inc_num, 2)
  expect_lt(b_inc_num[1], b_inc_num[2])
})

# 2) Single‑bound API -------------------------------------------------------

test_that("bounds_fuzzy returns single numeric when bounds = 'upper'", {
  dat <- test_data
  tc  <- post_counts_df(dat$x, cutoff, keep_frac = 0.9)

  upr <- bounds_fuzzy(dat$x, dat$y, dat$z, cutoff,
                      true_counts   = tc,
                      poly_order    = 1,
                      bounds        = "upper",
                      treat_direction = "increase")
  expect_true(is.numeric(upr) && length(upr) == 1)
})
