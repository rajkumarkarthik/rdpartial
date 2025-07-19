# 1) Collapse test: when the estimated non‑manipulated counts equal the
#    *observed* post‑cutoff histogram, the upper and lower Manski bounds should
#    coincide (width ≈ 0).

test_that("bounds_sharp collapses when no manipulation", {
  dat    <- test_data
  cutoff <- 16

  # Build exact post‑cutoff histogram in the shape expected by bounds_sharp()
  post_x <- dat$x[dat$x >= cutoff]
  tc <- as.data.frame(table(post_x))
  names(tc) <- c("x", "n_true")
  tc$x <- as.numeric(as.character(tc$x))  # ensure numeric

  b <- bounds_sharp(dat$x, dat$y,
                    cutoff      = cutoff,
                    true_counts = tc,
                    poly_order  = 1)

  expect_length(b, 2)
  expect_lt(diff(b), 1e-6)
})

# 2) Basic argument validation ------------------------------------------------

test_that("bounds_sharp argument validation", {
  expect_error(bounds_sharp(1:3, 1:2, 1, data.frame(x = 1, n_true = 1)),
               "length")
})
