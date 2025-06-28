library(testthat)
library(rdpartial)

# Small synthetic data for testing

test_that("bounds_sharp runs without dimension errors", {
  set.seed(1)
  x <- rep(14:17, each = 2)
  y0 <- 0.1 * x + rnorm(length(x))
  y1 <- y0 + 1
  cutoff <- 16
  y <- ifelse(x >= cutoff, y1, y0)
  tc <- data.frame(
    x = cutoff:max(x),
    n_true = tabulate(x[x >= cutoff] - cutoff + 1)
  )
  expect_error(bounds_sharp(x, y, cutoff, true_counts = tc), NA)
})


test_that("bounds_fuzzy runs without dimension errors", {
  set.seed(2)
  n <- 8
  x <- rep(14:17, length.out = n)
  z <- rbinom(n, 1, ifelse(x >= 16, 0.8, 0.2))
  y0 <- rnorm(n)
  y1 <- y0 + 1
  y <- ifelse(z == 1, y1, y0)
  tc <- data.frame(
    x = 16:max(x),
    n_true = tabulate(x[x >= 16] - 16 + 1)
  )
  expect_error(bounds_fuzzy(x, y, z, 16, true_counts = tc), NA)
})
