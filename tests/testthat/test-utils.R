library(testthat)

test_that(".assert_scalar_numeric errors and successes", {
  expect_error(.assert_scalar_numeric("a"), "single, non")
  expect_error(.assert_scalar_numeric(c(1, 2)), "single, non")
  expect_error(.assert_scalar_numeric(NA_real_), "single, non")
  expect_silent(.assert_scalar_numeric(2))
})

test_that(".poly_design matrix structure", {
  m0 <- .poly_design(1:3, 0)
  expect_equal(m0, matrix(1, 3, 1))
  m2 <- .poly_design(2, 2)
  expect_equal(as.numeric(m2), c(1, 2, 4))
})

test_that(".tricube weight properties", {
  d <- c(0, 1, 2)
  w <- .tricube(d)
  expect_true(all(w >= 0 & w <= 1))
  expect_equal(w[1], 1)
  expect_equal(w[3], 0)
})
