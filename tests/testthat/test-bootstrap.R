test_that("bootstrap_bounds shape and determinism", {
  manip_regions <- list(c(15.7, 16))
  set.seed(42)
  boot1 <- bootstrap_bounds(
    data          = test_data,
    running_var   = "x",
    outcome       = "y",
    treatment     = "z",
    cutoff        = 16,
    manip_regions = manip_regions,
    n_boot        = 10,
    parallel      = FALSE
  )
  expect_s3_class(boot1, "rdpartial_boot")
  expect_equal(dim(boot1$boot_array), c(10, 1, 2))

  set.seed(42)
  boot2 <- bootstrap_bounds(
    data          = test_data,
    running_var   = "x",
    outcome       = "y",
    treatment     = "z",
    cutoff        = 16,
    manip_regions = manip_regions,
    n_boot        = 10,
    parallel      = FALSE
  )
  expect_equal(boot1$ci, boot2$ci)
})

test_that("bootstrap_bounds parallel matches serial", {
  skip_on_cran()
  if (.Platform$OS.type == "windows") skip("Fork not available on Windows")

  manip_regions <- list(c(15.7, 16))
  set.seed(99)
  serial <- bootstrap_bounds(
    data          = test_data,
    running_var   = "x",
    outcome       = "y",
    treatment     = "z",
    cutoff        = 16,
    manip_regions = manip_regions,
    n_boot        = 6,
    parallel      = FALSE
  )
  set.seed(99)
  parallel <- bootstrap_bounds(
    data          = test_data,
    running_var   = "x",
    outcome       = "y",
    treatment     = "z",
    cutoff        = 16,
    manip_regions = manip_regions,
    n_boot        = 6,
    parallel      = TRUE,
    n_cores       = 2
  )
  expect_equal(serial$ci, parallel$ci)
})
