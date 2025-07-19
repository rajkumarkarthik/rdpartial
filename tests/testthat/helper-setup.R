# Global setup shared by all tests ---------------------------------------------------------

set.seed(123)
library(rdpartial)

# Reusable synthetic dataset for all tests ----------------------------------

test_data <- simulate_rdd_data(
  n            = 800,
  cutoff       = 16,
  design       = "fuzzy",
  manip_width  = 0.4,
  manip_prob   = 0.3,
  seed         = 1
)
