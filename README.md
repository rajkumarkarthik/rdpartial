# rdpartial

**⚠️ Alpha Release Notice:** This package is currently in active alpha development. It is not stable and should **NOT** be used for critical work yet. Functionality may change without notice.

`rdpartial` implements the partial identification approach for regression discontinuity designs with possible manipulation of the running variable.  It provides functions for estimating bounds under both sharp and fuzzy designs, utilities for density estimation of the running variable, a simulation generator and a parametric bootstrap helper.

## Installation

The package is not on CRAN.  You can install the development version from a GitHub repository with `devtools` or `remotes`:

```r
# install.packages("devtools")
devtools::install_github("rajkumarkarthik/rdpartial")
```

## Main functions

* `simulate_rdd_data()` – generate synthetic RDD examples with optional manipulation.
* `bounds_sharp()` / `bounds_fuzzy()` – compute Manski-style bounds for sharp or fuzzy designs given estimated non-manipulated counts.
* `bootstrap_bounds()` – parametric bootstrap for confidence intervals around the bounds.

The helper `.density_estimation()` (not exported) estimates the number of non-manipulated observations in a manipulation region.

## Minimal example

```r
library(rdpartial)

# Simulate a sharp design with manipulation
set.seed(42)
sim <- simulate_rdd_data(n = 2000, cutoff = 16, design = "sharp",
                         manip_width = 0.4, manip_prob = 0.25)

# Assume 90% of post-cutoff mass is genuine
post <- sim$x[sim$x >= 16]
n_bins <- max(post) - 16 + 1
true_counts <- data.frame(
  x      = 16:max(post),
  n_true = round(tabulate(post - 16 + 1, nbins = n_bins) * 0.9)
)

# Compute lower and upper bounds at the cutoff
bounds_sharp(sim$x, sim$y, cutoff = 16, true_counts = true_counts)
```

## Acknowledgments

Development of this package was assisted by OpenAI's ChatGPT o3 model and Anthropic's Claude 4 Sonnet model for code optimization and documentation refinement. All algorithmic design, mathematical implementation, and scientific validation remain the work of the authors.

## License

This package is provided under the GPL (>= 3) license. See LICENSE for details.
