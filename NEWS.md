# rdpartial 0.1.0

## Initial Release

* Initial implementation of partial identification bounds for regression discontinuity designs with possible manipulation of the running variable
* Implements methods described in Rosenman, Rajkumar, Gauriot & Slonim (2025) <arXiv:1910.02170>

### Main Functions

* `bounds_sharp()`: Compute Manski-style bounds for sharp RDD designs
* `bounds_fuzzy()`: Compute bounds for fuzzy RDD designs with treatment take-up
* `bootstrap_bounds()`: Parametric bootstrap for confidence intervals around bounds
* `simulate_rdd_data()`: Generate synthetic RDD data with optional manipulation

### Features

* Support for both sharp and fuzzy regression discontinuity designs
* Constrained optimization using CVXR for bound estimation
* Density estimation for detecting manipulation regions
* Bootstrap confidence intervals with parallel processing support
* Comprehensive input validation and error handling
* Extensive documentation with examples

### Dependencies

* Core optimization: CVXR
* Graphics: ggplot2
* Statistical functions: splines, stats
* Parallel processing: parallel

### License

* Released under GPL (>= 3)