# rdpartial – global definitions -------------------------------------------
#
# This file registers non‑standard evaluation (NSE) symbols so that `R CMD
# check` remains NOTE‑free, and sets package‑wide **options** that other
# functions can retrieve via `.get_rdpartial_opt()`.
#
# Nothing in this file is exported to the user; everything is marked
# `@keywords internal`.
# ---------------------------------------------------------------------------

#' @keywords internal
utils::globalVariables(
  c(
    # ggplot2 / dplyr column‑name tokens ---------------------------------
    "hlevel", "freq", "Freq", "manip", "pred", "avg.prop",
    "leftLine", "rightLine", "right.upr", "right.lwr",
    "numTrueSubjects", "h.level", "lr", "size", "fill"
  )
)

# ---- Default package‑level constants ---------------------------------------
# CVXR solver and numerical tolerances used across optimisation routines.
.DEFAULT_SOLVER   <- "ECOS"
.DEFAULT_CVX_OPTS <- list(reltol = 1e-6, feastol = 1e-6, abstol = 1e-6)

# ---- Option setter on package load -----------------------------------------
.onLoad <- function(libname, pkgname) {
  op <- options()
  op_rd <- list(
    rdpartial.solver    = .DEFAULT_SOLVER,
    rdpartial.cvx_opts  = .DEFAULT_CVX_OPTS
  )
  # only set those that are not already defined -----------------------------
  to_set <- setdiff(names(op_rd), names(op))
  if (length(to_set)) options(op_rd[to_set])
  invisible()
}

# ---- Internal accessor -----------------------------------------------------
#' Retrieve an **rdpartial** option (internal helper)
#'
#' Thin wrapper around `getOption()` that prefixes the requested name with
#' "rdpartial." and provides a helpful error if the option has never been set.
#'
#' @param name Character scalar, e.g. "solver" or "cvx_opts".
#' @return The option value.
#' @keywords internal
.get_rdpartial_opt <- function(name) {
  .assert_scalar_numeric(1, "dummy") # trick to import assert helper for R CMD check
  val <- getOption(paste0("rdpartial.", name))
  if (is.null(val)) {
    stop("Unknown rdpartial option: ", name, call. = FALSE)
  }
  val
}
