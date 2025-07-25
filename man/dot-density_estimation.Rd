% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/density_estimation.R
\name{.density_estimation}
\alias{.density_estimation}
\title{Estimate the number of \emph{non‑manipulated} observations to the right of the cutoff}
\usage{
.density_estimation(
  hist_df,
  manip_region,
  cutoff,
  knot_options = 3:15,
  mod_types = c("smooth", "spike_integer"),
  lambda = 100,
  eps = 1e-06,
  folds = NULL,
  num_folds = 5L,
  make_plot = FALSE,
  seed = NULL,
  solver = .get_rdpartial_opt("solver"),
  cvx_opts = .get_rdpartial_opt("cvx_opts")
)
}
\arguments{
\item{hist_df}{\code{data.frame} with columns \code{x} (numeric support points)
and \code{freq} (integer counts).}

\item{manip_region}{Numeric length‑2 vector giving the suspected manipulation
interval \emph{(lower, upper)} \strong{inclusive}.  Must contain the cutoff on its
\strong{upper} bound.}

\item{cutoff}{Numeric scalar – RDD threshold.}

\item{knot_options}{Integer vector of candidate spline knot counts;
default \code{3:15}.}

\item{mod_types}{Character vector specifying which model structures to
cross‑validate over.  Allowed values are:
\itemize{
\item "smooth" – plain spline;
\item "spike_integer" – spline + indicator for integer spikes;
\item "spike_half" – spline + integer + half‑integer spikes.
}}

\item{lambda}{Penalty weight on total predicted mass inside the manipulation
region (higher ⇒ stricter mass match).  Default \code{100} (paper default).}

\item{eps}{Numeric slack for the CVXR continuity constraints (default \code{1e-6}).}

\item{folds}{Optional integer vector assigning each \emph{donut} bin to a CV fold.
Length must equal \code{nrow(hist_df)}.  When \code{NULL}, a random \code{num_folds}‑way
split is generated.}

\item{num_folds}{Integer number of folds if \code{folds} is \code{NULL} (default \code{5}).}

\item{make_plot}{Logical; if \code{TRUE} return a \code{ggplot2} object visualising the
final fit.}

\item{seed}{Optional integer for reproducible fold assignment.}

\item{solver}{CVXR solver to use.  Defaults to
\code{.get_rdpartial_opt("solver")}.}

\item{cvx_opts}{List of solver options passed to \code{CVXR::solve()}.  Defaults to
\code{.get_rdpartial_opt("cvx_opts")}.}
}
\value{
A list with components:
\itemize{
\item \code{counts}       – \code{data.frame(x, n_true)} within \code{manip_region}.
\item \code{avg_cv_sse}   – cross‑validated SSE (per bin) for the chosen model.
\item \code{avg_sse}      – in‑sample SSE outside \code{manip_region}.
\item \code{knots}        – optimal knot count.
\item \code{model}        – optimal model type (one of \code{mod_types}).
\item \code{plot}         – \code{ggplot} object when \code{make_plot = TRUE}.
}
}
\description{
This helper replicates the constrained Poisson‑spline procedure from the
research code (Rosenman \emph{et al.}, 2025) but in a package‑friendly, testable
format.  It is \strong{not exported} – downstream functions access it internally
via \code{rdpartial:::.density_estimation()}.
}
\details{
The algorithm treats the running‑variable histogram as the realisation of a
Poisson process.  A flexible B‑spline (possibly augmented with indicator
spikes at integer or half‑integer values) is fit to the \emph{donut} region –
bins outside the user‑supplied manipulation interval.  Mass preservation and
continuity constraints are imposed through \code{CVXR}, and the optimal number of
knots / spike structure is chosen by $K$‑fold cross‑validation.  The final
step predicts the expected counts in the manipulation interval, which serve
as the "true" (non‑manipulated) totals fed into the partial‑bounds
estimators.  Bins with zero counts receive a tiny offset before logging so
that the inequality constraints remain finite.
}
\keyword{internal}
