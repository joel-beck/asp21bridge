% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mult_plot.R
\name{mult_plot}
\alias{mult_plot}
\title{Plot multiple Markov Chains together}
\usage{
mult_plot(
  samples,
  type = c("trace", "density", "both"),
  log = FALSE,
  free_scale = FALSE,
  latex = FALSE
)
}
\arguments{
\item{samples}{Model object of the class 'lmls', list or matrix
containing samples of the MCMC sampler.}

\item{type}{One of the values "trace", "density" or "both". \cr
Controls whether only trace plots or only density plots should
be included for all chains.
The option "both" stacks trace plots on top of the density plots.}

\item{log}{Logical. If TRUE, the y - axis in trace plots and the x - axis in
density plots are transformed to the logarithmic scale. \cr
Default: FALSE}

\item{free_scale}{Logical. If TRUE, the y - axis scale in trace plots is
chosen differently for all facets. \cr
Default: FALSE}

\item{latex}{Logical. If TRUE, facets in trace plots and ridgelines in density
plots can be labeled by rendered mathematical symbols such as
greek letters in combination with subscripts and superscripts.
This option requires the columns of the input matrices to be
labeled accordingly (as in the output of the \code{mcmc_ridge()}
function). \cr
Default: FALSE}
}
\value{
Plot object of the class "ggplot".
}
\description{
The \code{mult_plot()} function allows to combine trace and / or
density plots of multiple Markov Chains into a single graphical
illustration. \cr
Trace plots are displayed separately in facets, density plots
are shown in a ridgeline plot.
}
\examples{
fit <- lmls(
  location = y ~ x1 + x2 + z1 + z2, scale = ~ z1 + z2,
  data = toy_data, light = FALSE
) \%>\%
  mcmc_ridge(num_sim = 1000)

# list of 4 matrices with 1000 rows each
samples <- fit$mcmc_ridge$sampling_matrices

# plots all output list elements together for quick overview
mult_plot(samples, type = "both", free_scale = TRUE, latex = TRUE)

# trace plot of location coefficients
mult_plot(samples$location, type = "trace", free_scale = TRUE, latex = TRUE)

# log argument often useful for strictly positive variance parameters
mult_plot(
  list(samples$location_prior, samples$scale_prior),
  type = "both", log = TRUE, latex = TRUE
)

}
