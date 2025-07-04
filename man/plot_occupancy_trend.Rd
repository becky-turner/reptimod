% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting_toolkit.R
\name{plot_occupancy_trend}
\alias{plot_occupancy_trend}
\title{Plot annual occupancy trend from model output}
\usage{
plot_occupancy_trend(samples, years)
}
\arguments{
\item{samples}{A coda MCMC object (e.g. output from runMCMC with samplesAsCodaMCMC = TRUE).}

\item{years}{A vector of years corresponding to the psi.fs estimates.}
}
\value{
A ggplot object visualising posterior median and credible intervals.
}
\description{
Plots predicted mean annual occupancy (psi.fs) with credible intervals.
}
