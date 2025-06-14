% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model_building_toolkit.R
\name{extract_other_vars}
\alias{extract_other_vars}
\title{Extract detection covariates}
\usage{
extract_other_vars(extract_from, comparator = species_samples)
}
\arguments{
\item{extract_from}{List of full data frames.}

\item{comparator}{Reference list of species samples.}
}
\value{
A list of data frames containing only the additional covariates.
}
\description{
Removes columns present in the core species samples from each dataset, leaving detection covariates.
}
\keyword{internal}
