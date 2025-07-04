% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model_building_toolkit.R
\name{build_year_lookup}
\alias{build_year_lookup}
\title{Build year lookup table}
\usage{
build_year_lookup(species_datasets_list)
}
\arguments{
\item{species_datasets_list}{A named list of species datasets. Each must contain a \code{year} column.}
}
\value{
A data frame with columns \code{year} and \code{yearID}.
}
\description{
Generates a lookup table mapping each year to a unique numeric \code{yearID}.
}
\keyword{internal}
