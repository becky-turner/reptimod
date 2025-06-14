% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model_building_toolkit.R
\name{append_lookup_tables}
\alias{append_lookup_tables}
\title{Append site and year lookups to species datasets}
\usage{
append_lookup_tables(spDat, site_lookup, year_lookup)
}
\arguments{
\item{spDat}{A data frame representing a species dataset.}

\item{site_lookup}{A data frame of site IDs from \code{build_site_lookup()}.}

\item{year_lookup}{A data frame of year IDs from \code{build_year_lookup()}.}
}
\value{
A data frame with additional columns \code{siteID} and \code{yearID}, ordered by year, visit and site.
}
\description{
Adds \code{siteID} and \code{yearID} columns to a single species dataset by joining with the lookup tables.
}
\keyword{internal}
