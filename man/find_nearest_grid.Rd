% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spatial_data_toolkit.R
\name{find_nearest_grid}
\alias{find_nearest_grid}
\title{Find nearest grid cell by euclidean distance}
\usage{
find_nearest_grid(data1, data2)
}
\arguments{
\item{data1}{A data frame of unmatched sites with column \code{x1km_grid}.}

\item{data2}{A data frame of reference sites with column \code{x1km_grid}.}
}
\value{
A data frame of the nearest matching grid rows from \code{data2}.
}
\description{
Matches each site in \code{data1} to the nearest grid square in \code{data2}, using OSGB eastings/northings.
}
