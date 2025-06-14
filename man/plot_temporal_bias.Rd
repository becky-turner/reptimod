% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting_toolkit.R
\name{plot_temporal_bias}
\alias{plot_temporal_bias}
\title{Plot temporal bias in species datasets}
\usage{
plot_temporal_bias(species_datasets_list)
}
\arguments{
\item{species_datasets_list}{A named list of data frames, each containing columns \code{year} and \code{x1km_grid}.}
}
\value{
A ggplot object showing temporal variation in recording effort.
}
\description{
Generates a line plot showing the number of sites recorded per year for each dataset.
}
