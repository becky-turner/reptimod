% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model_building_toolkit.R
\name{build_z_init}
\alias{build_z_init}
\title{Build list of initial values for z matrix}
\usage{
build_z_init(
  species_datasets_list,
  site_id = "siteID",
  year_id = "yearID",
  value_var = "y1"
)
}
\arguments{
\item{species_datasets_list}{A list of species datasets.}

\item{site_id}{Column name for site ID.}

\item{year_id}{Column name for year ID.}

\item{value_var}{Column name for detection values.}
}
\value{
A matrix representing maximum detection per site-year.
}
\description{
Constructs a site-by-year occupancy matrix from pooled species datasets.
}
\keyword{internal}
