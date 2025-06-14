% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spatial_data_toolkit.R
\name{build_map_layers}
\alias{build_map_layers}
\title{Build map layers for selected spatial covariates}
\usage{
build_map_layers(covariates = character())
}
\arguments{
\item{covariates}{A character vector of covariate names to extract.
Polynomial terms can be specified by appending '2' (e.g., \code{"latcs2"}).}
}
\value{
A named list of data frames, each with columns \code{x1km_grid}, \code{x}, \code{y}, and \code{value},
where \code{value} is the selected covariate or its squared value.
}
\description{
Loads, prepares, and standardizes spatial environmental datasets
(e.g. temperature, land cover, elevation) for mapping. Returns a list of
data frames each containing \code{x1km_grid}, \code{x}, \code{y}, and \code{value} columns.
}
\details{
Optionally handles polynomial covariates by squaring the base variable.

Covariates are sourced from:
\itemize{
\item \code{import_climate_rasters()}
\item \code{import_spat_properties_df()}
\item \code{import_landcover_dfs()}
}
}
\examples{
\dontrun{layers <- build_map_layers(covariates = c("tas_MATcs", "latcs2"))}

}
