% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting_toolkit.R
\name{predict_map}
\alias{predict_map}
\title{Predict occupancy (psi_hat) from map layers and coefficients}
\usage{
predict_map(map_layers, state_covariates, pars)
}
\arguments{
\item{map_layers}{A wide-format data frame of environmental covariates for each grid cell.}

\item{state_covariates}{List of model covariate names split by dimension.}

\item{pars}{Data frame with columns \code{Parameter} and \code{Mean}.}
}
\value{
A data frame with added columns: \code{lp} and \code{psi_hat} (logit-scale prediction and plogis).
}
\description{
Predict occupancy (psi_hat) from map layers and coefficients
}
