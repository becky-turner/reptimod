
#' Extract posterior parameter estimates
#'
#' Converts MCMC samples into a tidy data frame of parameter estimates, including mean, SD, and quantiles.
#'
#' @param samples An `mcmc.list` object, typically the output of `runMCMC(..., samplesAsCodaMCMC = TRUE)`.
#' @param tag Optional character label to tag the resulting data frame (e.g., a model identifier).
#'
#' @return A data frame with columns: `Parameter`, `Mean`, `SD`, `Naive SE`, `Time-series SE`, and optionally `model_tag`.
#'
#' @examples
#' \dontrun{
#'   samples <- runMCMC(...)
#'   get_parameter_estimates(samples, tag = "Model1")
#' }
#'
#' @export
get_parameter_estimates <- function(samples, tag=NULL) {

  pars <- as.data.frame(summary(samples)$statistics)
  pars$Parameter <- rownames(pars)
  row.names(pars) <- NULL
  pars <- pars[,c(5,1:4)]

  if(!is.null(tag)){

    pars$model_tag <- tag

  }

  return(pars)

}

#' Plot MCMC diagnostics
#'
#' Generates diagnostic plots for a fitted NIMBLE model's MCMC output, including traceplots and density plots for state (`psi.fs`) and other parameters.
#'
#' @param samples An `mcmc.list` object (e.g., from `runMCMC(..., samplesAsCodaMCMC = TRUE)`).
#' @param greek Logical. If `TRUE`, Greek symbols will be used in plot labels (default: `TRUE`).
#'
#' @return A named list of `ggplot2` objects:
#'   \item{psi_trace}{Traceplots for the `psi.fs` parameter.}
#'   \item{psi_density}{Posterior density plots for `psi.fs`.}
#'   \item{monitors_trace}{Traceplots for all other monitored parameters.}
#'   \item{monitors_density}{Posterior density plots for all other monitored parameters.}
#'
#' @examples
#' \dontrun{
#'   diagnostics <- plot_mcmc_diagnostics(samples)
#'   diagnostics$psi_trace
#' }
#'
#' @importFrom ggmcmc ggs ggs_traceplot ggs_density
#' @importFrom ggplot2 ggtitle facet_wrap
#' @importFrom coda as.mcmc.list
#'
#' @export
plot_mcmc_diagnostics <- function(samples, greek = TRUE) {
  # samples: an mcmc.list (the output of runMCMC(..., samplesAsCodaMCMC=TRUE))

  psi = "psi.fs"

  # 1) traceplots for the psi family
  psi_trace <- samples %>%
    ggs(family = psi) %>%
    ggs_traceplot(greek = greek) +
    facet_wrap(~ Parameter) +
    ggtitle(paste0("Traceplots for '", psi, "'"))

  # 2) density for the psi family
  psi_density <- samples %>%
    ggs(family = psi) %>%
    ggs_density(greek = greek) +
    facet_wrap(~ Parameter) +
    ggtitle(paste0("Density plots for '", psi, "'"))

  # 3) extract everything *except* psi
  param_names <- dimnames(samples[[1]])[[2]]
  keep <- !grepl(psi, param_names)
  samples2 <- as.mcmc.list(
    lapply(samples, function(chain) chain[, keep, drop = FALSE])
  )

  # 4) traceplots for the other parameters
  monitors_trace <- samples2 %>%
    ggs() %>%
    ggs_traceplot(greek = greek) +
    facet_wrap(~ Parameter, scales = "free") +
    ggtitle("Traceplots for monitored parameters")

  # 5) posterior densities for the other parameters
  monitors_density <- samples2 %>%
    ggs() %>%
    ggs_density(greek = greek) +
    facet_wrap(~ Parameter, scales = "free") +
    ggtitle("Posterior densities for montired parameters")

  return(list(
    psi_trace = psi_trace,
    psi_density = psi_density,
    monitors_trace = monitors_trace,
    monitors_density = monitors_density
  ))
}
