#' defineModel.R
#'
#' Define a NIMBLE model for integrated species distribution model
#'
#' Constructs a flexible hierarchical Bayesian model using NIMBLE. The model includes latent
#' occupancy states, optional covariates (site, time, and site-year), detection sub-models with their
#' own covariates, and an optional phenology curve effects. The function builds and returns a `nimbleCode` object.
#'
#' @param args A named list of arguments used to customise model structure. Must include:
#' \describe{
#'   \item{`n_datasets`}{Integer. Number of separate detection sub-models (e.g., from different datasets).}
#'   \item{`covariates_i`}{Character vector. Names of 1D site-level covariates.}
#'   \item{`covariates_t`}{Character vector. Names of 1D year-level covariates.}
#'   \item{`covariates_it`}{Character vector. Names of 2D site-year covariates.}
#'   \item{`covariates_det`}{List of character vectors. One per dataset, naming detection covariates.}
#'   \item{`inclPhenology`}{Logical. Whether to include phenology terms (default = FALSE).}
#'   \item{`has_date`}{Logical vector of length `n_datasets`. Whether each detection dataset includes Julian dates.}
#' }
#'
#' @return A `nimbleCode` object defining the hierarchical state-space model, including:
#' \itemize{
#'   \item Priors for parameters governing occupancy and detection.
#'   \item A latent occupancy model indexed over sites and years.
#'   \item Detection models for each observational dataset (Bernoulli trials).
#'   \item Optionally, a phenology submodel (normalised gamma curve over Julian day).
#'   \item Derived quantities such as annual occupancy mean (`psi.fs[t]`).
#' }
#'
#' @details
#' This function dynamically generates model code as a character string using `glue()`, then evaluates
#' it into a NIMBLE-compatible model using `nimbleCode()`. Covariate names are automatically converted
#' to parameter names (e.g., `b_covname`) and included in the model formulae. Detection models are built
#' per dataset, each with its own detection covariates and phenology settings.
#'
#' @importFrom glue glue
#' @importFrom nimble nimbleCode
#'
#' @examples
#' \dontrun{
#' model_code <- defineModel(list(
#'   n_datasets = 2,
#'   covariates_i = c("elevationcs", "woodlandcs"),
#'   covariates_t = c("year_tas_anomaly"),
#'   covariates_it = c("annual_MAT"),
#'   covariates_det = list(c("mins_on_site", "list_length"), c("list_length")),
#'   inclPhenology = TRUE,
#'   has_date = c(TRUE, FALSE)
#' ))
#' }
#'
#' @export
defineModel <- function(
    args = list(n_datasets = 2, # define number of observation sub-models (n datasets)
                covariates_i = NULL,   # 1D site covs
                covariates_t = NULL,   # 1D time covs
                covariates_it = NULL, # 2D site-time covs
                covariates_det = NULL, # Observation covariates (list n_datasets)
                inclPhenology = FALSE,
                has_date = NULL)
) {# # Logical T/F vector length n_datasets
  # ---------------------------------------------------------------------
  # Build linear predictor
  lp <- "alpha.s + Trend*(t-1)"  # base terms

  # If we have site-level (1D) covariates (indexed by [i])
  if(!is.null(args$covariates_i)) {
    # Create one term per covariate: b_<name> * <name>[i]
    terms_i <- paste0("b_", args$covariates_i, " * ", args$covariates_i, "[i]")
    lp <- paste(lp, paste(terms_i, collapse = " + "), sep = " + ")
  }

  # If we have time-level (1D) covariates (indexed by [t])
  if(!is.null(args$covariates_t)) {
    terms_t <- paste0("b_", args$covariates_t, " * ", args$covariates_t, "[t]")
    lp <- paste(lp, paste(terms_t, collapse = " + "), sep = " + ")
  }

  # If we have spatio-temporal (2D) covariates (indexed by [i,t])
  if(!is.null(args$covariates_it)) {
    terms_it <- paste0("b_", args$covariates_it, " * ", args$covariates_it, "[i,t]")
    lp <- paste(lp, paste(terms_it, collapse = " + "), sep = " + ")
  }

  # ---------------------------------------------------------------------
  # Build state model priors strings for each covariate group ([i] + [t] + [i,t]).
  build_prior <- function(cov_names) {
    if(is.null(cov_names)) return("")
    paste(sapply(cov_names, function(var) {
      glue("b_{var} ~ dnorm(0, 0.1)") # All drawn from the same distribution
    }), collapse = "\n")
  }

  priors_i <- build_prior(args$covariates_i)
  priors_t <- build_prior(args$covariates_t)
  priors_it <- build_prior(args$covariates_it)

  #-----------------------------------------------------------------
  # Build observation model base priors
  # Create a string of detection priors for the n datasets as specified
  #obs_priors <- ""
  #for(j in 1:n_datasets) {
  #  obs_priors <- paste0(obs_priors,
  #                       glue("alpha.p{j} ~ dnorm(0, tau = 1)\n"))
  #}
  # Base observation model priors---
  obs_priors0 <- vector("character", length = args$n_datasets)
  for(j in seq_len(args$n_datasets)) {
    obs_priors0[j] <- glue("alpha.p{j} ~ dlogis(0, scale = 1)")
  }
  obs_priors <- paste(obs_priors0, collapse = "\n")

  #-----------------------------------------------------------------
  # Build observation model detection covariate priors---
  obs_cov_priors <- ""
  if(!is.null(args$covariates_det)){
    pri_strings <- character()
    for(cov in unlist(args$covariates_det, use.names = FALSE)){
      pri_strings <- c(pri_strings,
                       glue("b_{cov} ~ dnorm(0, 0.1)")
      )
    }
    # drop duplicates in case same cov appears more than once
    obs_cov_priors <- paste(unique(pri_strings), collapse = "\n")
  }


  #------------------------------------------------------
  # Build phenology curve block
  phen_block <- ""
  if (args$inclPhenology && any(args$has_date)) {
    phen_block <- glue('
      # build phenology priors ---
      beta1 ~ dunif(1,365)
      phShape ~ dunif(1.1,10)
      phScale ~ T(dnorm(1, tau=1),0,Inf)
      theta <- beta1/(phShape-1)

      # build  curve function ---
      for(d in 1:365) {{
        f0[d] <- dgamma(d, shape=phShape, scale=theta, log=FALSE)
        f_JD[d] <- f0[d]/max(f0[1:365])
      }}
    ')
  }


  # -------------------------------
  # Build the observation sub-models
  # For each observational dataset [j] we generate different sub-models
  # (i.e., separate likelihoods)
  obs_model_code <- ""
  for(j in seq_len(args$n_datasets)) {
    detection_terms <- ""
    covs_j <- args$covariates_det[[j]]
    if(length(covs_j)) {
      det_terms_j <- sapply(covs_j, function(cov) {
        paste0(" + b_", cov, " * ", cov, "[k]")
      })
      detection_terms <- paste(det_terms_j, collapse = "")
    }

    if(args$inclPhenology && args$has_date[j]) {
      # build with phenology
      obs_model_code <- paste0(obs_model_code, glue('
      for(k in 1:nvisit{j}) {{
        y{j}[k]  ~ dbern(Py{j}[k])
        Py{j}[k] <- z[site{j}[k], year{j}[k]] * p{j}[k]
        logit(p{j}[k]) <- alpha.p{j}{detection_terms} +
          phScale * (f_JD[JulDate{j}[k]] - max(f_JD[1:365]))
      }}
    '), "\n")
    } else {
      # build without phenology
      obs_model_code <- paste0(obs_model_code, glue('
      for(k in 1:nvisit{j}) {{
        y{j}[k]  ~ dbern(Py{j}[k])
        Py{j}[k] <- z[site{j}[k], year{j}[k]] * p{j}[k]
        logit(p{j}[k]) <- alpha.p{j}{detection_terms}
      }}
    '), "\n")
    }
  }


  # ---------------------------------------------------------------------
  # Build the full model code string
  # We use double braces {{...}} to include actual curly braces in nimbleCode.

  modelcode_str <- glue('
nimbleCode({{
  ###################### Priors -------------------------
  t.sd ~ dunif(0,5)
  t.var <- t.sd^2
  Trend ~ dnorm(0, t.var)

  mean.psi ~ dunif(0,1)
  alpha.s <- logit(mean.psi)

  {priors_i}
  {priors_t}
  {priors_it}

  ## Observation model priors
  {obs_priors}
  {obs_cov_priors}

  ## Phenology curve block
  {phen_block}

  ###################### Likelihood -------------------------

  ### State model
  for(i in 1:nsite) {{
    for(t in 1:nyear) {{
      linPred[i,t] <- {lp}
      logit(psi[i,t]) <- linPred[i,t]
      z[i,t] ~ dbern(psi[i,t])
    }}
  }}

  ## Observation sub-models
  {obs_model_code}

  ###################### Derived Quantities -------------------------
  for(t in 1:nyear) {{
    psi.fs[t] <- mean(z[1:nsite,t])
  }}

}})
')

  modelcode <- eval(parse(text = modelcode_str))
  return(modelcode)
}
