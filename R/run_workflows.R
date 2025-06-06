# run_workflows.R

utils::globalVariables(c(
  "2.5%", "25%", "50%", "75%", "97.5%",
  "Mean", "Parameter", "SD",
  "covariate", "dataset", "nsites", "psi_hat"
))


#' Run one full iSDM workflow
#'
#' This is the main workflow function to prepare, fit, and evaluate the integrated species distribution model (iSDM) using NIMBLE. It performs data validation, model input preparation, model compilation, MCMC execution, convergence diagnostics, and optional visualisation.
#'
#' @param species_datasets_list A named list of data frames, each containing species observation data including site IDs, years, and detection histories.
#' @param matching Logical. If `TRUE`, allows environmental data to be matched from the nearest 1km grid cell when unavailable at a site.
#' @param covariates_i Optional character vector. Names of spatial (site-level) covariates to include in the state model.
#' @param covariates_t Optional character vector. Names of temporal (year-level) covariates to include in the state model.
#' @param covariates_it Optional character vector. Names of spatiotemporal (site × year) covariates to include in the state model.
#' @param inclDet Logical. Whether to include detection covariates in the observation model. Defaults to `FALSE`.
#' @param inclPhenology Logical. Whether to include a seasonal phenology curve in the detection model (if date data is available).
#' @param nburnin Integer. Number of burn-in iterations for MCMC. Default is `3000`.
#' @param niter Integer. Total number of MCMC iterations. Default is `5000`.
#' @param plot_outputs Logical. Whether to plot the predicted occupancy time series and spatial map. Defaults to `FALSE`.
#' @param years Optional vector of years for plotting occupancy trends.
#' @param tag Optional character string used to tag the model estimates (e.g., for naming output models).
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{params}{A data frame of posterior parameter estimates (mean and SD), filtered to include interpretable outputs only.}
#'   \item{waic}{A data frame containing the WAIC value for model predictive accuracy.}
#' }
#'
#' @details
#' This function manages the full pipeline to run an iSDM using NIMBLE. It includes:
#' \itemize{
#'   \item Data filtering and site validation via `filter_valid_sites()`.
#'   \item Extraction and scaling of covariates via `format_state_covariates()` and `format_det_covariates()`.
#'   \item Model building with `build_imod()` and `defineModel()`.
#'   \item MCMC fitting via `run_imod()`.
#'   \item Convergence checks (Rhat, ESS) and WAIC calculation.
#'   \item Optional visualisation of predictions with `plot_occupancy_trend()` and `plot_occupancy_map()`.
#' }
#'
#' Assumes required helper functions and data processing utilities are defined and available in your package (e.g. `build_lookups()`, `build_map_layers()`, etc.).
#'
#' @examples
#' \dontrun{
#' run_model_workflow(
#'   species_datasets_list = list(dataset1 = df1, dataset2 = df2),
#'   covariates_i = c("elevationcs", "latcs"),
#'   covariates_t = c("tas_MATcs"),
#'   inclDet = TRUE,
#'   inclPhenology = TRUE,
#'   plot_outputs = TRUE,
#'   years = 2010:2020,
#'   tag = "run1"
#' )
#' }
#'
#' @importFrom dplyr filter mutate summarise bind_rows left_join select rename
#' @importFrom coda gelman.diag effectiveSize
#' @importFrom ggplot2 ggplot aes geom_line geom_point labs theme_bw theme_minimal scale_fill_viridis_c
#' @importFrom glue glue
#'
#' @export
run_model_workflow <- function(
    #-------------------------------------------------------------------------
    # List of the species observation datasets
    species_datasets_list,
    #-------------------------------------------------------------------------
    # Do we use data from the nearest grid cell if a site is missing env data?
    matching=FALSE,
    #------------------------------------------------------------------------
    # What are the names of the spatial[i], temporal[t] and spatio-temporal[it] covs?
    covariates_i = NULL,
    covariates_t = NULL,
    covariates_it = NULL,
    #-------------------------------------------------------------------------
    # Include detection covariates?
    inclDet = FALSE,
    #-------------------------------------------------------------------------
    # Should we include a phenology curve?
    inclPhenology = FALSE,
    #------------------------------------------------------------------------
    # MCMC setup
    nburnin = 3000,
    niter = 5000,
    #-----------------------------------------------------------------------
    # Should we plot the predicted occupancy trend and map?
    plot_outputs = FALSE,
    years = NULL, # If TRUE, give vector of years
    #-----------------------------------------------------------------------
    # How should we tag parameter estimate outputs (i.e., name the model)?
    tag = NULL
) {

  # Which species datasets shall we run with?-----#
  species_datasets_list <- species_datasets_list

  # Species data validation----------------------#
  ## Check we have environmental data for the sites,
  ## remove any sites without environmental data.
  ## env_filepath, env_var & env_value are pre-defined (edit if different).
  species_datasets_list <- filter_valid_sites(species_datasets_list,
                                              matching = matching)

  # Check effort over time
  print(plot_temporal_bias(species_datasets_list))

  #============================================================#
  # 2 SET DATA PARAMETERS ----
  #============================================================#

  ## imod lookups--------------------------------------
  lookups <- build_lookups(species_datasets_list)

  # Add site and year lookups to each species dataset:
  species_datasets_list <- lapply(species_datasets_list,
                                  append_lookup_tables,
                                  lookups$site_lookup,
                                  lookups$year_lookup)

  ## Metadata lookup-----------------------------------
  metadata <- get_metadata(lookups$species_samples)
  print(str(metadata)) #check

  # Check spatial overlap
  print(metadata$spatial_overlap)

  # Check detection rates summary
  print(metadata$detection_summary)
  ## naiveOcc is the fraction of sites ever occupied (≥1 detection).
  ## reportingRate is the fraction of all visits with a detection.
  ## reportingRate_z1 is the reportingRate restricted to visits at sites with historical occupancy.

  ## State covariates------------------------------------------
  ## Load in state covariates (matched to target sites/years) ----#
  all_state_covariates <- format_state_covariates(sitecodes = lookups$site_lookup,
                                                  yearcodes = lookups$year_lookup,
                                                  matching=matching,
                                                  types = c("t", "i", "it"))

  ## Print covariate names
  print(sort(names(all_state_covariates)))

  ## Select which covariates to use in the model---
  state_covariates <- build_state_covariates(all_state_covariates,
                                             covariates_i = covariates_i,
                                             covariates_t = covariates_t,
                                             covariates_it = covariates_it)

  # Remove the full suite to free up some memory---
  rm(all_state_covariates)

  ## Check model state covariates
  print(str(state_covariates))

  if (!is.null(print(str(state_covariates)))) {

    # Visualise covariate pairs plots
    plot_state_pairs(state_covariates)

  }

  ## Detection covariates----------------------------------------

  if(inclDet) {

    ### Extract all other (non-species_samples) variables from
    ### species_datasets_list
    det_covariates0 = extract_other_vars(species_datasets_list,
                                         lookups$species_samples)

    det_covariates0 <- drop_covs(det_covariates0,
                                 cols = c("survey_id", "visit_no"))

    # Format covariates
    det_covariates <- format_det_covariates(det_covariates0,
                                            impute_missing = TRUE,
                                            scale_cov = TRUE)
  } else {

    det_covariates = NULL

  }

  print(str(det_covariates)) #check

  #=============================================================#
  # 3 DEFINE MODELLING PARAMETERS----
  #=============================================================#
  imod <- build_imod(lookups$species_samples, metadata,
                     state_covariates, det_covariates,
                     inclPhenology)

  print(str(imod))

  #==========================================================#
  # 4 DEFINE MODEL CODE----
  #==========================================================#
  modelcode <- defineModel(args = imod$args)
  print(modelcode)

  #=========================================================#
  # 5 RUN MODEL----
  #=========================================================#
  # Run
  system.time(
    reptimod_out <- run_imod(code = modelcode,
                             imod$constants,
                             imod$data,
                             imod$inits,
                             imod$monitors,
                             nburnin = nburnin,
                             niter = niter)
  )


  #=======================================================#
  # 6 EVALUATION ----
  #=======================================================#

  ## Assess convergence--------------------------

  # Plot MCMC Traces
  model_traces <- plot_mcmc_diagnostics(reptimod_out$mcmc_samples)
  print(model_traces$psi_trace)
  model_traces$psi_density
  print(model_traces$monitors_density)

  # Rhat
  ### <1.1 = convergence
  rhat_vals <- gelman.diag(reptimod_out$mcmc_samples, multivariate = FALSE)$psrf[, "Point est."]

  # ESS
  ### How many “independent” draws are in the MCMC chain,
  ### Accounts for the fact that successive samples are correlated by down‐weighting for autocorrelation.
  ### ESS > 200 = posterior estimates will be precise and error will be small
  ### ESS < 100 = credible‐interval estimates have high Monte Carlo uncertainty, may need to run more iterations or improve mixing.
  ess_vals  <- effectiveSize(reptimod_out$mcmc_samples)

  diagnostics <- data.frame(Parameter = names(rhat_vals),
                            Rhat = rhat_vals,
                            ESS = ess_vals[names(rhat_vals)])

  message("Printing mcmc diagnostics...")
  print(diagnostics)

  ## Predictive power-----------------------------
  # WAIC
  ### Out-of-sample predictive accuracy that balances fit and complexity (penalises complexity)
  ### Lower WAIC = higher expected predictive accuracy
  waic <- reptimod_out$waicInfo$WAIC
  message(paste0("WAIC = ",waic))
  waic <- data.frame("WAIC" = waic)

  # 7 RESULTS -----------------------------------

  # Save estimates
  pars <- get_parameter_estimates(reptimod_out$mcmc_samples, tag = tag)

  # View estimates
  print(pars)

  covs = c(covariates_i, covariates_t, covariates_it)

  if (plot_outputs) {

    ### Plot trend----
    psi_plot <- plot_occupancy_trend(reptimod_out$mcmc_samples, years = years)
    print(psi_plot)

    ### Plot map----
    gb_covs <- build_map_layers(covariates = covs)

    map_sw <- plot_occupancy_map(gb_covs,
                                 state_covariates, pars,
                                 plot_title = NULL)
    print(map_sw)

  }

  ## Return results ---------------------------

  n_datasets = length(lookups$species_samples)
  base_params <- c("alpha.s", "Trend")
  alpha_p_params <- grep("^alpha\\.p\\d+$", pars$Parameter,
                         value = TRUE)
  psi_fs_params <- grep("^psi\\.fs", pars$Parameter, value = TRUE)
  beta_params <- paste0("b_", covs)

  if (inclPhenology && any(imod$args$has_date)) {
    phen_params = c("beta1", "phShape", "phScale")
  } else phen_params = character()

  params0 <- c(base_params, alpha_p_params, psi_fs_params,
               beta_params, phen_params)

  params <- pars %>%
    dplyr::filter(Parameter %in% params0) %>%
    dplyr::select(Parameter, Mean, SD) %>%
    dplyr::rename(parameter = "Parameter",
                  mean = "Mean",
                  sd = "SD")

  if(!is.null(tag)) {

    params$tag <- tag

  }

  return(list(params = params,
              waic = waic))

}


#' Run multiple iSDM workflow: run sensitivity analysis using custom dataset combinations
#'
#' Executes a parameter sensitivity analysis by running multiple integrated species distribution models (iSDMs), each with a different combination of input datasets. The function loops through combinations, fits models using `run_model_workflow()`, collects parameter estimates and WAIC, and saves intermediate and final results to file.
#'
#' @param spDat0 A named list of all available species observation datasets.
#' @param combinations A list of character vectors. Each element is a vector of dataset names to be used together in a model.
#' @param covariates_i Character vector of site-level covariates to include in the state model.
#' @param covariates_t Character vector of year-level covariates to include in the state model.
#' @param covariates_it Character vector of site-year covariates to include in the state model.
#' @param inclDet Logical. Whether to include detection covariates in the observation model. Defaults to `FALSE`.
#' @param inclPhenology Logical. Whether to include a phenology curve in the detection model. Defaults to `FALSE`.
#' @param nburnin Integer. Number of burn-in iterations for MCMC.
#' @param niter Integer. Total number of MCMC iterations.
#' @param plot_outputs Logical. Whether to generate and display plots during the workflow.
#' @param save_filepath String. Path prefix (including trailing slash if needed) for saving partial and final result files.
#' @param version Character string to append to output filenames (e.g., `"v1"` or `"2025-06-05"`).
#'
#' @return A list containing:
#' \describe{
#'   \item{params}{A data frame of posterior parameter estimates across all models.}
#'   \item{waic}{A data frame of WAIC values across all models.}
#' }
#'
#' @details
#' For each combination of datasets, this function:
#' \itemize{
#'   \item Subsets the datasets,
#'   \item Runs the full model via `run_model_workflow()`,
#'   \item Tags the output with a model number and dataset combination name,
#'   \item Saves intermediate `.rds` files for parameter and WAIC results,
#'   \item Saves final combined outputs once all models have run.
#' }
#'
#' Models that fail to run (e.g., due to convergence issues) are skipped and logged with a warning.
#'
#' @examples
#' \dontrun{
#' combs <- list(c("narrs", "arcarg"), c("narrs", "bto"))
#' results <- run_sensitivity(
#'   spDat0 = species_data_list,
#'   combinations = combs,
#'   covariates_i = c("elevationcs", "woodlandcs"),
#'   covariates_t = c("tas_MATcs"),
#'   covariates_it = NULL,
#'   inclDet = TRUE,
#'   inclPhenology = TRUE,
#'   nburnin = 3000,
#'   niter = 5000,
#'   plot_outputs = FALSE,
#'   save_filepath = "outputs/",
#'   version = "v1"
#' )
#' }
#'
#' @importFrom dplyr bind_rows mutate
#' @importFrom stats setNames
#' @export
run_sensitivity <- function(spDat0,
                            combinations,
                            covariates_i,
                            covariates_t,
                            covariates_it,
                            inclDet=FALSE,
                            inclPhenology=FALSE,
                            nburnin,
                            niter,
                            plot_outputs,
                            save_filepath,
                            version) {
  params_list <- list()
  waic_list   <- list()
  for(i in seq_along(combinations)) {
    set <- combinations[[i]]
    ds_list <- spDat0[set]
    model_no <- i
    model_name <- paste(set, collapse = "_")

    cat("Running model", i, "for datasets:", model_name, "…\n")

    out <- tryCatch({
      run_model_workflow(
        species_datasets_list = ds_list,
        matching = FALSE,
        covariates_i = covariates_i,
        covariates_t = covariates_t,
        covariates_it = covariates_it,
        inclDet = inclDet,
        inclPhenology = inclPhenology,
        nburnin = nburnin,
        niter = niter,
        plot_outputs = plot_outputs,
        tag = model_no
      )
    }, error = function(e) {
      warning("❌ Model ", i, " (",model_name,") failed: ", e$message)
      return(NULL)
    })

    if (is.null(out)) next

    # add tags
    params_i <- out$params %>%
      mutate(model_no = model_no,
             datasets = model_name)
    waic_i <- out$waic %>%
      mutate(model_no = model_no,
             datasets = model_name)

    cat("✅ Model", i, "(", model_name, ") completed\n\n")

    # store
    params_list[[length(params_list) + 1]] <- params_i
    waic_list[[length(waic_list) + 1]] <- waic_i

    # save partial
    partial_params <- bind_rows(params_list)
    saveRDS(partial_params,
            file = paste0(save_filepath, "sensitivity_partial_params_",version,".rds"))
    partial_waic   <- bind_rows(waic_list)
    saveRDS(partial_waic,
            file = paste0(save_filepath, "sensitivity_partial_waic_", version,".rds"))
  }

  # final combine & save
  all_params <- bind_rows(params_list)
  all_waic   <- bind_rows(waic_list)

  saveRDS(all_params,
          file = paste0(save_filepath, "sensitivity_final_params_",version,".rds"))
  saveRDS(all_waic,
          file = paste0(save_filepath, "sensitivity_final_waic_",version,".rds"))

  return(list(params = all_params,
              waic   = all_waic))
}

