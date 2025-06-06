#' run_imod.R
#'
#' Run integrated occupancy-detection model using NIMBLE
#'
#' Compiles and runs a Bayesian model using NIMBLE, based on a supplied model code object, constants,
#' data, initial values, and monitored parameters. Returns posterior samples and WAIC for model evaluation.
#'
#' @param code A NIMBLE model definition object (as returned by `nimbleCode()`).
#' @param constants A named list of model constants (e.g., number of sites, years).
#' @param data A named list of observed data for the model (e.g., detection/non-detection matrix).
#' @param inits A named list of initial values for latent states and parameters.
#' @param monitors A character vector naming the parameters to monitor in MCMC sampling.
#' @param nburnin Number of iterations for burn-in (default = 3000).
#' @param niter Total number of MCMC iterations (default = 5000).
#'
#' @return A list with two elements:
#' \describe{
#'   \item{`mcmc_samples`}{An object of class `mcmc.list` containing posterior samples for the monitored parameters.}
#'   \item{`waicInfo`}{WAIC statistics for model comparison and assessment.}
#' }
#'
#' @details
#' This function performs the following steps:
#' \enumerate{
#'   \item Builds the NIMBLE model from the specified code, constants, data, and inits.
#'   \item Constructs an MCMC object with the requested monitors and no conjugacy.
#'   \item Compiles both the model and the MCMC.
#'   \item Runs the MCMC sampler with 3 chains, specified burn-in, and total iterations.
#'   \item Returns posterior samples as a `coda::mcmc.list`, with WAIC if enabled.
#' }
#'
#' Any model build or run failure will stop with an informative error message.
#'
#' @importFrom nimble nimbleModel buildMCMC compileNimble runMCMC
#' @importFrom coda mcmc.list
#'
#' @examples
#' \dontrun{
#' results <- run_imod(code = my_model_code,
#'                     constants = model_constants,
#'                     data = model_data,
#'                     inits = model_inits,
#'                     monitors = c("beta0", "psi", "p"),
#'                     nburnin = 2000,
#'                     niter = 4000)
#' }
#'
#' @export
run_imod <- function(code,
                     constants, # List of model constants
                     data, # List of y1 obs data
                     inits, # List of initial values
                     monitors, # Names of parameters to monitor in MCMC
                     nburnin = 3000 , # nburnin number value for MCMC
                     niter = 5000) { # niter number value for MCMC


  #─── 0) CHECK INPUTS ───────────────────────────────────────
  cat("▶ Length constants:", length(constants),  "\n")
  cat("▶ Names(constants):", paste(names(constants), collapse = ", "), "\n\n")
  cat("▶ Length data:", length(data),      "\n")
  cat("▶ Names(data):", paste(names(data), collapse = ", "), "\n\n")
  cat("▶ Length inits:", length(inits),    "\n")
  cat("▶ Names(inits):", paste(names(inits), collapse = ", "), "\n\n")

  #------------------------------------------------------------
  # Step 1: Create operational model

  model <- tryCatch({
    nimbleModel(code = code,
                constants = constants,
                data = data,
                inits = inits)
  }, error = function(e) {
    stop("nimbleModel() failed:\n", e$message)
  })
  cat("✅ nimbleModel() succeeded\n\n")
  model$initializeInfo()


  #----------------------------------------------------------
  # Step 2: Build an MCMC object using buildMCMC().

  occMCMC <- tryCatch({
    buildMCMC(model,
              monitors = monitors,
              thin = 3,
              useConjugacy = FALSE,
              enableWAIC = TRUE)
  }, error = function(e) {
    stop("buildMCMC() failed:\n", e$message)
  })
  #about 25 seconds
  cat("✅ nimbleModel() constructed.\n\n")


  #----------------------------------------------------------
  # Step 3: Before compiling the MCMC object, we need to
  # compile the model first
  Cmodel <- tryCatch({
    compileNimble(model)
  }, error = function(e) {
    stop("Cmodel failed:\n", e$message)
  })
  # 25 seconds (less for fewer nodes)

  #--------------------------------------------------------
  # Step 4: Compile the MCMC (project = NIMBLE model already associated with a project)
  CoccMCMC <- tryCatch({
    compileNimble(occMCMC, project = Cmodel)
  }, error = function(e) {
    stop("compileNimble MCMC failed:\n", e$message)
  })
  cat("✅ Build and compile MCMC OK\n\n")
  # 5 mins?

  cat("✅ MCMC compiled successfully.\n\n")

  #-------------------------------------------------------
  # Step 5: Run the MCMC. either $run or runMCMC() on the compiled model object.
  runMCMC_samples <- tryCatch({
    runMCMC(CoccMCMC,
            nburnin = nburnin,
            niter = niter,
            nchains = 3,
            samplesAsCodaMCMC = T,
            WAIC = TRUE)
  }, error = function(e) {
    cat("❌ MCMC run failed with error:\n", e$message, "\n\n")
    stop(e)
  })

  # Remove large objects
  rm(model, occMCMC, Cmodel, CoccMCMC)
  gc()
  return(list(mcmc_samples = runMCMC_samples$samples,
              waicInfo = runMCMC_samples$WAIC))

}
