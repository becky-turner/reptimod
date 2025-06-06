# model_building.R
#
# This script contains functions for constructing the full suite of model inputs required to
# run an integrated occupancy-detection model in NIMBLE. These include generation of:
# - site and year lookup tables
# - formatted input datasets
# - covariate constants
# - initial values
# - arguments and monitors for model specification
#
# These functions assume that the input data have been pre-processed to include standard
# fields (e.g., 'x1km_grid', 'year', 'y1') and, optionally, a 'date' column if phenology is used.
#
# Exported functions:
# - build_lookups(): Prepares site/year lookup tables and species_samples list
# - get_metadata(): Extracts peak-day, overlap, and detection summaries
# - build_model_constants(): Assembles constants list for nimble model
# - build_model_data(): Extracts and formats detection data for model input
# - build_model_inits(): Compiles all initial values for MCMC
# - build_model_arguments(): Constructs arguments list for use in defineModel()
# - build_monitors(): Specifies parameters to monitor during MCMC
# - build_imod(): Wrapper to return all required components for running NIMBLE iSDM


#' @importFrom dplyr %>% left_join arrange select group_by summarise filter pull all_of any_of bind_rows n_distinct
#' @importFrom magrittr %>%
#' @importFrom tibble tibble
#' @importFrom reshape2 acast
#' @importFrom stats as.formula rnorm runif sd setNames
#' @importFrom utils globalVariables

utils::globalVariables(c(
  "x1km_grid", "year", "y1", "siteID", "yearID", "visit_no",
  "ever_detected", "total_count", "species_samples", "metadata"
))

#' Build site lookup table
#'
#' Generates a unique lookup table mapping each 1km grid cell (site) to a unique numeric `siteID`.
#'
#' @param species_datasets_list A named list of species datasets. Each must contain a `x1km_grid` column.
#'
#' @return A data frame with columns `x1km_grid` and `siteID`.
#' @keywords internal
build_site_lookup <- function(species_datasets_list) {
  all_sites <- unique(unlist(lapply(species_datasets_list, function(df) df$x1km_grid)))
  site_lookup <- data.frame(x1km_grid = all_sites, siteID = seq_along(all_sites), stringsAsFactors = FALSE)
  return(site_lookup)
}

#' Build year lookup table
#'
#' Generates a lookup table mapping each year to a unique numeric `yearID`.
#'
#' @param species_datasets_list A named list of species datasets. Each must contain a `year` column.
#'
#' @return A data frame with columns `year` and `yearID`.
#' @keywords internal
build_year_lookup <- function(species_datasets_list) {
  all_years <- unique(unlist(lapply(species_datasets_list, function(df) df$year)))
  all_years <- sort(all_years)
  year_lookup <- data.frame(year = all_years, yearID = seq_along(all_years), stringsAsFactors = FALSE)
  return(year_lookup)
}

#' Append site and year lookups to species datasets
#'
#' Adds `siteID` and `yearID` columns to a single species dataset by joining with the lookup tables.
#'
#' @param spDat A data frame representing a species dataset.
#' @param site_lookup A data frame of site IDs from `build_site_lookup()`.
#' @param year_lookup A data frame of year IDs from `build_year_lookup()`.
#'
#' @importFrom dplyr left_join arrange
#' @importFrom magrittr %>%
#'
#' @return A data frame with additional columns `siteID` and `yearID`, ordered by year, visit and site.
#' @keywords internal
append_lookup_tables <- function(spDat, site_lookup, year_lookup) {
  spDat <- spDat %>%
    dplyr::left_join(site_lookup, by = "x1km_grid") %>%
    dplyr::left_join(year_lookup, by = "year")

  if ("visit_no" %in% colnames(spDat)) {
    spDat <- spDat %>% dplyr::arrange(yearID, visit_no, siteID)
  } else {
    spDat <- spDat %>% dplyr::arrange(yearID, siteID)
  }
  return(spDat)
}

#' Build list of species samples data
#'
#' Extracts relevant fields for model input from each species dataset, returning a clean list.
#'
#' @param species_datasets_list A list of species datasets with `x1km_grid`, `year`, `y1`, `siteID`, `yearID`, and optionally `date`.
#' @param date_id The name of the date column (default = "date").
#'
#' @importFrom dplyr select all_of
#'
#' @return A list of cleaned species samples suitable for input to model functions.
#' @keywords internal
build_species_samples <- function(species_datasets_list, date_id = "date") {
  lapply(species_datasets_list, function(df) {
    if (date_id %in% colnames(df)) {
      df <- df %>% select(x1km_grid, year, all_of(date_id), y1, siteID, yearID)
    } else {
      df <- df %>% select(x1km_grid, year, y1, siteID, yearID)
    }
  })
}

#' Extract detection covariates
#'
#' Removes columns present in the core species samples from each dataset, leaving detection covariates.
#'
#' @param extract_from List of full data frames.
#' @param comparator Reference list of species samples.
#'
#' @importFrom dplyr select any_of
#'
#' @return A list of data frames containing only the additional covariates.
#' @keywords internal
extract_other_vars <- function(extract_from, comparator = species_samples) {
  excl_cols = unique(unlist(lapply(comparator, names)))
  lapply(extract_from, function(df) df %>% select(-any_of(excl_cols)))
}

#' Build NIMBLE model lookups
#'
#' Prepares and returns a set of lookup tables and cleaned species samples for model input.
#'
#' @param species_datasets_list A named list of species datasets.
#'
#' @return A list containing `site_lookup`, `year_lookup`, and `species_samples`.
#' @export
build_lookups <- function(species_datasets_list) {
  site_lookup <- build_site_lookup(species_datasets_list)
  year_lookup <- build_year_lookup(species_datasets_list)
  datasets_with_lookups <- lapply(species_datasets_list, function(df) {
    append_lookup_tables(df, site_lookup, year_lookup)
  })
  species_samples <- build_species_samples(datasets_with_lookups)
  return(list(site_lookup = site_lookup,
              year_lookup = year_lookup,
              species_samples = species_samples))
}

#' Check spatial overlap
#'
#' Computes summary statistics on the spatial overlap between the species datasets.
#'
#' Returns total number of unique sites, number of shared sites across all datasets,
#' and a pairwise matrix of sites overlap.
#'
#' @param species_datasets_list A named list of species datasets.
#'
#' @return A list with elements: `total_unique`, `sites_in_all`, and `pairwise` matrix.
#' @keywords internal
check_spatial_overlap <- function(species_datasets_list) {
  # Extract the unique x1km_grid values from each species dataset.
  site_list <- lapply(species_datasets_list, function(df) unique(df$x1km_grid))
  total_unique <- length(unique(unlist(site_list)))
  sites_in_all <- length(Reduce(intersect, site_list))
  n <- length(site_list)
  pairwise <- matrix(0, nrow = n, ncol = n)
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      pairwise[i, j] <- length(intersect(site_list[[i]], site_list[[j]]))
    }
  }
  if (!is.null(names(species_datasets_list))) {
    rownames(pairwise) <- names(species_datasets_list)
    colnames(pairwise) <- names(species_datasets_list)
  }
  return(list(total_unique = total_unique,
              sites_in_all = sites_in_all,
              pairwise = pairwise))
}

#' Extract datasets metadata
#'
#' Computes metadata summaries including peak detection date, spatial overlap, and detection metrics.
#'
#' @param species_datasets_list A named list of species datasets.
#'
#' @importFrom dplyr group_by summarise slice pull
#' @importFrom stats sd
#' @importFrom tibble tibble
#'
#' @return A list with keys: `totalSites`, `totalYears`, `peak_per_year`, `peakDay_sd`, `peakDay`,
#'   `datasets_for_peakStats`, `spatial_overlap`, `detection_summary`.
#' @export
get_metadata <- function(species_datasets_list) {
  # Identify which datasets have a date column
  has_date <- sapply(species_datasets_list, function(df) "date" %in% colnames(df))
  used_datasets <- names(species_datasets_list)[has_date]
  if (is.null(used_datasets) && any(has_date)) {
    used_datasets <- which(has_date)
  }
  # Compute total sites & years across *all* datasets
  combined <- dplyr::bind_rows(species_datasets_list)
  #-----------------------------------------------------#
  # Total sites
  totalSites <- dplyr::n_distinct(combined$siteID)
  #-----------------------------------------------------#
  # Total years
  totalYears <- dplyr::n_distinct(combined$yearID)

  #------------------------------------------------#
  # Calculate peak day per year
  # If none of the datasets have dates, set peaks to NA and note it in a message
  if (length(used_datasets) == 0) {
    peak_per_year <- tibble()
    peakDay_sd <- NA_real_
    peakDay <- NA_real_
    msg <- "No datasets contained a 'date' column; peak‐day metrics unavailable."
  } else {
    # Combine only datasets which have date column
    dated_list  <- species_datasets_list[has_date]
    all_dat_dated <- dplyr::bind_rows(dated_list) %>% dplyr::select(year, date, y1)
    peak_per_year <- all_dat_dated %>%
      group_by(year) %>%
      summarise(peak_day = date[which.max(y1)], .groups = "drop")
    #------------------------------------------------#
    # Calculate sd of peak days
    peakDay_sd <- sd(peak_per_year$peak_day)
    #------------------------------------------------#
    # Calculate overall peak day
    peakDay <- all_dat_dated %>%
      group_by(date) %>%
      summarise(total_count = sum(y1, na.rm = TRUE), .groups = "drop") %>%
      slice(which.max(total_count)) %>%
      pull(date)
    #------------------------------------------------#
    # Save info on what datasets were used
    msg <- paste0("Peak‐day metrics calculated using dataset",
                  if (length(used_datasets)>1) "s " else " ",
                  paste(used_datasets, collapse = ", "))
  }

  #------------------------------------------------#
  # Check spatial overlap
  spatial_overlap = check_spatial_overlap(species_datasets_list)
  #------------------------------------------------#
  # Summarise datasets
  detection_summary = summarise_data(species_datasets_list)
  print(msg)
  return(list(totalSites = totalSites,
              totalYears = totalYears,
              peak_per_year = peak_per_year,
              peakDay_sd = peakDay_sd,
              peakDay = peakDay,
              datasets_for_peakStats = used_datasets,
              spatial_overlap = spatial_overlap,
              detection_summary = detection_summary))
}

#' Build model constants
#'
#' Constructs a named list of constants required for the nimble iSDM.
#'
#' Includes sample dimensions, visit counts, site/year IDs, and optional phenology and covariate constants.
#'
#' @param species_datasets_list A list of species datasets with site and year information.
#' @param state_covariates_list Optional list of state covariates.
#' @param ob_covariates_list Optional list of observation covariates.
#' @param inclPhenology Logical; should phenology constants be included?
#'
#' @return A named list of constants.
#' @export
build_model_constants <- function(species_datasets_list,
                                  state_covariates_list = NULL,
                                  ob_covariates_list = NULL,
                                  inclPhenology) {
  # 1) Full sample dimensions
  all_sites <- unlist(lapply(species_datasets_list, `[[`, "siteID"))
  all_years <- unlist(lapply(species_datasets_list, `[[`, "yearID"))
  sample_constants <- list(
    nsite = length(unique(all_sites)),
    nyear = length(unique(all_years))
  )

  # 2) per‐dataset visit counts, dates, and site/year indices
  for (j in seq_along(species_datasets_list)) {

    df <- species_datasets_list[[j]]

    if(inclPhenology == TRUE) {

      if("date" %in% colnames(df)){
        sample_constants[[paste0("JulDate", j)]] <- as.integer(df$date)
      } else {
        warning(
          "inclPhenology = TRUE but none of the species datasets has a 'date' column.\n",
          "  Ignoring phenology and will not add JulDate to constants."
        )
      }
    }

    sample_constants[[paste0("nvisit", j)]] <- nrow(df)
    sample_constants[[paste0("site",j)]] <- as.integer(df$siteID)
    sample_constants[[paste0("year",j)]] <- as.integer(df$yearID)
  }

  # 3) state‐covariate constants
  state_constants <- if(is.null(state_covariates_list)) list() else {
    sc <- unlist(state_covariates_list, recursive = FALSE, use.names = TRUE)
    names(sc) <- sub("^[^.]*\\.", "", names(sc))
    sc
  }

  # 4) obs‐cov constants
  obs_constants <- if(is.null(ob_covariates_list)) list() else {
    oc <- unlist(ob_covariates_list, recursive = FALSE, use.names = TRUE)
    names(oc) <- sub("^[^.]*\\.", "", names(oc))
    oc
  }

  c(sample_constants, state_constants, obs_constants)
}

#' Build model data
#'
#' Constructs a named list of response variables (y1) for NIMBLE input.
#'
#' @param species_datasets_list A list of species dataframes with y1 values.
#'
#' @importFrom stats setNames
#'
#' @return A named list where each element is a numeric vector of y1 values.
#' @export
build_model_data <- function(species_datasets_list) {
  setNames(lapply(species_datasets_list, function(df) as.numeric(df$y1)),
           paste0("y", seq_along(species_datasets_list)))
}

#' Build list of initial values for covariates
#'
#' Generates initial values for covariate effects in the occupancy-detection model.
#'
#' @param state_covariates_list Optional list of state covariates.
#' @param det_covariates_list Optional list of detection covariates.
#'
#' @importFrom stats rnorm
#'
#' @return A named list of initial values for covariate coefficients.
#' @keywords internal
build_b_inits <- function(state_covariates_list = NULL,
                          det_covariates_list = NULL) {

  # helper to flatten and strip weird prefixes
  flatten_names <- function(lst) {
    if(is.null(lst)) return(character())
    # lst is a named list of vectors,
    # so names(lst) are the covariate names already:
    #names(lst)
    covs <- unlist(lapply(lst, names), use.names = FALSE)
    return(covs)
  }

  # grab the names out of each list (or empty if NULL)
  state_inits  <- flatten_names(state_covariates_list)
  obs_inits  <- flatten_names(det_covariates_list)
  all_cov_inits <- c(state_inits, obs_inits)

  # if there are no covariates, return an empty list
  if (length(all_cov_inits) == 0) {
    return(list())
  }

  # for each covariate name draw a single rnorm(1,0,1)
  b_inits <- lapply(all_cov_inits, function(x) rnorm(1, 0, 1))

  # name them b_<covariate>
  names(b_inits) <- paste0("b_", all_cov_inits)

  return(b_inits)
}

#' Build list of initial values for detection intercept(s)
#'
#' Creates a named list of initial values for detection intercepts per dataset.
#'
#' @param n_datasets Number of datasets.
#' @param init_value Initial value to assign to each alpha.p parameter.
#'
#' @return A named list of initial values.
#' @keywords internal
build_alpha_p_inits <- function(n_datasets, init_value = 0.5) {
  vals <- rep(init_value, n_datasets)
  names(vals) <- paste0("alpha.p", seq_len(n_datasets))
  as.list(vals)
}

#' Build list of initial values for z matrix
#'
#' Constructs a site-by-year occupancy matrix from pooled species datasets.
#'
#' @param species_datasets_list A list of species datasets.
#' @param site_id Column name for site ID.
#' @param year_id Column name for year ID.
#' @param value_var Column name for detection values.
#'
#' @importFrom reshape2 acast
#' @importFrom stats as.formula
#'
#' @return A matrix representing maximum detection per site-year.
#' @keywords internal
build_z_init <- function(species_datasets_list,
                         site_id = "siteID",
                         year_id = "yearID",
                         value_var = "y1") {
  all_species_samples <- do.call(rbind, lapply(species_datasets_list, function(df) {
    df[, c(site_id, year_id, value_var)]
  }))
  # create matrix of max(y) per site.x.year, filling missing with 0
  z.init <- acast(all_species_samples,
                  formula = as.formula(paste0(site_id, " ~ ", year_id)),
                  value.var = value_var, fun = max, fill = 0)
  return(z.init)
}

#' Build list of initial values for model
#'
#' Compiles initial values for latent state, core parameters, phenology, and covariate coefficients.
#'
#' @param species_datasets_list List of species datasets.
#' @param state_covariates Optional list of state covariates.
#' @param det_covariates Optional list of detection covariates.
#' @param peakDay Peak day used for phenology initialisation.
#' @param n_datasets Number of species datasets.
#' @param inclPhenology Logical indicating whether to include phenology parameters.
#'
#' @importFrom stats rnorm runif
#'
#' @return A list of named initial values.
#' @export
build_model_inits <- function(species_datasets_list,
                              state_covariates=NULL,
                              det_covariates=NULL,
                              peakDay = metadata$peakDay,
                              n_datasets = length(species_datasets_list),
                              inclPhenology = FALSE) {
  # z matrix
  z.init <- build_z_init(species_datasets_list)

  # start building fixed inits
  fixed <- list(z = z.init,
                t.sd = 1,
                mean.psi = runif(1,0,1),
                alpha.s = 0,
                Trend = rnorm(1, 0, 0.2))

  # add alpha.ps
  fixed <- c(fixed, build_alpha_p_inits(n_datasets, init_value = 0.2))

  # optionally add phenology inits IF we have a date column
  has_date = vapply(species_datasets_list,
                    function(df) "date" %in% names(df),
                    logical(1))

  if (inclPhenology && any(has_date)) {
    fixed <- c(fixed, list(beta1 = peakDay, phShape = 2, phScale = 1))
  }

  if (inclPhenology && !any(has_date)) {
    warning(
      "inclPhenology = TRUE but none of the species datasets has a 'date' column.\n",
      "  Ignoring phenology and building inits without 'beta1', 'phShape', 'phScale'."
    )
  }

  # covariate slopes
  b_inits <- build_b_inits(state_covariates, det_covariates)

  # combine and return
  c(fixed, b_inits)
}

#' Build arguments for defineModel()
#'
#' Prepares a list of arguments for defining the iSDM modelcode.
#'
#' @param species_samples List of formatted species datasets.
#' @param state_covariates Named list of state covariates.
#' @param det_covariates Named list of detection covariates.
#' @param inclPhenology Logical indicating whether to include phenology terms.
#'
#' @return A list of arguments used in defineModel.
#' @export
build_model_arguments <- function(species_samples,
                                  state_covariates,
                                  det_covariates,
                                  inclPhenology=FALSE) {
  # number of species datasets
  n_datasets = length(species_samples)

  # state model covariates (or NULL if none)
  covariates_i  = {
    x <- names(state_covariates$covariates_i)
    if (length(x)) x else NULL
  }
  covariates_t  = {
    x <- names(state_covariates$covariates_t)
    if (length(x)) x else NULL
  }
  covariates_it = {
    x <- names(state_covariates$covariates_it)
    if (length(x)) x else NULL
  }

  # detection covariates (or NULL if none)
  covariates_det = {
    x <- get_cov_names(det_covariates, format = "list")
    if (length(x)) x else NULL
  }

  # phenology
  inclPhenology = inclPhenology

  # build first part of arguments
  out <- list(n_datasets = n_datasets,
              covariates_i = covariates_i,
              covariates_t = covariates_t,
              covariates_it = covariates_it,
              covariates_det = covariates_det,
              inclPhenology = inclPhenology)

  if(inclPhenology == TRUE) {

    has_date_check <- vapply(species_samples,
                             function(df) "date" %in% names(df), logical(1))
    out$inclPhenology <- TRUE
    # which datasets have a date?
    out$has_date = sapply(species_samples, function(spDat) "date" %in% names(spDat))

    if(!any(has_date_check)) {
      warning(
        "inclPhenology = TRUE but one or more species datasets does not have a 'date' column.")
    }
  }

  # output without phenology curve arguments
  return(out)
}

#' Build list of parameters to monitor for NIMBLE
#'
#' Generates a list of model parameters to monitor during MCMC.
#'
#' @param state_covariates_list Optional list of state covariates.
#' @param det_covariates_list Optional list of detection covariates.
#' @param n_datasets Number of datasets.
#' @param inclPhenology Logical indicating whether to include phenology parameters.
#'
#' @return A character vector of parameter names to monitor.
#' @export
build_monitors <- function(state_covariates_list = NULL,
                           det_covariates_list = NULL,
                           n_datasets,
                           inclPhenology=FALSE) {
  #---------------------------------------------------#
  # Add core state-model monitors
  monitors <- c("psi.fs", "Trend", "alpha.s")

  #---------------------------------------------------#
  # Add detection intercepts: alpha.p
  monitors <- c(monitors, paste0("alpha.p", seq_len(n_datasets)))

  flatten_cov_names <- function(lst) {
    if (is.null(lst)) return(character())
    # lst is either a list of named vectors or list of lists.
    # unlist(names(.)) will give us exactly the cov names.
    covs <- unlist(lapply(lst, names), use.names = FALSE)
    # drop any empty strings
    covs[covs != ""]
  }

  #--------------------------------------------------#
  # State covariates (if we have any)
  state_covs <- flatten_cov_names(state_covariates_list)
  if (length(state_covs)) {
    monitors <- c(monitors, paste0("b_", state_covs))
  }

  # ------------------------------------------------#
  # Detection covs (if we have any)
  det_covs <- flatten_cov_names(det_covariates_list)
  if (length(det_covs)) {
    monitors <- c(monitors, paste0("b_", det_covs))
  }

  #------------------------------------------------
  # Phenology parameters, if used
  if (inclPhenology) {
    monitors <- c(monitors, "phScale", "phShape", "beta1")
  }

  return(monitors)
}

#' Extract covariate names
#'
#' Utility function to extract and clean covariate names from a list of covariate inputs.
#'
#' @param cov_list A named list of covariates.
#' @param format Output format, either "vector" or "list".
#'
#' @return A vector or list of cleaned covariate names.
#' @keywords internal
get_cov_names <- function(cov_list,
                          format = c("vector", "list")) {
  output_format <- match.arg(format)
  clean <- function(x) sub("^[^.]*\\.", "", x)
  if (output_format == "vector") {
    # flatten one level, keep names
    covs <- unlist(cov_list, recursive = FALSE, use.names = TRUE)
    nm <- clean(names(covs))
    return(unique(nm))
  } else if (output_format == "list") {
    # for each element of cov_list, return cleaned names
    out <- lapply(cov_list, function(sub) {
      nms <- names(sub)
      if (is.null(nms)) {
        character(0)
      } else {
        clean(nms)
      }
    })
    if (!is.null(names(cov_list))) {
      names(out) <- names(cov_list)
    }
    return(out)
  }
}

#' Build state-level covariate inputs
#'
#' Given a full list of pre-processed state covariates (site, time, or spatiotemporal),
#' this function extracts and organises them into structured lists for inclusion in the
#' NIMBLE iSDM. It checks for variable presence, resolves `.value` vs. matrix formats,
#' and handles polynomial terms where the covariate name ends in `2`.
#'
#' Polynomial character strings (e.g., `"tas_MATcs2"`) are resolved by stripping the trailing `2`
#' and accessing the `value2` slot in the corresponding covariate list.
#'
#' @param all_state_covariates A named list of covariate objects (usually the output of
#'   \code{\link{format_state_covariates}}). Each object can be:
#'   - a list with a `$value` and optional `$value2` (for polynomials), or
#'   - a matrix (for `it` covariates).
#' @param covariates_i Character vector of spatial covariate names to include.
#' @param covariates_t Character vector of temporal covariate names to include.
#' @param covariates_it Character vector of spatiotemporal covariate names to include.
#'
#' @return A named list with elements:
#' \describe{
#'   \item{covariates_i}{Named list of spatial covariate vectors.}
#'   \item{covariates_t}{Named list of temporal covariate vectors.}
#'   \item{covariates_it}{Named list of spatiotemporal covariate matrices.}
#' }
#'
#' @seealso \code{\link{format_state_covariates}}
#' @export
build_state_covariates <- function(all_state_covariates,
                                   covariates_i  = character(),
                                   covariates_t  = character(),
                                   covariates_it = character()) {
  # -------------------------------------------------------------------#
  # helper to pull out the base covariate vector
  pull_cov <- function(covname) {
    if (! covname %in% names(all_state_covariates)) {
      available <- sort(names(all_state_covariates))
      bullets <- paste0(" • ", available)
      stop(sprintf("Unknown covariate ‘%s’.
                   \nNOTE1: Variables are not calculated on-the-fly and must already be available in the input data.
                   \nNOTE2: Polynomial character strings are constructed in the function and exist in the data as list$covname$value2.
                   \nAvailable names are:\n%s",
                   covname, paste(bullets, collapse = "\n")
      )
      )
    }
    elt <- all_state_covariates[[covname]]
    if (is.list(elt) && !is.null(elt$value)) {
      return(elt$value)
    }
    # else: it must be already the matrix
    if (! is.matrix(elt)) {
      stop("Covariate ‘", covname,
           "’ is neither a list with $value nor a matrix.")
    }
    elt
  }

  # -------------------------------------------------------------------#
  # helper to pull out the polynomial term(s)
  pull_cov2 <- function(covname2) {
    base <- sub("2$", "", covname2)
    if (! base %in% names(all_state_covariates)) {
      stop("Requested polynomial ‘", covname2,
           "’ but base covariate ‘", base, "’ not found.
           \nPolynomial character strings are constructed in the function and should exist in the data as list$covname$value2.")
    }
    elt <- all_state_covariates[[base]]
    if (! is.list(elt) || is.null(elt$value2)) {
      stop("No second‐order (value2) component for ‘", base, "’.")
    }
    elt$value2
  }

  # -------------------------------------------------------------------#
  # build site‐only covs
  cov_i <- lapply(covariates_i, function(covname) {
    if (grepl("2$", covname)) pull_cov2(covname) else pull_cov(covname)
  })
  names(cov_i) <- covariates_i

  # -------------------------------------------------------------------#
  # build time‐only covs
  cov_t <- lapply(covariates_t, function(covname) {
    if (grepl("2$", covname)) pull_cov2(covname) else pull_cov(covname)
  })
  names(cov_t) <- covariates_t

  # -------------------------------------------------------------------#
  # build 2D covs
  cov_it <- lapply(covariates_it, function(covname) {
    if (! covname %in% names(all_state_covariates)) {
      available <- sort(names(all_state_covariates))
      bullets <- paste0(" • ", available)
      stop(sprintf("Unknown covariate ‘%s’.
                   \nNOTE1: Variables are not calculated on-the-fly and must already be available in the input data.
                   \nNOTE2: Polynomial character strings are constructed in the function and should exist in the data as list$covname$value2.
                   \nAvailable names are:\n%s",
                   covname, paste(bullets, collapse = "\n")
      )
      )
    }
    elt <- all_state_covariates[[covname]]
    if (is.list(elt) && !is.null(elt$value)) {
      mat <- elt$value
    } else {
      mat <- elt
    }
    if (! is.matrix(mat)) {
      stop("Spatio‐temporal covariate ‘", covname,
           "’ must be provided as a matrix or a list($value) that is a matrix.")
    }
    mat
  })
  names(cov_it) <- covariates_it

  # -------------------------------------------------------------------#
  # Return
  list(covariates_i = cov_i,
       covariates_t = cov_t,
       covariates_it = cov_it)
}

#' Format detection covariates for use in NIMBLE models
#'
#' This function processes a list of detection covariate datasets by optionally imputing missing values,
#' scaling the covariates, and renaming each column to a unique name for NIMBLE compatibility.
#'
#' @param covariates_det_list A named list of data frames. Each data frame should contain one or more detection covariates
#'   for a species or dataset. Names of the list will be used in diagnostic messages.
#' @param impute_missing Logical. If `TRUE`, missing values in each covariate are imputed using the mean of the non-missing values. Default is `FALSE`.
#' @param scale_cov Logical. If `TRUE`, all covariates are scaled using `scale()`. Default is `TRUE`.
#' @param center Logical. Passed to `scale()`. If `TRUE`, variables are centered before scaling. Default is `FALSE`.
#'
#' @return A list of lists. Each sub-list contains numeric vectors named by covariate and dataset index
#'   (e.g., `LL1`, `mins_on_site2`, etc.), ready for use in NIMBLE model input.
#'
#' @examples
#' # Example with two datasets, each with 2 covariates
#' cov_list <- list(
#'   species1 = data.frame(temp = c(20, 22, NA), effort = c(5, 7, 8)),
#'   species2 = data.frame(temp = c(18, 19, 21), effort = c(6, 7, 7))
#' )
#' formatted <- format_det_covariates(cov_list, impute_missing = TRUE)
#'
#' @export
format_det_covariates <- function(covariates_det_list,
                                  impute_missing = FALSE,
                                  scale_cov = TRUE,
                                  center = FALSE) {
  # -------------------------------------------------------------------#
  # Empty output list
  out <- vector("list", length(covariates_det_list))
  names(out) <- names(covariates_det_list)
  if (is.null(names(out))) names(out) <- paste0("ds", seq_along(out))

  # Loop over each data-frame
  for (i in seq_along(covariates_det_list)) {
    df <- covariates_det_list[[i]]
    ds_name <- names(covariates_det_list)[i]
    cov_names <- names(df)
    formatted <- list()

    for (cov in cov_names) {
      vec <- df[[cov]]
      # impute if requested
      if (impute_missing) {
        n_na <- sum(is.na(vec))
        if (n_na > 0) {
          mean_val <- mean(vec, na.rm = TRUE)
          vec[is.na(vec)] <- mean_val
          message(sprintf(
            "Imputed %d missing values in dataset '%s' for covariate '%s'",
            n_na, ds_name, cov))
        }
      }
      # scale if requested
      if (scale_cov) {
        vec <- as.numeric(scale(vec, center = center))
      }
      # store under names like "LL1", "mins_on_site2", etc
      formatted[[ paste0(cov, i) ]] <- vec
    }
    out[[i]] <- formatted
  }

  out
}

#' Build primary input list for reptimod iSDM
#'
#' Assembles all components needed to run the NIMBLE occupancy-detection model, including constants, data, initial values, arguments, and monitors.
#'
#' @param species_samples A list of formatted species detection data per species.
#' @param metadata Output from `get_metadata()`, including peakDay if phenology is used.
#' @param state_covariates Optional list of state covariates (e.g. site/year/site×year variables).
#' @param det_covariates Optional list of detection covariates.
#' @param inclPhenology Logical indicating whether to include phenology in the model.
#'
#' @return A list with named elements: constants, data, inits, args, and monitors.
#' @export
build_imod <- function(species_samples,
                       metadata,
                       state_covariates,
                       det_covariates,
                       inclPhenology=FALSE) {

  ## Define model constants list --------------
  constants <- build_model_constants(species_samples,
                                     state_covariates, det_covariates,
                                     inclPhenology = inclPhenology)

  ## Define model data list --------------------------
  data <- build_model_data(species_samples)

  ## Define model inits list----------------------
  inits <- build_model_inits(species_samples,
                             state_covariates, det_covariates,
                             peakDay = metadata$peakDay,
                             inclPhenology = inclPhenology)

  ## defineModel arguments-------------
  args <- build_model_arguments(species_samples,
                                state_covariates,
                                det_covariates,
                                inclPhenology = inclPhenology)

  ## Parameters to monitor----
  monitors <- build_monitors(state_covariates, det_covariates,
                             args$n_datasets, args$inclPhenology)

  # Remove phenology parameters from monitors if we have no inits for them
  if (!("beta1" %in% names(inits))) {
    monitors <- setdiff(monitors, c("beta1", "phShape", "phScale"))
  }

  return(list(constants = constants,
              data = data,
              inits = inits,
              args = args,
              monitors = monitors))
}

