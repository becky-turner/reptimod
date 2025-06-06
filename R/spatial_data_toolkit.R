# format_spatial_data.R
#
# Environmental covariate formatting utilities for iSDMs
# -------------------------------------------------------
# This script provides a suite of helper functions for handling, filtering,
# transforming, and visualising environmental covariate data to be used in
# integrated species distribution models (iSDMs) built in NIMBLE.

#' @importFrom magrittr %>%
#' @importFrom dplyr filter select left_join anti_join rename distinct arrange all_of any_of
#' @importFrom readr read_csv
#' @importFrom sf st_as_sf st_distance st_intersects
#' @importFrom GGally ggpairs wrap
#' @importFrom ggplot2 aes theme theme_minimal element_rect element_text
#' @importFrom reshape2 acast
#' @importFrom BRCmap gr_let2num gr_num2let reformat_gr OSGridstoLatLong
NULL

utils::globalVariables(c(
  "x", "y", "tas_mean", "tas_mean_CS", "woodcs", "heathercs", "UKelvcs",
  "latcs", "geometry", "survey_siteGrid", "value", "lc__bl_wood",
  "lc__heather", "UKelv", "climate_data_xy", "get_monad", "CRS", "str"
))

#' Create lookup of valid GB monads
#'
#' Loads background environmental data and applies the `get_monads()` function to extract the spatial lookup.
#'
#' @param env_filepath Path to the `.rds` file containing the environmental raster stack.
#' @param env_var_name Name of the variable to extract from the list (default is "tas_MAT").
#'
#' @return A data frame of 1km monads from the specified environmental variable.
#'
#' @keywords internal
valid_gb_monads <- function(env_filepath, env_var_name = "tas_MAT") {
  climate_data <- readRDS(env_filepath)
  env_var <- climate_data[[env_var_name]]
  env_monads_lookup <- get_monads(env_var)
  rm(climate_data, env_var)
  return(env_monads_lookup)
}

#' Filter species datasets by valid GB monads
#'
#' Filters each dataset in the input list to retain only those sites for which valid environmental data exists.
#'
#' This ensures compatibility between species records and environmental predictors.
#'
#' @param species_datasets_list A named list of species datasets (data frames).
#' @param env_filepath Path to the environmental `.rds` file. Default is ERA5 winter temperature.
#' @param env_var_name Name of the variable in the RDS to use. Default is "tas_MAT".
#' @param env_value The column in the environmental dataset to filter by (default = "tas_mean").
#' @param region_filepath Path to a CSV file listing valid grid references to retain.
#' @param gr_col Column name in the region file with grid references. Default is "grid_ref".
#' @param dims Number of dimensions in environmental data (1 = spatial).
#' @param matching Whether to match nearest grid cell if exact match not available.
#'
#' @importFrom magrittr %>%
#'
#' @return A filtered list of species datasets with only valid `x1km_grid` entries retained.
#'
#' @export
filter_valid_sites <- function(species_datasets_list,
                               env_filepath = "/data/env-data/ERA5-land-tas/winter_tas_2007to2020_2025_03_08.rds",
                               env_var_name = "tas_MAT",
                               env_value = "tas_mean",
                               region_filepath = "/data/env-data/gr_ref.csv",
                               gr_col = "grid_ref",
                               dims = 1, matching = FALSE) {

  env_monads_lookup <- valid_gb_monads(env_filepath = env_filepath, env_var_name = env_var_name)

  filtered_species_list <- lapply(species_datasets_list, function(species_df) {
    valid_sites0 <- filter_site_auxillary(data = env_monads_lookup,
                                          values = env_value,
                                          sitecodes = species_df,
                                          yearcodes = NULL,
                                          dims = dims,
                                          matching = matching)

    if (!is.null(region_filepath)) {
      region_file <- file.path(region_filepath)
      if (!file.exists(region_file)) {
        stop("Cannot find region file: ", region_file)
      }
      region_df <- readr::read_csv(region_file, na = "")
      if (! gr_col %in% names(region_df)) {
        stop("Column '", gr_col, "' not found in your region file.")
      }
      keep_refs <- region_df[[gr_col]]
      valid_sites0 <- dplyr::filter(valid_sites0, x1km_grid %in% keep_refs)
    }

    valid_sites <- valid_sites0$x1km_grid
    species_filtered <- species_df %>% dplyr::filter(x1km_grid %in% valid_sites)
    return(species_filtered)
  })

  return(filtered_species_list)
}

#' Find nearest grid cell by euclidean distance
#'
#' Matches each site in `data1` to the nearest grid square in `data2`, using OSGB eastings/northings.
#'
#' @param data1 A data frame of unmatched sites with column `x1km_grid`.
#' @param data2 A data frame of reference sites with column `x1km_grid`.
#'
#' @return A data frame of the nearest matching grid rows from `data2`.
#' @importFrom sf st_as_sf st_distance st_intersects
#' @importFrom dplyr distinct select
#' @export
find_nearest_grid <- function(data1, data2){

  # if no missing sites to match, quit function:
  if (nrow(data1) == 0) {
    # data2[0,] has the same columns but zero rows
    return(data2[0, ])
  }

  #get easting/northing
  easnordat <- BRCmap::gr_let2num(data1$x1km_grid)

  data1 <- cbind(data1, easnordat) # data1 = sampling sites data

  data2 <- data2 %>% distinct(x1km_grid, .keep_all = TRUE)

  easnordata2 <- BRCmap::gr_let2num(data2$x1km_grid)

  data2 <- cbind(data2, easnordata2) # data2 = basemap data (e.g., temperature raster)

  # Convert Eastings and Northings to coordinates (assuming the datasets are in OSGB projection)

  missingsite_coords <- st_as_sf(data1, coords = c("EASTING", "NORTHING"), crs = 27700)

  basemap_coords <- st_as_sf(data2, coords = c("EASTING", "NORTHING"), crs = 27700)

  # Calculate distances between each pair of grid squares
  distances <- st_distance(missingsite_coords$geometry, basemap_coords$geometry)
  str(distances)

  # Find the closest grid square in data2 for each grid square in data1
  closest_indices <- apply(distances, 1, which.min)
  str(closest_indices)

  # Check for overlap between boundaries
  overlaps <- st_intersects(basemap_coords$geometry, missingsite_coords$geometry)

  # Print the closest grid squares and overlapping grid squares

  closest_basemap_grids <- basemap_coords[closest_indices, ]

  closest_basemap_grids <- closest_basemap_grids %>% as.data.frame() %>% select(-geometry)

  return(closest_basemap_grids)
}

#' Convert eastings/northings to OSGB 1km grid squares
#'
#' Uses `BRCmap` functions to convert x/y coordinates to 1km OS grid references (monads).
#'
#' @param data A data frame with columns `x` and `y` (eastings/northings).
#' @param crs Character. Currently only "OSGB" supported.
#'
#' @return A data frame with new `x1km_grid` column.
#' @importFrom BRCmap gr_num2let reformat_gr
#' @importFrom dplyr distinct
#' @export
get_monads <- function(data, crs="OSGB"){

  gridrefs <- BRCmap::gr_num2let(data$x, data$y,
                                 OSgrid="OSGB",
                                 keep_precision = TRUE,
                                 min_10km = FALSE)
  data$x1km_grid <- BRCmap::reformat_gr(gridrefs, prec_out = 1000)

  # Check if `year` exists in the dataset
  if ("year" %in% colnames(data)) {
    data <- data %>% distinct(x1km_grid, year, .keep_all = TRUE)  # Keep all years
  } else {
    data <- data %>% distinct(x1km_grid, .keep_all = TRUE)  # Only filter by grid
  }

  return(data)

  if(CRS != "OSGB") {

    print("Error! Only eastings and northings be re-projected")

  }

}

#' Filter environmental data to match species survey sites
#'
#' Filters auxiliary environmental data to match the survey site grid squares, optionally by year.
#'
#' @param data Environmental data with `x1km_grid` (and `year` if dims=2).
#' @param values Column name of variable to retain.
#' @param sitecodes Data frame with `x1km_grid` column.
#' @param yearcodes Optional. Data frame with `year` and `yearID` columns.
#' @param dims 1 for site-only data, 2 for site-year data.
#' @param matching Logical. If TRUE, calls `match_nearest_env()` to fill gaps.
#'
#' @return Filtered environmental data.
#' @importFrom dplyr filter select left_join
#' @export
filter_site_auxillary <- function(data, values, sitecodes, yearcodes = NULL, dims = 1, matching=FALSE){
  # Ensure x1km_grid exists
  if (!"x1km_grid" %in% colnames(data) | !"x1km_grid" %in% colnames(sitecodes)) {
    stop("Both the environmental dataset and sitecodes dataset must contain 'x1km_grid' column.")
  }

  # Ensure values exists in dataset
  if (!values %in% colnames(data)) {
    stop(paste("Error in filter_site_auxillary(). Column", values, "not found in dataset."))
  }

  # **Handle 1D datasets (site)**
  if (dims == 1) {
    DataxSurveySites <- data %>%
      filter(x1km_grid %in% sitecodes$x1km_grid) %>%
      select(x1km_grid, all_of(values))

    # **Handle 2D datasets (site-year)**
  } else if (dims == 2) {
    if (is.null(yearcodes)) {
      stop("Error: yearcodes must be provided for dims='2'.")
    }

    if (!"year" %in% colnames(yearcodes)) {
      stop("yearcodes dataset must contain 'year' column for dims='2'.")
    }

    DataxSurveySites <- data %>%
      filter(x1km_grid %in% sitecodes$x1km_grid, year %in% yearcodes$year) %>%
      select(x1km_grid, year, all_of(values))

    # Match `year` to `yearID`
    DataxSurveySites <- DataxSurveySites %>%
      left_join(yearcodes, by = "year") %>%
      select(-year) %>%
      select(x1km_grid, yearID, all_of(values))
  } else {
    stop("Error: dims argument must be either '1' or '2'.")
  }

  # **Optional: Find nearest available environmental data relative to a given location**
  if (matching == TRUE) {
    print("Matching missing x1km_grid survey sites to nearest available x1km_grid-year square with environmental data.")
    DataxSurveySites <- match_nearest_env(data, values, sitecodes, yearcodes, dims)
  } else {
    print("Filtered auxiliary data and removed sites without matching environmental x1km_grid squares.")
  }

  return(DataxSurveySites)
}

#' Match missing environmental data nearest available OSGB square
#'
#' Identifies missing x1km grids (and optionally years) and substitutes the nearest available environmental information.
#'
#' @param data Environmental data.
#' @param values Variable to extract from data.
#' @param sitecodes Data frame of site names / grid references.
#' @param yearcodes Optional year names to match.
#' @param dims 1 (site) or 2 (site-year).
#'
#' @return An environmental dataset matched to survey sites and/or years.
#' @importFrom dplyr filter select left_join anti_join rename
#' @export
match_nearest_env <- function(data, values, sitecodes, yearcodes = NULL, dims = 1){
  print("Matching x1km_grids and years to environmental data where needed")

  # **1D: Match by x1km_grid only**
  if (dims == 1) {
    temp <- data %>%
      filter(x1km_grid %in% sitecodes$x1km_grid) %>%
      select(x1km_grid, all_of(values))

    # Find missing sites
    missing_sites <- sitecodes %>%
      filter(!x1km_grid %in% temp$x1km_grid) %>%
      distinct(x1km_grid)

    # Match missing sites
    subbedGrids <- find_nearest_grid(missing_sites, data)
    subbedGrids$survey_siteGrid <- missing_sites$x1km_grid
    subbedGrids <- subbedGrids %>%
      select(survey_siteGrid, all_of(values)) %>%
      rename(x1km_grid = survey_siteGrid)

    # Combine matched and substituted environmental data
    DataxSurveySites <- rbind(temp, subbedGrids)

    # **2D: Match by x1km_grid & year**
  } else if (dims == 2) {
    if (is.null(yearcodes)) {
      stop("Error: yearcodes must be provided for dims='2'.")
    }

    # Build the expected site-year combinations from sitecodes and yearcodes
    # (Assuming sitecodes has a column 'x1km_grid' and yearcodes has a column 'year')
    expected <- expand.grid(x1km_grid = unique(sitecodes$x1km_grid),
                            year = unique(yearcodes$year),
                            stringsAsFactors = FALSE)

    # Get the actual environmental data for these combinations
    actual <- data %>%
      filter(x1km_grid %in% unique(sitecodes$x1km_grid),
             year %in% unique(yearcodes$year)) %>%
      select(x1km_grid, year, all_of(values))

    # Identify the missing site-year combinations
    missing <- anti_join(expected, actual, by = c("x1km_grid", "year"))

    # For each missing combination, use find_nearest_grid() to substitute data
    if(nrow(missing) > 0) {
      subbedGrids <- find_nearest_grid(missing, data)
      # Attach the missing combination's x1km_grid and year to the substituted data
      subbedGrids$survey_siteGrid <- missing$x1km_grid
      subbedGrids$year <- missing$year
      subbedGrids <- subbedGrids %>%
        select(survey_siteGrid, year, all_of(values)) %>%
        rename(x1km_grid = survey_siteGrid)

      # Convert the 'year' to a yearID by joining with yearcodes
      subbedGrids <- subbedGrids %>%
        left_join(yearcodes, by = "year") %>%
        select(-year) %>%
        select(x1km_grid, yearID, all_of(values))
    } else {
      # If none missing, create an empty data frame with the correct columns
      subbedGrids <- data.frame(x1km_grid = character(), yearID = integer(),
                                stringsAsFactors = FALSE)
      subbedGrids[[values]] <- numeric()
    }

    # Convert the actual environmental data to include yearID:
    actual_joined <- actual %>%
      left_join(yearcodes, by = "year") %>%
      select(-year) %>%
      select(x1km_grid, yearID, all_of(values))

    # Combine the actual data and the substituted data
    DataxSurveySites <- rbind(actual_joined, subbedGrids)

    return(DataxSurveySites)

  }

}

#' Format environmental covariate data for NIMBLE (1D, site)
#'
#' Matches grid references to siteID and returns a site-level covariate vector.
#'
#' @param data Environmental data with `x1km_grid`.
#' @param dataset_name Name of the covariate dataset.
#' @param values Name of variable to extract.
#' @param sitecodes Data frame with `x1km_grid` and `siteID`.
#'
#' @return A data frame with columns `siteID` and `value`.
#' @export
nimblefy_i <- function(data, dataset_name, values, sitecodes){

  print(colnames(data))
  print("Matching to nimble sitecodes")
  # Match x1km_grid to NIMBLE site codes
  result <- data %>%
    left_join(sitecodes, by = "x1km_grid") %>%
    rename(value = all_of(values)) %>%
    select(siteID, value) %>%
    arrange(siteID)

  # Return the result with as named list
  return(result)
  #return(list(result))
  #return(setNames(list(result), dataset_name))

}

#' Format environmental covariate data for NIMBLE (1D, year)
#'
#' Matches years to yearID and returns a year-level covariate vector.
#'
#' @param data Environmental data with `year`.
#' @param dataset_name Name of the covariate dataset.
#' @param values Name of variable to extract.
#' @param yearcodes Data frame with `year` and `yearID`.
#'
#' @return A data frame with columns `yearID` and `value`.
#' @export
nimblefy_t <- function(data, dataset_name, values, yearcodes){

  print("Matching to nimble yearcodes")

  # Match x1km_grid to NIMBLE site codes
  result <- data %>%
    left_join(yearcodes, by = "year") %>%
    rename(value = all_of(values)) %>%
    select(yearID, value) %>%
    arrange(yearID)

  # Return the result with as named list
  return(result)

}

#' Format environmental covariate data for NIMBLE (2D, site/year)
#'
#' Reshapes a site-by-year covariate matrix from long to wide using siteID and yearID.
#'
#' @param data Environmental data with `x1km_grid` and `year`.
#' @param dataset_name Name of the covariate dataset.
#' @param values Name of variable to extract.
#' @param sitecodes Data frame with `x1km_grid` and `siteID`.
#' @param yearcodes Data frame with `year` and `yearID`.
#'
#' @return A matrix with dimensions siteID × yearID.
#' @importFrom reshape2 acast
#' @export
nimblefy_it <- function(data, dataset_name, values, sitecodes, yearcodes){

  print("Matching x1km_grids to NIMBLE sitecodes and formatting for 2D site/year structure")
  # Match site codes
  data <- data %>%
    left_join(sitecodes, by = "x1km_grid") %>%
    rename(value = all_of(values)) %>%
    select(-x1km_grid)

  # Convert to siteID x yearID matrix
  data_matrix <- reshape2::acast(data, siteID ~ yearID, value.var = "value", fill = 0, fun = max)

  return(data_matrix)

}

#' Wrapper to prepare environmental covariates for NIMBLE
#'
#' Handles filtering, matching, and formatting of covariate data for site, year, or site-year structures.
#'
#' @param data_list Named list of environmental data frames.
#' @param dataset_name Name of the target dataset within `data_list`.
#' @param values Name of variable to extract.
#' @param sitecodes Data frame with `x1km_grid` and `siteID`.
#' @param yearcodes Data frame with `year` and `yearID`.
#' @param dims Either 1 or 2 for site/year or site-year data.
#' @param matching Logical. Whether to use nearest-match fallback.
#' @param nimblefy_method "i", "t", or "it" for site, time, or site-year structure.
#'
#' @return Covariate formatted for input to NIMBLE.
#' @export
prep_nimble_cov <- function(data_list, dataset_name, values,
                            sitecodes=NULL, yearcodes = NULL,
                            dims = 1, matching=FALSE,
                            nimblefy_method = "i"){

  print("Starting preprocessing for NIMBLE")

  # **Step 1: Check inputs**
  if (!dataset_name %in% names(data_list)) {
    stop(paste("Dataset", dataset_name, "not found in environmental data list object."))
  }

  if (!values %in% colnames(data_list[[dataset_name]])) {
    stop(paste("Error: Column", values, "not found in dataset", dataset_name))
  }

  print("Datasets and columns passed naming checks.")

  # **Step 2: Extract dataset**
  data <- data_list[[dataset_name]]

  if (nimblefy_method != "t") {

    # **Step 3: Filter & match auxiliary data**
    data <- filter_site_auxillary(data, values, sitecodes, yearcodes, dims, matching)

    # **Step 4: Format for NIMBLE**
    if (nimblefy_method == "i") {
      nimble_cov <- nimblefy_i(data, dataset_name, values, sitecodes)
    } else if (nimblefy_method == "it") {
      nimble_cov <- nimblefy_it(data, dataset_name, values, sitecodes, yearcodes)
    } else {
      stop("Error: nimblefy_method must be 'i' or 'it'.")
    }

    print("Environmental data has been prepped for NIMBLE model.")

  } else {
    if (is.null(yearcodes)) {
      stop(paste("Error: yearcodes = NULL. Must supply yearcodes"))
    }

    nimble_cov <- nimblefy_t(data, dataset_name, values, yearcodes)

  }

  return(nimble_cov)

}

#' Import climate covariates
#'
#' Imports processed raster-derived climate covariates. Can extract 1km grid references if `get_monads = TRUE`.
#'
#' @param get_monads Logical. Whether to compute x1km_grid from coordinates.
#' @param region_filepath Optional CSV of grid refs to filter datasets.
#' @param gr_col Column in region file containing grid references.
#'
#' @return A list of climate datasets: `climate_nonspatial`, `climate_xy`, `soil_data`.
#' @importFrom readr read_csv
#' @importFrom dplyr filter
#' @export
import_climate_rasters <- function(get_monads = FALSE,
                                   region_filepath = NULL,
                                   gr_col = "grid_ref") {

  #### 1) Import raw climate data
  climate_data <- readRDS("/data/env-data/ERA5-land-tas/winter_tas_2007to2020_2025_03_08.rds")

  # **Extract Non-Spatial & Spatial Datasets**
  climate_nonspatial <- climate_data[sapply(climate_data, function(df) !all(c("x", "y") %in% colnames(df)))]

  if (get_monads == TRUE) {
    # Extract spatial datasets
    # tas
    climate_xy <- climate_data[sapply(climate_data, function(df) all(c("x", "y") %in% colnames(df)))]
    #soil
    soil_data <- readRDS("/data/env-data/soil_data/winter_soil_water_2007to2020_2025_03_08.rds")

    # Apply get_monad() to spatial datasets
    #tas
    climate_xy <- lapply(climate_data_xy, get_monad)
    #soil
    soil_data <- lapply(soil_data, get_monad)

  } else {
    # Load pre-processed climate data if monads are not needed
    climate_xy <- readRDS("/data/env-data/ERA5-land-tas/climate_data_xy_2007to2020_2025-03-09.rds")
    soil_data <- readRDS("/data/env-data/soil_data/winter_soil_water_2007to2020_2025_03_10.rds")

  }

  #-----------------------------------------------------------------#
  # If a region file is supplied, read its grid_refs and filter
  if(!is.null(region_filepath)) {
    region_file <- file.path(region_filepath)
    if(!file.exists(region_file)) {
      stop("Cannot find region file: ", region_file)
    }
    region_df <- readr::read_csv(region_file, na = "")
    if(! gr_col %in% names(region_df)) {
      stop("Column '", gr_col, "' not found in your region file.")
    }
    refs <- region_df[[gr_col]]

    # filter every spatial dataframe in climate_xy and soil_data
    climate_xy <- lapply(climate_xy, function(df) {
      df %>% filter(x1km_grid %in% refs)
    })
    soil_data <- lapply(soil_data, function(df) {
      df %>% filter(x1km_grid %in% refs)
    })
  }

  return(list(climate_nonspatial=climate_nonspatial,
              climate_xy=climate_xy,
              soil_data=soil_data))

}

#' Import and process landcover covariates
#'
#' Imports UKCEH landcover data and performs centre + scaling transformations.
#'
#' @param filepath Path to CSV file of landcover data.
#' @param region_filepath Optional CSV file with grid references for filtering.
#' @param gr_col Name of grid ref column in region file.
#'
#' @return A list with covariate data frames.
#' @export
import_landcover_dfs <- function(filepath = "/data/env-data/ukceh_landcover_all_aux_robB.csv",
                                 region_filepath = NULL,
                                 gr_col = "grid_ref") {

  # Import data
  landcover <- readr::read_csv(paste0(filepath))
  landcover <- landcover %>% rename("x1km_grid" = "monad")

  # if region filtering is requested, load the CSV of grid‐refs
  if (!is.null(region_filepath)) {
    region_file <- file.path(region_filepath)
    if (!file.exists(region_file)) {
      stop("Cannot find region file: ", region_file)
    }
    region_df <- readr::read_csv(region_file, na = "")
    if (! gr_col %in% names(region_df)) {
      stop("Column '", gr_col, "' not found in your region file.")
    }
    keep_refs <- region_df[[gr_col]]
    landcover <- filter(landcover, x1km_grid %in% keep_refs)
  }

  #------------------------------------------------------#
  # Broad-leaved woodland
  woodland <- landcover %>% select(x1km_grid, x, y, lc__bl_wood)
  # Centre & scale
  woodland$woodcs <- as.vector(scale(woodland$lc__bl_wood))

  #------------------------------------------------------#
  # Heathland
  heathland <- landcover %>% select(x1km_grid, x, y, lc__heather)
  # Centre & scale
  heathland$heathercs <- as.vector(scale(heathland$lc__heather))

  rm(landcover)

  return(list(woodland=woodland,
              heathland=heathland))

}

#' Import and process spatial properties (elevation, latitude)
#'
#' Loads and processes spatial covariates from landcover dataset, including elevation and latitude.
#'
#' @param filepath Path to input CSV.
#' @param region_filepath Optional path to grid-ref CSV for filtering.
#' @param gr_col Column name in region file containing grid references.
#'
#' @return A list with `elevation` and `latitude` data frames.
#' @export
import_spat_properties_df <- function(filepath = "/data/env-data/ukceh_landcover_all_aux_robB.csv",
                                      region_filepath = NULL,
                                      gr_col = "grid_ref") {

  # Import data
  landcover <- readr::read_csv(filepath)
  landcover <- landcover %>% rename("x1km_grid" = "monad")

  # if region filtering is requested, load the CSV of grid‐refs
  if (!is.null(region_filepath)) {
    region_file <- file.path(region_filepath)
    if (!file.exists(region_file)) {
      stop("Cannot find region file: ", region_file)
    }
    region_df <- readr::read_csv(region_file, na = "")
    if (! gr_col %in% names(region_df)) {
      stop("Column '", gr_col, "' not found in your region file.")
    }
    keep_refs <- region_df[[gr_col]]
    landcover <- filter(landcover, x1km_grid %in% keep_refs)
  }

  #------------------------------------------------------#
  # Elevation
  elevation <- landcover %>% select(x1km_grid, x, y, UKelv)
  # Centre & scale
  elevation$UKelvcs <- as.vector(scale(elevation$UKelv))

  #------------------------------------------------------#
  # Latitude
  latitude <- landcover %>% select(x1km_grid, x, y)
  # Calculate latitude
  latlon_dat <- BRCmap::OSGridstoLatLong(latitude$x, latitude$y)
  latitude <- cbind(latitude, latlon_dat)
  latitude <- latitude[,c(1:4)]
  latitude <- latitude %>% rename(lat = "LATITUDE")

  # Centre & scale
  latitude$latcs <- as.vector(scale(latitude$lat))

  rm(landcover)

  return(list(elevation=elevation, latitude=latitude))

}

#' Format state-level environmental covariates for NIMBLE
#'
#' Combines and processes spatial, temporal, and spatiotemporal covariate datasets
#' (climate, landcover, and spatial properties) to match model site and year structures
#' for use in integrated species distribution models (iSDMs) built in NIMBLE.
#'
#' The function pulls in pre-processed datasets using dedicated import functions, and
#' applies `prep_nimble_cov()` with `nimblefy_*()` helpers to match covariates to the
#' expected model dimensions (`i`, `t`, `it`). It also applies centering, scaling, and
#' optional second-order polynomial transformation.
#'
#' @param sitecodes A data frame with columns `x1km_grid` and `siteID` mapping sites to IDs.
#' @param yearcodes A data frame with columns `year` and `yearID` for temporal matching.
#' @param matching Logical. If TRUE, attempts to fill missing environmental records with
#'   the nearest available value.
#' @param types A character vector of which covariate types to process. Must be a subset of:
#'   `"i"` (spatial), `"t"` (temporal), `"it"` (spatiotemporal).
#'
#' @return A named list of covariate inputs formatted for inclusion in the NIMBLE model.
#'
#' @seealso \code{\link{prep_nimble_cov}}, \code{\link{import_climate_rasters}}, \code{\link{import_landcover_dfs}}, \code{\link{import_spat_properties_df}}
#' @export
format_state_covariates <- function(sitecodes, yearcodes, matching=FALSE, types = c()) {

  #------------------------------------------------------#
  # LOAD DATASETS

  # climate data: temperature, rainfall, soil
  climate <- import_climate_rasters(region_filepath = "/data/env-data/gr_ref.csv",
                                    gr_col = "grid_ref")

  # landcover data: woodland, heathland %
  landcover <- import_landcover_dfs(region_filepath = "/data/env-data/gr_ref.csv",
                                    gr_col = "grid_ref")

  # properties data: elevation, latitude
  properties <- import_spat_properties_df(region_filepath = "/data/env-data/gr_ref.csv",
                                          gr_col = "grid_ref")

  #------------------------------------------------------#
  # RUN COVARIATE FORMATTING

  # Combine all datasets into one "covs" list.
  # Ensure all variables are centred + scaled, and add polynomial term.

  covs <- list()

  #------------------------------------------------------#
  ## **1D TEMPORAL DATA**
  if ("t" %in% types) {
    covs$year_tas_anomaly <- prep_nimble_cov(data_list = climate$climate_nonspatial,
                                             dataset_name = "year_tas_anomaly",
                                             values = "mean_anomaly",
                                             yearcodes = yearcodes,
                                             dims = 1, matching=matching,
                                             nimblefy_method = "t")
  }

  #------------------------------------------------------#
  ## **1D SPATIAL DATA**
  if ("i" %in% types) {
    #temperature---
    covs$tas_MAT <- prep_nimble_cov(data_list = climate$climate_xy,
                                    dataset_name = "tas_MAT",
                                    values = "tas_mean", sitecodes = sitecodes,
                                    dims = 1, matching=matching,
                                    nimblefy_method = "i")

    covs$tas_MATcs <- prep_nimble_cov(data_list = climate$climate_xy, dataset_name = "tas_MATcs",
                                      values = "tas_mean_CS", sitecodes = sitecodes,
                                      dims = 1, matching=matching,
                                      nimblefy_method = "i")

    covs$tas_MATcs$value2 <- covs$tas_MATcs$value^2

    covs$tas_ManTA <- prep_nimble_cov(data_list = climate$climate_xy, dataset_name = "tas_ManTA",
                                      values = "mean_anTA", sitecodes = sitecodes,
                                      dims = 1, matching=matching,
                                      nimblefy_method = "i")

    #soil---
    covs$soil_watercs <- prep_nimble_cov(data_list = climate$soil_data, dataset_name = "soil_watercs",
                                         values = "mean", sitecodes = sitecodes,
                                         dims = 1, matching=matching,
                                         nimblefy_method = "i")

    covs$soil_watercs$value2 <- covs$soil_watercs$value^2

    #landcover---

    #woodland
    covs$woodlandcs <- prep_nimble_cov(data_list = landcover, dataset_name = "woodland",
                                       values = "woodcs", sitecodes = sitecodes,
                                       dims = 1, matching=TRUE,
                                       nimblefy_method = "i")
    covs$woodlandcs$value2 <- covs$woodlandcs$value^2

    #heathland
    covs$heathlandcs <- prep_nimble_cov(data_list = landcover, dataset_name = "heathland",
                                        values = "heathercs", sitecodes = sitecodes,
                                        dims = 1, matching=TRUE,
                                        nimblefy_method = "i")
    covs$heathlandcs$value2 <- covs$heathlandcs$value^2

    # spatial properties---

    #elevation
    covs$elevationcs <- prep_nimble_cov(data_list = properties, dataset_name = "elevation",
                                        values = "UKelvcs", sitecodes = sitecodes,
                                        dims = 1, matching=TRUE,
                                        nimblefy_method = "i")
    covs$elevationcs$value2 <- covs$elevationcs$value^2

    #latitude
    covs$latcs <- prep_nimble_cov(data_list = properties, dataset_name = "latitude",
                                  values = "latcs", sitecodes = sitecodes,
                                  dims = 1, matching=TRUE,
                                  nimblefy_method = "i")
    covs$latcs$value2 <- covs$latcs$value^2


  }

  #------------------------------------------------------#
  ## **2D SPATIOTEMPORAL DATA**

  if ("it" %in% types) {
    covs$annual_MAT <- prep_nimble_cov(data_list = climate$climate_xy, dataset_name = "annual_MAT",
                                       values = "tas_MAT",
                                       sitecodes = sitecodes, yearcodes = yearcodes,
                                       dims = 2, matching=matching, nimblefy_method = "it")

    covs$annual_MAT_dev <- prep_nimble_cov(data_list = climate$climate_xy, dataset_name = "annual_MAT_dev",
                                           values = "tas_dev",
                                           sitecodes = sitecodes, yearcodes = yearcodes,
                                           dims = 2, matching=matching, nimblefy_method = "it")

    covs$tas_anTA <- prep_nimble_cov(data_list = climate$climate_xy, dataset_name = "tas_anTA",
                                     values = "mean_anTA",
                                     sitecodes = sitecodes, yearcodes = yearcodes,
                                     dims = 2, matching=matching, nimblefy_method = "it")
  }

  return(covs)
}

#' Pairwise Plot of 1D State Covariates
#'
#' Produces a matrix of pairwise plots for all one-dimensional covariates
#' (spatial and temporal) contained within a `state_covariates` object.
#' Correlations are displayed in the upper triangle and scatterplots in the lower
#' triangle (or either, as selected).
#'
#' Covariates are drawn from the `covariates_i` and `covariates_t` slots of the
#' input list. Covariates must be stored as named vectors.
#'
#' @param state_covariates A list of covariates (typically the output of
#'   \code{\link{build_state_covariates}}), containing elements \code{$covariates_i}
#'   and \code{$covariates_t}.
#' @param panels Character string indicating which panels to show:
#'   \code{"both"} (default), \code{"lower"}, or \code{"upper"}.
#'
#' @return A \code{ggpairs} object for plotting.
#'
#' @seealso \code{\link{build_state_covariates}}, \code{\link[GGally]{ggpairs}}
#'
#' @importFrom GGally ggpairs wrap
#' @importFrom ggplot2 aes theme_minimal element_rect element_text
#' @export
plot_state_pairs <- function(state_covariates,
                             panels = c("both")) {
  panels <- match.arg(panels, c("lower","upper","both"))
  # pull out all the 1-D covariates
  vecs <- c(state_covariates$covariates_i,
            state_covariates$covariates_t)
  if (length(vecs) == 0) {
    stop("No 1D covariates found in state_covariates$covariates_i or $covariates_t")
  }
  # make a data.frame where each list-element becomes one column
  df <- as.data.frame(vecs)
  names(df) <- names(vecs)

  # build the ggpairs call
  low_list <- if (panels %in% c("lower","both")) {
    list(continuous = GGally::wrap("points",
                                   shape = 1,
                                   color = "steelblue4",
                                   alpha=1))
  } else NULL
  up_list <- if (panels %in% c("upper","both")) {
    list(continuous = GGally::wrap("cor", size = 4))
  } else NULL

  diag_list <- list(continuous = GGally::wrap("densityDiag"))

  GGally::ggpairs(df,
                  mapping = aes(),
                  lower = low_list,
                  upper = up_list,
                  diag = diag_list
  ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      strip.background = element_rect(fill = "steelblue", color = NA),
      strip.text = element_text(colour = "white")
    )
}

#' Drop specified covariate columns from a list of data frames
#'
#' Removes unwanted columns from each element (data frame) in a list of covariate datasets.
#' Useful for pruning unused or redundant covariates prior to formatting or modelling.
#'
#' @param cov_list A list of data frames, each containing covariate columns.
#' @param cols A character vector of column names to drop from each data frame in the list.
#'
#' @return A list of data frames with the specified columns removed.
#'
#' @examples
#' covs <- list(
#'   df1 = data.frame(a = 1:3, b = 4:6, c = 7:9),
#'   df2 = data.frame(a = 10:12, b = 13:15, c = 16:18)
#' )
#' cleaned <- drop_covs(covs, cols = c("b", "c"))
#'
#' @export
drop_covs <- function(cov_list, cols = c()) {
  ### Remove any unwanted columns
  lapply(cov_list, function(df) {
    df %>% select(-any_of(cols))
  })
}

#' Build map layers for selected spatial covariates
#'
#' Loads, prepares, and standardizes spatial environmental datasets
#' (e.g. temperature, land cover, elevation) for mapping. Returns a list of
#' data frames each containing `x1km_grid`, `x`, `y`, and `value` columns.
#'
#' Optionally handles polynomial covariates by squaring the base variable.
#'
#' @param covariates A character vector of covariate names to extract.
#'                   Polynomial terms can be specified by appending '2' (e.g., `"latcs2"`).
#'
#' @return A named list of data frames, each with columns `x1km_grid`, `x`, `y`, and `value`,
#'         where `value` is the selected covariate or its squared value.
#'
#' @details
#'
#' Covariates are sourced from:
#' - `import_climate_rasters()`
#' - `import_spat_properties_df()`
#' - `import_landcover_dfs()`
#'
#' @examples
#' \dontrun{layers <- build_map_layers(covariates = c("tas_MATcs", "latcs2"))}
#'
#' @export
build_map_layers <- function(covariates = character()) {
  # 1) load environmental inputs
  clim  <- import_climate_rasters(region_filepath = "/data/env-data/gr_ref.csv",
                                  gr_col = "grid_ref")
  spat <- import_spat_properties_df(region_filepath = "/data/env-data/gr_ref.csv",
                                    gr_col = "grid_ref")
  lc <- import_landcover_dfs(region_filepath = "/data/env-data/gr_ref.csv",
                             gr_col = "grid_ref")

  # 2) pull into one named list, always renaming to "value"
  spatial_env <- list(
    # climate_xy
    tas_MAT = clim$climate_xy$tas_MAT %>%
      select(x1km_grid, x, y, value = tas_mean),
    tas_MATcs = clim$climate_xy$tas_MATcs %>%
      select(x1km_grid, x, y, value = tas_mean_CS),
    # soil
    soil_watercs = clim$soil_data$soil_watercs %>%
      select(x1km_grid, x, y, value = mean),
    # landcover
    woodlandcs  = lc$woodland %>% select(x1km_grid, x, y, value = woodcs),
    heathlandcs = lc$heathland %>% select(x1km_grid, x, y, value = heathercs),
    # elevation/lat
    elevationcs = spat$elevation %>% select(x1km_grid, x, y, value = UKelvcs),
    latcs = spat$latitude %>% select(x1km_grid, x, y, value = latcs)
  )

  out <- list()
  for(cov in covariates) {
    if(grepl("2$", cov)) {
      # polynomial: strip the trailing "2"
      base <- sub("2$", "", cov)
      if(! base %in% names(spatial_env))
        stop("No base covariate called '", base, "' for polynomial '", cov, "'")
      df <- spatial_env[[base]]
      df$value <- df$value^2
      out[[cov]] <- df
    } else {
      if(! cov %in% names(spatial_env))
        stop("Unknown covariate '", cov, "'")
      out[[cov]] <- spatial_env[[cov]]
    }
  }
  out
}

