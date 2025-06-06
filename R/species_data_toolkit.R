# format_species_data.R
#
# Functions for preparing and validating input data for integrated species distribution models (iSDMs)
# using the reptimod workflow. This includes:
# - Importing pre-cleaned species datasets
# - Filtering observations based on environmental data availability
# - Summarising detection histories
# - Creating pooled datasets

#' Import pre-processed species datasets
#'
#' Reads all `.rds` files from a specified directory and returns them as a named list of data frames,
#' each representing a separate species dataset.
#'
#' This function assumes that each file in the directory is a pre-cleaned RDS object with a consistent structure.
#'
#' @param datadir A character string specifying the directory containing the `.rds` files.
#'   Default is `"/data/model-data/"`.
#'
#' @return A named list of species datasets (data frames).
#'
#' @examples
#' \dontrun{
#' species_data <- import_species_data(datadir = "myfolder/cleaned-data/")
#' }
#'
#' @export
import_species_data <- function(datadir = "/data/model-data/") {
  file_list <- list.files(path = datadir, pattern = "\\.rds$", full.names = TRUE)
  names(file_list) <- tools::file_path_sans_ext(basename(file_list))
  out_list <- lapply(file_list, readRDS)
  return(out_list)
}

#' Create pooled ataset
#'
#' Combines a list of individual species observation datasets into a single pooled dataset.
#'
#' Only columns specified in `keep_cols` will be retained. This function is useful for preparing
#' pooled detection histories for plotting or combined model input.
#'
#' @param species_datasets_list A named list of species observation data frames.
#' @param keep_cols A character vector of column names to retain. Defaults to `c("x1km_grid", "date", "year", "y1")`.
#'
#' @importFrom dplyr %>% select any_of arrange bind_rows
#'
#' @return A single data frame combining all input datasets with selected columns only.
#'
#' @examples
#' \dontrun{
#' pooled <- create_pooled(species_datasets_list)
#' }
#'
#' @export
create_pooled <- function(species_datasets_list,
                          keep_cols = c("x1km_grid", "date", "year", "y1")) {
  pooled <- lapply(species_datasets_list, function(df) {
    df %>%
      select(any_of(keep_cols)) %>%
      select(c("x1km_grid", "date", "year", "y1")) %>%
      arrange("x1km_grid", "year", "date")
  })
  pooled <- do.call(bind_rows, pooled)
  return(pooled)
}

#' Summarise detection histories
#'
#' Computes a set of basic summary statistics describing detection rates for each species dataset.
#'
#' This includes naive occupancy (proportion of sites ever visited that had ≥1 detection),
#' overall reporting rate (detections per visit), and reporting rate conditional on historical occupancy.
#'
#' @param species_datasets_list A named list of species datasets. Each data frame must contain `x1km_grid`, `y1`, and `siteID`.
#'
#' @importFrom dplyr %>% group_by summarise filter pull
#'
#' @return A data frame with one row per dataset and columns: `dataset`, `naiveOcc`, `reportingRate`, and `reportingRate_z1`.
#'
#' @examples
#' \dontrun{
#' summary_df <- summarise_data(species_datasets_list)
#' }
#'
#' @keywords internal
summarise_data <- function(species_datasets_list) {
  out <- lapply(names(species_datasets_list), function(nm) {
    df <- species_datasets_list[[nm]]
    z1_sites <- df %>%
      group_by(siteID) %>%
      summarise(ever_detected = any(y1 == 1, na.rm = TRUE), .groups = "drop") %>%
      filter(ever_detected) %>%
      pull(siteID)

    data.frame(
      dataset = nm,
      naiveOcc = mean(df$y1 > 0, na.rm = TRUE),
      reportingRate = mean(df$y1, na.rm = TRUE),
      reportingRate_z1 = mean(df$y1[df$siteID %in% z1_sites], na.rm = TRUE)
    )
  })
  return(dplyr::bind_rows(out))
}


