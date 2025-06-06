

#' Plot temporal bias in species datasets
#'
#' Generates a line plot showing the number of sites recorded per year for each dataset.
#'
#' @param species_datasets_list A named list of data frames, each containing columns `year` and `x1km_grid`.
#'
#' @return A ggplot object showing temporal variation in recording effort.
#'
#' @importFrom dplyr group_by summarise mutate bind_rows n_distinct
#' @importFrom ggplot2 ggplot aes geom_line geom_point theme_bw scale_x_continuous labs scale_color_brewer
#' @export
plot_temporal_bias <- function(species_datasets_list) {

  # Compute sites‐per‐year for each dataset and add a dataset column
  plotData <- lapply(names(species_datasets_list), function(nm) {
    species_datasets_list[[nm]] %>%
      group_by(year) %>%
      summarise(nsites = n_distinct(x1km_grid), .groups = "drop") %>%
      mutate(dataset = nm)
  }) %>%
    bind_rows()

  years <- c(min(plotData$year):max(plotData$year))

  # Plot
  ggplot(plotData, aes(x = year, y = nsites, color = dataset)) +
    geom_line(size = 1) +
    geom_point(size = 2) +
    theme_bw() +
    scale_x_continuous(limits = c(min(years), max(years)),
                       breaks = c(min(years):max(years)),
                       labels=years) +
    labs(x = "Year", y = "Number of sites",
         color = "Dataset",
         title = "Annual recording effort") +
    scale_color_brewer(palette = "Paired")
}

#' Plot annual occupancy trend from model output
#'
#' Plots predicted mean annual occupancy (psi.fs) with credible intervals.
#'
#' @param samples A coda MCMC object (e.g. output from runMCMC with samplesAsCodaMCMC = TRUE).
#' @param years A vector of years corresponding to the psi.fs estimates.
#'
#' @return A ggplot object visualising posterior median and credible intervals.
#'
#' @importFrom ggplot2 ggplot aes geom_line geom_ribbon scale_x_continuous scale_y_continuous ylab xlab theme_light ggtitle
#' @export
plot_occupancy_trend <- function(samples, years) {

  # 1) get the trend
  # first coerce the samples into a format with year, mean and CI
  qntiles <- summary(samples)$quantiles
  #mulam <-  as.data.frame(output[grepl("mu.lambda", dimnames(output)[[1]]),])
  psi.fs <- as.data.frame(qntiles[grepl("psi.fs", dimnames(qntiles)[[1]]),])

  # 2) then plot
  trend_plot <- psi.fs %>%
    ggplot(aes(x=1:length(years))) +
    geom_line(aes(y=`50%`), color = "steelblue4") +
    geom_ribbon(aes(ymin=`25%`, ymax=`75%`), fill = "steelblue", alpha=0.5) +
    geom_ribbon(aes(ymin=`2.5%`, ymax=`97.5%`), fill = "lightblue4", alpha=0.2) +
    scale_x_continuous(limits = c(1, length(years)), breaks = c(1:length(years)),
                       labels=years) +
    scale_y_continuous(limits = c(0,1), breaks = c(0.2, 0.4, 0.6, 0.8, 1.0)) +
    ylab("Occupancy") +
    xlab("Year") +
    theme_light() +
    ggtitle("Predicted annual occupancy")

  return(trend_plot)

}

#' Convert spatial covariate list to a wide-format data frame
#'
#' @param spatial_cov_list A named list of data frames with `x1km_grid`, `value`, `x`, and `y` columns.
#'
#' @return A single wide-format data frame, one row per grid cell.
#'
#' @importFrom dplyr bind_rows select mutate distinct left_join
#' @importFrom tidyr pivot_wider
#' @export
list_to_wide_df <- function(spatial_cov_list) {
  # stack into long first
  long_df <- bind_rows(
    lapply(names(spatial_cov_list), function(cov) {
      spatial_cov_list[[cov]] %>%
        select(x1km_grid, value) %>%
        mutate(covariate = cov)
    })
  )
  # then pivot to wide: one row per grid cell, cols = covariates
  wide_df <-  long_df %>%
    pivot_wider(names_from = covariate, values_from = value)

  # retrieve the x y coordinates from first list element
  coords <- spatial_cov_list[[1]] %>% select(x1km_grid, x, y) %>% distinct()
  map_layers <- coords %>% left_join(wide_df, by = "x1km_grid")

}

#' Predict occupancy (psi_hat) from map layers and coefficients
#'
#' @param map_layers A wide-format data frame of environmental covariates for each grid cell.
#' @param state_covariates List of model covariate names split by dimension.
#' @param pars Data frame with columns `Parameter` and `Mean`.
#'
#' @return A data frame with added columns: `lp` and `psi_hat` (logit-scale prediction and plogis).
#'
#'@importFrom stats plogis
#'
#' @export
predict_map <- function(map_layers, state_covariates, pars) {
  # intercept
  alpha.s <- pars$Mean[pars$Parameter == "alpha.s"]
  # names of the covariates in the model

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

  covs <- c(covariates_i, covariates_t, covariates_it)

  # slope estimates
  slope_names <- paste0("b_", covs)
  betas <- pars$Mean[match(slope_names, pars$Parameter)]
  # check no missing
  if (any(is.na(betas))) {
    stop("Some coefficient names not found in pars$Parameter: ",
         paste(slope_names[is.na(betas)], collapse = ", "))
  }
  # build the linear predictor = α + X %*% β
  Xmat <- as.matrix(map_layers[ , covs])
  lp <- alpha.s + drop(Xmat %*% betas)
  map_layers %>%
    mutate(lp = lp,
           psi_hat = plogis(lp)
    )
}

#' Plot occupancy prediction map
#'
#' Combines prediction and plotting in one function. Uses `ggplot2` and `viridis`.
#'
#' @param map_layers Raw spatial covariate layers.
#' @param state_covariates List of covariates (i, t, it).
#' @param pars Data frame with model output.
#' @param legend_name Expression for legend title.
#' @param col_palette Color palette (default = 'inferno').
#' @param rev_palette Reverse palette direction (logical).
#' @param plot_title Optional title.
#'
#' @return A ggplot object with spatial tiles of psi_hat.
#'
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_viridis_c coord_equal theme_minimal labs
#' @importFrom dplyr mutate
#' @export
plot_occupancy_map <- function(map_layers, state_covariates, pars,
                               legend_name = expression(hat(psi)),
                               col_palette = "inferno",
                               rev_palette = FALSE,
                               plot_title=NULL) {
  wide <- map_layers %>% list_to_wide_df()
  pred <- predict_map(wide, state_covariates, pars)
  map <- pred %>%
    ggplot(aes(x = x, y = y, fill = psi_hat)) +
    geom_tile() +
    scale_fill_viridis_c(option = col_palette,
                         direction = if (rev_palette) -1 else 1,
                         na.value = "transparent",
                         name = legend_name) +
    coord_equal() +
    theme_minimal() +
    labs(x="Easting", y="Northing")

  if(!is.null(plot_title)) {
    map <- map + labs(x="Easting", y="Northing", title = plot_title)
  }

  return(map)

}
