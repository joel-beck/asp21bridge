# Packages
pacman::p_load(
  dplyr,
  tidyr,
  purrr
)


### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### Creates data frame with results for a single matrix                     ####

collect_results <- function(samples, include_plot = FALSE) {
  full_data <- samples %>%
    as.data.frame() %>%
    tidyr::pivot_longer(
      cols = dplyr::everything(), names_to = "Parameter", values_to = "values"
    ) %>%
    dplyr::nest_by(Parameter, .key = "Data") %>%
    dplyr::mutate(
      "5% Quantile" = stats::quantile(Data$values, probs = 0.05),
      "Posterior Mean" = mean(Data$values),
      "Posterior Median" = stats::median(Data$values),
      "95% Quantile" = stats::quantile(Data$values, probs = 0.95),
      "Standard Deviation" = stats::sd(Data$values)
    )

  full_data <- full_data %>%
    dplyr::ungroup() %>%
    dplyr::mutate(Data = purrr::map2(
      .x = Data, .y = Parameter,
      .f = ~ dplyr::rename_at(.tbl = .x, .vars = 1, .funs = function(z) paste0(.y, ""))
    ))

  if (include_plot) {
    full_data <- full_data %>%
      dplyr::rowwise() %>%
      dplyr::mutate(Plot = list(diagnostic_plots(
        samples = Data, lag_max = 30, latex = TRUE
      )))
  }

  output <- full_data %>%
    dplyr::ungroup() %>%
    dplyr::select(-Data)

  # drop Parameter column, if Input is numeric vector
  if (is.vector(samples, mode = "numeric")) {
    output <- output %>%
      dplyr::select(-Parameter)
  }

  # if no column names of input matrix are provided
  if (is.matrix(samples) && is.null(colnames(samples))) {
    output$Parameter <- paste("col", 1:nrow(output), sep = "_")
  }

  return(output)
}



### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### Combines results for all matrices in the input list into one data frame ####

summary_complete <- function(samples, include_plot = FALSE) {
  if (class(samples) == "lslm") {
    # if input is model object
    purrr::map_dfr(
      .x = samples$mcmc_ridge[[1]], .f = collect_results,
      include_plot = include_plot
    )
  } else if (is.matrix(samples) || is.vector(samples, mode = "numeric")) {
    # if input is a matrix or a numeric vector
    collect_results(samples = samples, include_plot = include_plot)
  } else if (is.list(samples)) {
    # if input is a list but not a model object
    purrr::map_dfr(.x = samples, .f = collect_results, include_plot = include_plot)
  } else {
    stop(paste(
      "Input 'samples' must be a numeric vector,",
      "matrix or list of vectors and matrices!"
    ))
  }
}
