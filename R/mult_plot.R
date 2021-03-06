#   ____________________________________________________________________________
#   Helper Functions                                                        ####

#' @importFrom rlang .data

mult_trace <- function(data, log, free_scale, latex) {
  free_scale <- ifelse(test = free_scale, yes = "free_y", no = "fixed")

  trace_p <- data %>%
    ggplot2::ggplot(mapping = ggplot2::aes(
      x = .data$time, y = .data$vals, color = .data$Parameter
    )) +
    ggplot2::geom_line() +
    ggplot2::labs(x = "Iterations", y = NULL, title = "Trace Plots") +
    ggplot2::facet_wrap(facets = ggplot2::vars(.data$Parameter), scales = free_scale) +
    ggplot2::guides(color = "none") +
    ggplot2::theme_minimal(base_size = 9) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

  if (log) {
    trace_p <- trace_p +
      ggplot2::scale_y_log10()
  }

  if (latex) {
    trace_p <- trace_p +
      ggplot2::facet_wrap(
        facets = ggplot2::vars(.data$Parameter), scales = free_scale,
        labeller = ggplot2::label_parsed
      )
  }
  return(trace_p)
}


#' @importFrom rlang .data

mult_density <- function(data, log, latex) {
  density_p <- data %>%
    ggplot2::ggplot(mapping = ggplot2::aes(x = .data$vals, fill = .data$Parameter)) +
    ggridges::geom_density_ridges(
      mapping = ggplot2::aes(y = .data$Parameter %>% forcats::fct_rev()), alpha = 0.5
    ) +
    ggplot2::labs(x = "Values", y = NULL, title = "Density Estimates") +
    ggplot2::guides(fill = "none") +
    ggplot2::theme_minimal(base_size = 9) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5),
      axis.title = ggplot2::element_text(size = 10)
    )

  if (log) {
    density_p <- density_p +
      ggplot2::scale_x_log10()
  }

  if (latex) {
    density_p <- density_p +
      ggplot2::scale_y_discrete(
        labels = parse(text = levels(data$Parameter %>% forcats::fct_rev()))
      )
  }

  return(density_p)
}



#   ____________________________________________________________________________
#   Main Function                                                           ####

#' @title Plot multiple Markov Chains together
#'
#' @description The `mult_plot()` function allows to combine trace and / or
#'              density plots of multiple Markov Chains into a single graphical
#'              illustration. \cr
#'              Trace plots are displayed separately in facets, density plots
#'              are shown in a ridgeline plot.
#'
#' @param samples Model object of the class 'lmls', list or matrix
#'                containing samples of the MCMC sampler.
#'
#' @param type One of the values "trace", "density" or "both". \cr
#'             Controls whether only trace plots or only density plots should
#'             be included for all chains.
#'             The option "both" stacks trace plots on top of the density plots.
#'
#' @param log Logical. If TRUE, the y - axis in trace plots and the x - axis in
#'            density plots are transformed to the logarithmic scale. \cr
#'            Default: FALSE
#'
#' @param free_scale Logical. If TRUE, the y - axis scale in trace plots is
#'                   chosen differently for all facets. \cr
#'                   Default: FALSE
#'
#' @param latex Logical. If TRUE, facets in trace plots and ridgelines in density
#'              plots can be labeled by rendered mathematical symbols such as
#'              greek letters in combination with subscripts and superscripts.
#'              This option requires the columns of the input matrices to be
#'              labeled accordingly (as in the output of the `mcmc_ridge()`
#'              function). \cr
#'              Default: FALSE
#'
#' @returns Plot object of the class "ggplot".
#'
#' @examples
#' fit <- lmls(
#'   location = y ~ x1 + x2 + z1 + z2, scale = ~ z1 + z2,
#'   data = toy_data, light = FALSE
#' ) %>%
#'   mcmc_ridge(num_sim = 1000)
#'
#' # list of 4 matrices with 1000 rows each
#' samples <- fit$mcmc_ridge$sampling_matrices
#'
#' # plots all output list elements together for quick overview
#' mult_plot(samples, type = "both", free_scale = TRUE, latex = TRUE)
#'
#' # trace plot of location coefficients
#' mult_plot(samples$location, type = "trace", free_scale = TRUE, latex = TRUE)
#'
#' # log argument often useful for strictly positive variance parameters
#' mult_plot(
#'   list(samples$location_prior, samples$scale_prior),
#'   type = "both", log = TRUE, latex = TRUE
#' )
#'
#' @export
#'
#' @importFrom rlang .data
#' @import patchwork

mult_plot <- function(samples, type = c("trace", "density", "both"),
                      log = FALSE, free_scale = FALSE, latex = FALSE) {

  # validate input ----------------------------------------------------------

  type <- rlang::arg_match(type)

  if (!(is.logical(free_scale) && is.logical(log) && is.logical(latex))) {
    stop(paste(
      "Inputs 'log', 'free_scale' and 'latex'",
      "must be either 'TRUE' or 'FALSE'!"
    ))
  }

  if (is.matrix(samples)) {
    data <- as.data.frame(samples)
  } else if (class(samples) == "lmls") {
    # if input is model object
    data <- purrr::map_dfc(.x = samples$mcmc_ridge[[1]], .f = as.data.frame)
  } else if (is.list(samples)) {
    # if input is a list but not a model object
    if (any(names(samples) == "sampling_matrices")) {
      # if input is mcmc_ridge() output in list format with two entries
      data <- purrr::map_dfc(.x = samples[["sampling_matrices"]], .f = as.data.frame)
    } else {
      # if input is list with just the sampling matrices
      data <- purrr::map_dfc(.x = samples, .f = as.data.frame)
    }
  } else {
    stop(paste(
      "Input 'samples' must be a model object of the class 'lmls',",
      "a numeric matrix or a list of vectors and matrices!"
    ))
  }

  # common data for all output types ----------------------------------------

  if (latex) {
    # replace underscores in column names with brackets
    colnames(data) <- stringr::str_replace(
      string = colnames(data), pattern = "_(.*)", replacement = "[\\1]"
    )
  }

  data <- data %>%
    dplyr::mutate(time = seq_len(nrow(data))) %>%
    tidyr::pivot_longer(
      cols = -.data$time, names_to = "Parameter", values_to = "vals"
    ) %>%
    dplyr::mutate(Parameter = factor(.data$Parameter) %>% forcats::fct_inorder())

  # trace plot ---------------------------------------------------------------

  trace_p <- mult_trace(
    data = data, log = log, free_scale = free_scale, latex = latex
    )

  # density plot ------------------------------------------------------------

  density_p <- mult_density(data = data, log = log, latex = latex)

  # output depending on type ------------------------------------------------

  if (type == "both") {
    plot <- trace_p / density_p
  } else if (type == "trace") {
    plot <- trace_p
  } else if (type == "density") {
    plot <- density_p
  }

  # avoid console message of gg_ridges
  suppressMessages(print(plot))
}
