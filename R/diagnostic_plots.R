#' @title Plot values of a single Markov Chain over time
#'
#' @description The `trace_plot()` function creates a line graph of Markov Chain
#'              samples, where the sampled values are mapped to the y - axis
#'              and the iteration / index is placed on the x - axis in
#'              chronological order. \cr
#'              This type of plot is typically used to evaluate the convergence
#'              of the chain, possibly towards a posterior distribution.
#'
#' @param samples Numeric vector, matrix with a single column or list
#'                containing such a vector or matrix as the only element.
#'
#' @param log Logical. If TRUE, the y - axis is transformed to the
#'            logarithmic scale. \cr
#'            Default: FALSE
#'
#' @param robust Logical. If TRUE, the first and last percentile of the samples'
#'               distribution are omitted to avoid a strong influence of
#'               outliers on the y - axis scale. \cr
#'               Default: FALSE
#'
#' @param latex Logical. If TRUE, mathematical symbols such as greek letters
#'              in the plot title with subscripts and superscripts are properly
#'              rendered.
#'              This option requires the column of the input matrix to be
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
#' # extract first chain of the location parameter
#' diagnostic_plots(samples$location[, 1, drop = FALSE], latex = TRUE)
#'
#' # log argument often useful for strictly positive variance parameters
#' trace_plot(samples$scale_prior, log = TRUE, latex = TRUE)
#'
#' @export
#'
#' @importFrom rlang .data

trace_plot <- function(samples, log = FALSE, robust = FALSE, latex = FALSE) {

  # validate input ----------------------------------------------------------

  # check for data frame included to display correct column names in
  # summary_complete(., include_plot = TRUE)
  if (is.data.frame(samples)) {
    data <- samples
  } else if (is.list(samples) && length(samples) == 1) {
    data <- as.data.frame(samples[[1]])
  } else if ((is.matrix(samples) && ncol(samples) == 1) ||
    is.vector(samples, mode = "numeric")) {
    data <- as.data.frame(samples)
  } else {
    stop(paste(
      "Input 'samples' must be a numeric vector, a matrix with one column",
      "or a list with a single element of these types!"
    ))
  }

  # include check !is.data.frame() for latex labels in summary_complete()
  if ((!is.list(samples) && is.null(colnames(samples)) ||
       is.list(samples) && !is.data.frame(samples) && is.null(colnames(samples[[1]])))
      && latex) {
    stop("Set 'latex = TRUE' only if column names are provided.")
  }

  # basic plot -------------------------------------------------------------

  name <- colnames(data)
  data$values <- data[[1]]
  data$time <- seq_len(nrow(data))

  plot <- ggplot2::ggplot(
    data = data, mapping = ggplot2::aes(x = .data$time, y = .data$values)
  ) +
    ggplot2::geom_line() +
    ggplot2::labs(x = "Iterations", y = NULL, title = "Trace Plot") +
    ggplot2::theme_minimal(base_size = 9) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5),
      axis.title = ggplot2::element_text(size = 10)
    )

  # plot options ------------------------------------------------------------

  if (log) {
    plot <- plot +
      ggplot2::scale_y_log10()
  }

  if (robust) {
    plot <- plot +
      ggplot2::coord_cartesian(
        ylim = stats::quantile(data$values, probs = c(0.01, 0.99))
      )
  }

  if (latex) {
    name <- stringr::str_glue("$\\{name}$")
    plot <- plot +
      ggplot2::labs(title = latex2exp::TeX(paste("Trace Plot for", name)))
  }

  return(plot)
}





#' @title Plot kernel density estimate for values of a single Markov Chain
#'
#' @description The `density_plot()` function creates a kernel density plot of
#'              Markov Chain samples with the sampled values mapped to the
#'              x - axis.
#'              The individual data points are displayed by gray points on the
#'              x - axis to get a more detailed overview of the distribution.
#'
#' @param samples Numeric vector, matrix with a single column or list
#'                containing such a vector or matrix as the only element.
#'
#' @param log Logical. If TRUE, the x - axis is transformed to the
#'            logarithmic scale. \cr
#'            Default: FALSE
#'
#' @param robust Logical. If TRUE, the first and last percentile of the samples'
#'               distribution are omitted to avoid a strong influence of
#'               outliers on the x - axis scale. \cr
#'               Default: FALSE
#'
#' @param latex Logical. If TRUE, mathematical symbols such as greek letters
#'              in the plot title with subscripts and superscripts are properly
#'              rendered.
#'              This option requires the column of the input matrix to be
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
#' # extract first chain of the scale parameter
#' density_plot(samples$scale[, 1], robust = TRUE)
#'
#' # log argument often useful for strictly positive variance parameters
#' density_plot(samples$scale_prior, log = TRUE, latex = TRUE)
#'
#' @export
#'
#' @importFrom rlang .data

density_plot <- function(samples, log = FALSE, robust = FALSE, latex = FALSE) {

  # validate input ----------------------------------------------------------

  if (is.list(samples) && length(samples) == 1) {
    data <- as.data.frame(samples[[1]])
  } else if ((is.matrix(samples) && ncol(samples) == 1) ||
    is.vector(samples, mode = "numeric")) {
    data <- as.data.frame(samples)
  } else {
    stop(paste(
      "Input 'samples' must be a numeric vector, a matrix with one column",
      "or a list with a single element of these types!"
    ))
  }

  if ((!is.list(samples) && is.null(colnames(samples)) ||
       is.list(samples) && is.null(colnames(samples[[1]])))
      && latex) {
    stop("Set 'latex = TRUE' only if column names are provided.")
  }

  # basic plot --------------------------------------------------------------

  name <- colnames(data)
  data$values <- data[[1]]

  plot <- ggplot2::ggplot(
    data = data, mapping = ggplot2::aes(x = .data$values)
  ) +
    ggplot2::geom_density() +
    ggplot2::geom_point(
      mapping = ggplot2::aes(y = 0), size = 0.5, shape = 4,
      color = "grey50"
    ) +
    ggplot2::labs(x = "Values", y = NULL, title = "Density Estimate") +
    ggplot2::theme_minimal(base_size = 9) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5),
      axis.title = ggplot2::element_text(size = 10)
    )

  # plot options ------------------------------------------------------------

  if (log) {
    plot <- plot +
      ggplot2::scale_x_log10()
  }

  if (robust) {
    plot <- plot +
      ggplot2::coord_cartesian(
        xlim = stats::quantile(data$values, probs = c(0.01, 0.99))
      )
  }

  if (latex) {
    name <- stringr::str_glue("$\\{name}$")
    plot <- plot +
      ggplot2::labs(title = latex2exp::TeX(paste("Density Estimate for", name)))
  }

  return(plot)
}







#' @title Plot autocorrelations of a single Markov Chain
#'
#' @description The `acf_plot()` function creates a graphical display of
#'              autocorrelations for various lags of a single Markov Chain.
#'              The lags are mapped to the x - axis in increasing order with the
#'              corresponding autocorrelation values mapped to the y - axis. \cr
#'              This type of plot can be used to analyze the dependence of the
#'              collected samples i.e. from a posterior distribution.
#'
#' @param samples Numeric vector, matrix with a single column or list
#'                containing such a vector or matrix as the only element.
#'
#' @param lag_max Positive integer representing the maximum number of lags
#'                that are shown on the x - axis. \cr
#'                Default: 30
#'
#' @param latex Logical. If TRUE, mathematical symbols such as greek letters
#'              in the plot title with subscripts and superscripts are properly
#'              rendered.
#'              This option requires the column of the input matrix to be
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
#' # uses default lag_max = 30
#' acf_plot(samples$scale_prior, latex = TRUE)
#'
#' # value of lag_max should be adapted depending on the correlation structure
#' acf_plot(samples$scale[, 1, drop = FALSE], lag_max = 100, latex = TRUE)
#'
#' @export
#'
#' @importFrom rlang .data

acf_plot <- function(samples, lag_max = 30, latex = FALSE) {

  # validate input ----------------------------------------------------------

  if (is.list(samples) && length(samples) == 1) {
    name <- colnames(as.data.frame(samples[[1]]))
    data <- data.frame(
      acl = drop(stats::acf(x = samples[[1]], lag.max = lag_max, plot = FALSE)$acf)
    )
  } else if ((is.matrix(samples) && ncol(samples) == 1) ||
    is.vector(samples, mode = "numeric")) {
    name <- colnames(as.data.frame(samples))
    data <- data.frame(
      acl = drop(stats::acf(x = samples, lag.max = lag_max, plot = FALSE)$acf)
    )
  } else {
    stop(paste(
      "Input 'samples' must be a numeric vector, a matrix with one column",
      "or a list with a single element of these types!"
    ))
  }

  if ((!is.list(samples) && is.null(colnames(samples)) ||
       is.list(samples) && is.null(colnames(samples[[1]])))
      && latex) {
    stop("Set 'latex = TRUE' only if column names are provided.")
  }

  # basic plot --------------------------------------------------------------

  data$time <- 0:(length(data$acl) - 1)

  plot <- ggplot2::ggplot(
    data = data, mapping = ggplot2::aes(x = .data$time, y = .data$acl)
  ) +
    ggplot2::geom_point(size = 1.5) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    ggplot2::geom_segment(mapping = ggplot2::aes(xend = .data$time, yend = 0)) +
    ggplot2::coord_cartesian(xlim = c(1, min(length(data$time), lag_max))) +
    ggplot2::labs(x = "Lag", y = NULL, title = "Autocorrelation") +
    ggplot2::theme_minimal(base_size = 9) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5),
      axis.title = ggplot2::element_text(size = 10)
    )

  # plot options ------------------------------------------------------------

  if (latex) {
    name <- stringr::str_glue("$\\{name}$")
    plot <- plot +
      ggplot2::labs(title = latex2exp::TeX(paste("Autocorrelation for", name)))
  }

  return(plot)
}






#' @title Combine trace, density and autocorrelation plots for a single
#'        Markov Chain into one graphical display
#'
#' @description The `diagnostic_plots()` function stacks a trace plot on top
#'              of side-by-side density and autocorrelation plots for a single
#'              Markov Chain. \cr
#'              This functionality is particularly convenient to analyze both
#'              the convergence as well as the dependence of the collected
#'              samples.
#'
#' @param samples Numeric vector, matrix with a single column or list
#'                containing such a vector or matrix as the only element.
#'
#' @param lag_max Positive integer representing the maximum number of lags
#'                that are shown on the x - axis of the autocorrelation plot. \cr
#'                Default: 30
#'
#' @param log Logical. If TRUE, the y - axis in the trace plot and the x - axis
#'            in the density plot are transformed to the logarithmic scale. \cr
#'            Default: FALSE
#'
#' @param robust Logical. If TRUE, the first and last percentile of the samples'
#'               distribution are omitted from the trace plot as well as the
#'               density plot to avoid a strong influence of outliers. \cr
#'               Default: FALSE
#'
#' @param latex Logical. If TRUE, mathematical symbols such as greek letters
#'              in the plot title with subscripts and superscripts are properly
#'              rendered.
#'              This option requires the column of the input matrix to be
#'              labeled accordingly (as in the output of the `mcmc_ridge()`
#'              function). \cr
#'              Default: FALSE
#'
#' @returns Plot object of the classes "patchwork" and "ggplot".
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
#' # extract first chain of the location parameter
#' diagnostic_plots(samples$location[, 1, drop = FALSE], latex = TRUE)
#'
#' # log argument often useful for strictly positive variance parameters
#' diagnostic_plots(samples$scale_prior, log = TRUE, robust = TRUE, latex = TRUE)
#'
#' @export
#'
#' @import patchwork

diagnostic_plots <- function(samples, lag_max = 30, log = FALSE,
                             robust = FALSE, latex = FALSE) {
  trace_plot <- trace_plot(samples = samples, log = log, robust = robust, latex = latex)
  density_plot <- density_plot(samples = samples, log = log, robust = robust)
  acf_plot <- acf_plot(samples = samples, lag_max = lag_max)

  trace_plot / (density_plot | acf_plot)
}
