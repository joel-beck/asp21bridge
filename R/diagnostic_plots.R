
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param samples PARAM_DESCRIPTION
#' @param robust PARAM_DESCRIPTION, Default: FALSE
#' @param latex PARAM_DESCRIPTION, Default: FALSE
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname time_plot
#' @export

time_plot <- function(samples, log = FALSE, robust = FALSE, latex = FALSE) {

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

  name <- colnames(data)
  data$values <- data[[1]]
  data$time <- seq_len(nrow(data))


  # Error Messages
  if (name == "samples" && latex) {
    stop("Set 'latex = TRUE' only if column names are provided.")
  }


  plot <- ggplot2::ggplot(data = data, mapping = ggplot2::aes(x = time, y = values)) +
    ggplot2::geom_line() +
    ggplot2::labs(x = "Iterations", y = NULL, title = "Time Plot") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5),
      axis.title = ggplot2::element_text(size = 10)
    )

  if (log) {
    plot <- plot +
      ggplot2::scale_y_log10()
  }


  # Cut off Tails of Distribution
  if (robust) {
    plot <- plot +
      ggplot2::coord_cartesian(ylim = stats::quantile(data$values, probs = c(0.01, 0.99)))
  }


  # Format Plot Label in LaTeX
  if (latex) {
    name <- stringr::str_glue("$\\{name}$")
    plot <- plot +
      ggplot2::labs(title = latex2exp::TeX(paste("Time Plot for", name)))
  }


  return(plot)
}





#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param samples PARAM_DESCRIPTION
#' @param robust PARAM_DESCRIPTION, Default: FALSE
#' @param latex PARAM_DESCRIPTION, Default: FALSE
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname density_plot
#' @export

density_plot <- function(samples, log = FALSE, robust = FALSE, latex = FALSE) {

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

  name <- colnames(data)
  data$values <- data[[1]]


  # Error Messages
  if (name == "samples" && latex) {
    stop("Set 'latex = TRUE' only if column names are provided.")
  }


  plot <- ggplot2::ggplot(data = data, mapping = ggplot2::aes(x = values)) +
    ggplot2::geom_density() +
    ggplot2::geom_point(mapping = ggplot2::aes(y = 0), size = 0.5, shape = 4, color = "grey50") +
    ggplot2::labs(x = "Values", y = NULL, title = "Density Estimate") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5),
      axis.title = ggplot2::element_text(size = 10)
    )


  if (log) {
    plot <- plot +
      ggplot2::scale_x_log10()
  }

  # Cut off Tails of Distribution
  if (robust) {
    plot <- plot +
      ggplot2::coord_cartesian(xlim = stats::quantile(data$values, probs = c(0.01, 0.99)))
  }

  # Format Plot Label in LaTeX
  if (latex) {
    name <- stringr::str_glue("$\\{name}$")
    plot <- plot +
      ggplot2::labs(title = latex2exp::TeX(paste("Density Estimate for", name)))
  }

  return(plot)
}







#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param samples PARAM_DESCRIPTION
#' @param lag_max PARAM_DESCRIPTION, Default: 30
#' @param latex PARAM_DESCRIPTION, Default: FALSE
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname acl_plot
#' @export

acl_plot <- function(samples, lag_max = 30, latex = FALSE) {

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

  data$time <- 0:(length(data$acl) - 1)


  # Error Messages
  if (name == "samples" && latex) {
    stop("Set 'latex = TRUE' only if column names are provided.")
  }


  plot <- ggplot2::ggplot(data = data, mapping = ggplot2::aes(x = time, y = acl)) +
    ggplot2::geom_point(size = 1.5) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    ggplot2::geom_segment(mapping = ggplot2::aes(xend = time, yend = 0)) +
    ggplot2::coord_cartesian(xlim = c(1, min(length(data$time), lag_max))) +
    ggplot2::labs(x = "Lag", y = NULL, title = "Autocorrelation") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5),
      axis.title = ggplot2::element_text(size = 10)
    )


  # Format Plot Label in LaTeX
  if (latex) {
    name <- stringr::str_glue("$\\{name}$")
    plot <- plot +
      ggplot2::labs(title = latex2exp::TeX(paste("Autocorrelation for", name)))
  }

  return(plot)
}






#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param samples PARAM_DESCRIPTION
#' @param lag_max PARAM_DESCRIPTION, Default: 30
#' @param robust PARAM_DESCRIPTION, Default: FALSE
#' @param latex PARAM_DESCRIPTION, Default: FALSE
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname diagnostic_plots
#' @export

diagnostic_plots <- function(samples, lag_max = 30, log = FALSE,
                             robust = FALSE, latex = FALSE) {
  time_plot <- time_plot(samples = samples, log = log, robust = robust, latex = latex)
  density_plot <- density_plot(samples = samples, log = log, robust = robust)
  acl_plot <- acl_plot(samples = samples, lag_max = lag_max)

  time_plot / (density_plot | acl_plot)
}
