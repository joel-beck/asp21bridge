
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param samples PARAM_DESCRIPTION
#' @param type PARAM_DESCRIPTION, Default: c("time", "density", "both")
#' @param robust PARAM_DESCRIPTION, Default: FALSE
#' @param free_scale PARAM_DESCRIPTION, Default: FALSE
#' @param latex PARAM_DESCRIPTION, Default: FALSE
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # EXAMPLE1
#' }
#' }
#' @rdname mult_plot
#' @export

mult_plot <- function(samples, type = c("time", "density", "both"),
                      log = FALSE, robust = FALSE, free_scale = FALSE,
                      latex = FALSE) {
  type <- rlang::arg_match(type)

  # Error Messages
  if (!is.list(samples) && !is.matrix(samples)) {
    stop("Input 'samples' must be a numeric matrix or a list of matrices!")
  }

  if (!(is.logical(robust) && is.logical(free_scale) && is.logical(log) &&
    is.logical(latex))) {
    stop(paste(
      "Inputs 'log', 'robust', 'free_scale' and 'latex'",
      "must be either 'TRUE' or 'FALSE'!"
    ))
  }

  # rows_matrix <- purrr::map_dbl(.x = samples, .f = nrow)
  # if (any(rows_matrix != rows_matrix[1])) {
  #   stop("The number of simulations must be equal for all parameters!")
  # }

  # Plotting Data for both types
  if (is.list(samples)) {
    data <- purrr::map_dfc(.x = samples, .f = as.data.frame)
  }

  if (is.matrix(samples)) {
    data <- as.data.frame(samples)
  }


  if (latex) {
    # replace underscores in column names with brackets
    colnames(data) <- stringr::str_replace(
      string = colnames(data), pattern = "_(.*)", replacement = "[\\1]"
    )
  }

  data <- data %>%
    dplyr::mutate(time = seq_len(nrow(data))) %>%
    tidyr::pivot_longer(
      cols = -time, names_to = "Parameter", values_to = "vals"
    ) %>%
    dplyr::mutate(Parameter = factor(Parameter) %>% forcats::fct_inorder())



  # Time Plot
  free_scale <- ifelse(test = free_scale, yes = "free_y", no = "fixed")

  time_p <- data %>%
    ggplot2::ggplot(mapping = ggplot2::aes(x = time, y = vals, color = Parameter)) +
    ggplot2::geom_line() +
    ggplot2::labs(x = "Iterations", y = NULL, title = "Time Plots") +
    ggplot2::facet_wrap(facets = ggplot2::vars(Parameter), scales = free_scale) +
    ggplot2::guides(color = FALSE) +
    ggplot2::theme_minimal() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

  if (log) {
    time_p <- time_p +
      ggplot2::scale_y_log10()
  }

  if (latex) {
    time_p <- time_p +
      ggplot2::facet_wrap(
        facets = ggplot2::vars(Parameter), scales = free_scale,
        labeller = ggplot2::label_parsed
      )
  }



  # Density Plot
  density_p <- data %>%
    ggplot2::ggplot(mapping = ggplot2::aes(x = vals, fill = Parameter)) +
    ggridges::geom_density_ridges(
      mapping = ggplot2::aes(y = Parameter %>% forcats::fct_rev()), alpha = 0.5
    ) +
    ggplot2::labs(x = "Values", y = NULL, title = "Density Estimates") +
    ggplot2::guides(fill = FALSE) +
    ggplot2::theme_minimal() +
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
      ggplot2::scale_y_discrete(labels = parse(text = levels(data$Parameter)))
  }


  if (robust) {
    bounds <- data %>%
      dplyr::summarise(bounds = stats::quantile(vals, probs = c(0.01, 0.99)))
    density_p <- density_p +
      ggplot2::coord_cartesian(xlim = bounds$bounds)
  }



  # return plot object depending on input "type"
  if (type == "both") {
    plot <- time_p / density_p
  }

  if (type == "time") {
    plot <- time_p
  }

  if (type == "density") {
    plot <- density_p
  }

  return(plot)
}
