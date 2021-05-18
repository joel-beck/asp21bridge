burn_element <- function(samples, num_burn) {
  if (is.matrix(samples)) {
    samples[(num_burn + 1):nrow(samples), , drop = FALSE]
  } else if (is.vector(samples, mode = "numeric")) {
    samples[(num_burn + 1):length(samples)]
  } else {
    stop(paste(
      "Input 'samples' must be a numeric vector,",
      "matrix or list of vectors and matrices!"
    ))
  }
}


#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param samples PARAM_DESCRIPTION
#' @param num_burn PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname burnin
#' @export

burnin <- function(samples, num_burn) {
  if (is.matrix(samples) || is.vector(samples, mode = "numeric")) {
    burn_element(samples = samples, num_burn = num_burn)
  } else if (is.list(samples)) {
    lapply(X = samples, burn_element, num_burn = num_burn)
  } else {
    stop(paste(
      "Input 'samples' must be a numeric vector,",
      "matrix or list of vectors and matrices!"
    ))
  }
}
