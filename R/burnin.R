burn_element <- function(samples, num_burn) {
  # check for validity of 'num_burn'
  if ((as.integer(num_burn) != num_burn) || num_burn < 0) {
    stop("'num_burn' must be a nonnegative integer!")
  }

  if (is.matrix(samples)) {
    if (num_burn >= nrow(samples)) {
      stop("'num_burn' must be less than the number of samples!")
    }
    samples[(num_burn + 1):nrow(samples), , drop = FALSE]
  } else if (is.vector(samples, mode = "numeric")) {
    if (num_burn >= length(samples)) {
      stop("'num_burn' must be less than the number of samples!")
    }
    samples[(num_burn + 1):length(samples)]
  } else {
    stop(paste(
      "Input 'samples' must be a numeric vector,",
      "matrix or list of vectors and matrices!"
    ))
  }
}


#' @title Remove first `num_burn` samples of a Markov Chain
#'
#' @description
#' `burnin()` takes a numeric vector, matrix or list of matrices as input and
#' returns a modified copy of the input object with a removed first part of
#' the original entries.
#' The function is primarily used in the context of Markov-Chain Monte-Carlo
#' sampling to obtain samples of the posterior distribution after a converging
#' phase.
#'
#' @param samples Numeric vector, matrix or list of matrices containing the
#'                samples of the MCMC sampler.
#'
#' @param num_burn Nonnegative integer which controls the number of removed
#'                 samples, e.g. `num_burn = 10` removes the first 10 samples.
#'
#' @returns
#' Same data type as the provided input, i.e. a numeric vector, matrix or list
#' of matrices with removed elements at the beginning.
#'
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
