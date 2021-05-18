thin_element <- function(samples, freq) {
  # check for validity of 'freq'
  if ((as.integer(freq) != freq) || freq <= 0) {
    stop("'freq' must be a positive integer!")
  }

  if (is.matrix(samples)) {
    if (freq > nrow(samples)) {
      stop("'freq' must not be larger than the number of samples!")
    }
    samples[seq(from = 1, to = nrow(samples), by = freq), , drop = FALSE]
  } else if (is.vector(samples, mode = "numeric")) {
    if (freq > length(samples)) {
      stop("'freq' must not be larger than the number of samples!")
    }
    samples[seq(from = 1, to = length(samples), by = freq)]
  } else {
    stop(paste(
      "Input 'samples' must be a numeric vector,",
      "matrix or list of vectors and matrices!"
    ))
  }
}



#' @title Remove samples of a Markov Chain by a constant frequency `freq`
#'
#' @description
#' `thinning()` takes a numeric vector, matrix or list of matrices as input and
#' returns a modified copy of the input object with only a fraction of the
#' original entries.
#' The function is primarily used in the context of Markov-Chain Monte-Carlo
#' sampling to obtain uncorrelated samples of the posterior distribution.
#'
#' @param samples Numeric vector, matrix or list of matrices containing the
#'                samples of the MCMC sampler.
#'
#' @param freq Positive integer which controls the degree to which the
#'             Markov-Chain is thinned out, e.g. `freq = 10` only keeps every
#'             10th sample.
#'
#' @returns
#' Same data type as the provided input, i.e. a numeric vector, matrix or list
#' of matrices with thinned out elements.
#'
#' @export

thinning <- function(samples, freq) {
  if (is.matrix(samples) || is.vector(samples, mode = "numeric")) {
    thin_element(samples = samples, freq = freq)
  } else if (is.list(samples)) {
    lapply(X = samples, thin_element, freq = freq)
  } else {
    stop(paste(
      "Input 'samples' must be a numeric vector,",
      "matrix or list of vectors and matrices!"
    ))
  }
}
