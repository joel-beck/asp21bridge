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



#' @title Remove samples of a Markov Chain by a constant frequency
#'
#' @description `thinning()` takes a numeric vector, matrix or list of matrices
#'              as input and returns a modified copy of the input object with
#'              only a fraction of the original entries. \cr
#'              The function is primarily used in the context of Markov-Chain
#'              Monte-Carlo sampling to obtain uncorrelated samples of the
#'              posterior distribution.
#'
#' @param samples Numeric vector, matrix or list of matrices containing the
#'                samples of the MCMC sampler.
#'
#' @param freq Positive integer which controls the degree to which the
#'             Markov-Chain is thinned out, e.g. `freq = 10` only keeps every
#'             10th sample.
#'
#' @returns Same data type as the provided input, i.e. a numeric vector,
#'          matrix or list of matrices with thinned out elements.
#'
#' @examples
#' fit <- lslm(
#'   location = y ~ x1 + x2 + z1 + z2, scale = ~ z1 + z2,
#'   data = toy_data, light = FALSE
#' ) %>%
#'   gibbs_sampler(num_sim = 1000)
#'
#' # list of 4 matrices with 1000 rows each
#' samples <- fit$mcmc_ridge$coefficient_sampling_matrices
#'
#' # list input applies thinning() to all matrices
#' thinning(samples, freq = 10)
#'
#' # extract matrices separately for different thinning frequencies
#' thinning(samples$location, freq = 5)
#' thinning(samples$scale, freq = 20)
#'
#' # works in many different contexts
#' thinning(rnorm(50), freq = 5)
#'
#' # works well in combination with burnin() and in
#' # pipe workflows
#' samples$location %>%
#'   burnin(num_burn = 10) %>%
#'   thinning(freq = 5) %>%
#'   mult_plot(type = "time", latex = TRUE, free_scale = TRUE)
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
