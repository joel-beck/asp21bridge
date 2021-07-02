#   ____________________________________________________________________________
#   Helper Functions                                                        ####

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


#   ____________________________________________________________________________
#   Main Function                                                           ####

#' @title Remove beginning samples of a Markov Chain
#'
#' @description `burnin()` takes a numeric vector, matrix or list of matrices as
#'              input and returns a modified copy of the input object with a
#'              removed first part of the original entries. \cr
#'              The function is primarily used in the context of Markov-Chain
#'              Monte-Carlo sampling to obtain samples of the posterior
#'              distribution after a converging phase.
#'
#' @param samples Numeric vector, matrix or list of matrices containing the
#'                samples of the MCMC sampler.
#'
#' @param num_burn Nonnegative integer which controls the number of removed
#'                 samples, e.g. `num_burn = 10` removes the first 10 samples.
#'
#' @returns Same data type as the provided input, i.e. a numeric vector,
#'          matrix or list of matrices with removed elements at the beginning.
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
#' # list input applies burnin() to all matrices
#' burnin(samples, num_burn = 100)
#'
#' # extract matrices separately for different length of Burn-In phase
#' burnin(samples$location, num_burn = 50)
#' burnin(samples$scale, num_burn = 200)
#'
#' # works in many different contexts
#' burnin(rnorm(50), num_burn = 1)
#'
#' # works well in combination with thinning() and in
#' # pipe workflows
#' samples$location %>%
#'   burnin(num_burn = 100) %>%
#'   thinning(freq = 5) %>%
#'   mult_plot(type = "time", latex = TRUE, free_scale = TRUE)
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
