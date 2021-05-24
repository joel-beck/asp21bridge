#<<<<<<< HEAD
#' @title Gibbs sampling algorithm
#=======
# helper function for gibbs_sampler
# detects intercept in design matrix
# later use: intercept_x <- includes_intercept(X)
#            intercept_z <- includes_intercept(Z)

includes_intercept <- function(mat) {
  any(
    apply(
      X = mat,
      MARGIN = 2,
      FUN = function(col) all(col == col[1])
    )
  )
}





#' @title Gibbs sampling algorithm
#>>>>>>> 35ceac7d5515364d2bd6325b2e8b71b801b232ff
#'
#' @description The `gibbs_sampler()` function is a Markov chain Monte Carlo (MCMC)
#'              algorithm allowing to simulate from parameters of a location scale
#'              regression model or from matrices and parameters that are assigned by
#'              hand. \cr
#'              The function is primarly used in a Bayesian setting, where
#'              ridge-regularization can be obtained through the variance of the
#'              coefficient distributions.
#'
#'
#' @param m Model object containing a beta and gamma predictor, as well as the according
#'          model matrices. The model is used in second stage, if one of `X`, `Z`, `y`,
#'          `beta_start` or `gamma_start` is not given, Default: NULL
#' @param X Matrix containing the data according to the beta_start coefficient,
#'          Default: NULL
#' @param Z Matrix containing the data according to the gamma_start coefficient,
#'          Default: NULL
#' @param y Response vector; `y = X * beta + Z * gamma`, Default: NULL
#' @param beta_start Starting vector for simulation, where `beta` is a linear predictor
#'                   for the mean (= the location), that is normally distributed,
#'                   Default: NULL
#' @param gamma_start Starting vector for simulation, where `gamma` is a linear predictor
#'                    for the standard deviation (= the scale), that is normally
#'                    distributed, Default: NULL
#' @param tau_start Starting value for the variance of the normal distribution of the beta
#'                  parameter, where `tau_start` is inverse gamma (IG) distributed.
#'                  Regularization parameter in a Bayesian ridge setting, Default: 1
#' @param xi_start Starting value for the variance of the normal distribution of the gamma
#'                 parameter, where `xi_start` is inverse gamma (IG) distributed.
#'                 Regularization parameter in a Bayesian ridge setting, Default: 1
#' @param a_tau Fix shape parameter of the IG distribution of `tau_start`, Default: 1
#' @param b_tau Fix scale parameter of the IG distribution of `tau_start`, Default: 3
#' @param a_xi Fix shape parameter of the IG distribution of `xi_start`, Default: 1
#' @param b_xi Fix scale parameter of the IG distribution of `xi_start`, Default: 3
#' @param prop_var Variance of proposal distribution for gamma sampling, Default: 3
#' @param num_sim Number of simulations, Default: 1000
#'
#' @return Depending on the input structure, two different output types are available.
#'         In each case a a Markov Chain of samples for the parameters is generated.
#'         If matrices and parameters are only assigned by hand, a list containing the
#'         parameters with their according iterations is available via the `print()`
#'         command.
#'         In case a model object is used at any stage, the output extends the model by
#'         the same list, that can be obtained via `model$mcmc_ridge` or outputs the last
#'         iteration as model parameters via `summary(model, type = "mcmc_ridge)`.
#'
#' @examples
#'
#' # Gibbs sampling with lslm model input.
#' # Output extends model by Ridge-part.
#' fit <- lslm(
#'   location = y ~ x1 + x2 + z1 + z2, scale = ~ z1 + z2,
#'   data = toy_data, light = FALSE
#' )
#' fit <- gibbs_sampler(
#'   m = fit,
#'   num_sim = 1000
#' )
#' summary(fit, type = "mcmc_ridge")
#'
#' # Gibbs sampling with input by hand.
#' # Outputs list with each parameter and iteration.
#' beta <- c(2, -4, 5, 1)
#' gamma <- c(0.003, 0.002)
#' fit <- gibbs_sampler(
#'   X = cbind(toy_data$x1, toy_data$x2, toy_data$z1, toy_data$z2),
#'   Z = cbind(toy_data$z1, toy_data$z2), y = toy_data$y,
#'   beta_start = beta, gamma_start = gamma,
#'   tau_start = 3, xi_start = 0.5,
#'   prop_var = 2.3
#' )
#' print(fit)
#'
#' # Gibbs sampling with mixed input. The function uses matrices assigned by
#' # hand on first stage and takes missing parameters out of the model.
#' # Output extends model by Ridge-part, since model objects were used.
#' gamma <- c(0.003, 0.002)
#' fit <- lslm(
#'   location = y ~ x1 + x2 + z1 + z2, scale = ~ z1 + z2,
#'   data = toy_data, light = FALSE
#' )
#' fit <- gibbs_sampler(
#'   m = fit,
#'   Z = cbind(toy_data$z1, toy_data$z2),
#'   gamma_start = gamma,
#'   num_sim = 10000
#' )
#' summary(fit, type = "mcmc_ridge")
#'
#' @export

gibbs_sampler <- function(m = NULL, X = NULL, Z = NULL, y = NULL,
                          beta_start = NULL, gamma_start = NULL, tau_start = 1, xi_start = 1,
                          a_tau = 1, b_tau = 3, a_xi = 1, b_xi = 3, prop_var = 3,
                          num_sim = 1000) {
  mod <- FALSE
  mcmc_ridge_m <- m

  if (is.null(m) &
    (is.null(X) | is.null(Z) | is.null(y) | is.null(beta_start) | is.null(gamma_start))) {
    stop("At least either all model matrices (X, Z, y) and coefficients (beta_start, gamma_start)
         or a model object (m) must be given.")
  } else if (is.null(X) | is.null(Z) | is.null(y) | is.null(beta_start) | is.null(gamma_start)) {
    input_list <- list(
      X = X,
      Z = Z,
      y = y,
      beta_start = beta_start,
      gamma_start = gamma_start
    )
    mod_name_list <- list(
      X = mcmc_ridge_m$x,
      Z = mcmc_ridge_m$z,
      Y = mcmc_ridge_m$y,
      beta_start = mcmc_ridge_m$coefficients$location, # hard-gecoded
      gamma_start = mcmc_ridge_m$coefficients$scale # hard-gecoded
    )
    for (l in 1:length(input_list)) {
      if (is.null(input_list[[l]])) {
        mod <- TRUE
        assign(names(input_list[l]), mod_name_list[[l]])
      }
    }

    if ((ncol(X) != length(beta_start)) | (ncol(Z) != length(gamma_start))) {
      stop("Dimensions of design matrices do not match with length of coefficients.")
    }
  } else if ((ncol(X) != length(beta_start)) | (ncol(Z) != length(gamma_start))) {
    stop("Dimensions of design matrices do not match with length of coefficients.")
  }


  n <- length(y) # number of observations
  acceptance_count <- 0

  beta_samples <- matrix(rep(0, times = num_sim * length(beta_start)),
    nrow = num_sim, ncol = length(beta_start)
  )
  colnames(beta_samples) <- paste("beta", 0:(ncol(beta_samples) - 1), sep = "_")

  gamma_samples <- matrix(rep(0, times = num_sim * length(gamma_start)),
    nrow = num_sim, ncol = length(gamma_start)
  )
  colnames(gamma_samples) <- paste("gamma", 0:(ncol(gamma_samples) - 1), sep = "_")

  tau_samples <- matrix(rep(0, times = num_sim * length(tau_start)),
    nrow = num_sim, ncol = length(tau_start)
  )
  colnames(tau_samples) <- "tau^2"

  xi_samples <- matrix(rep(0, times = num_sim * length(xi_start)),
    nrow = num_sim, ncol = length(xi_start)
  )
  colnames(xi_samples) <- "xi^2"

  # Eintragen der Startwerte in die erste Zeile
  beta_samples[1, ] <- beta_start
  gamma_samples[1, ] <- gamma_start
  tau_samples[1, ] <- tau_start
  xi_samples[1, ] <- xi_start

  # Berechnung der Variablen für die Full Conditionals
  K <- length(beta_start)
  J <- length(gamma_start)
  W <- matrix(rep(0, times = n * K), nrow = n)
  u <- numeric(length = n)


  for (i in 2:num_sim) {


    # sampling beta
    g_gamma <- vector(mode = "numeric", length = n)
    for (k in 1:n) {
      g_gamma[k] <- exp(sum(Z[k, ] * gamma_samples[i - 1, ]))
      W[k, ] <- X[k, ] / g_gamma[k]
      u[k] <- y[k] / g_gamma[k]
    }


    beta_var <- solve(crossprod(W) + (1 / tau_start^2) * diag(K))
    beta_mean <- beta_var %*% crossprod(W, u)
    beta_samples[i, ] <- mvtnorm::rmvnorm(n = 1, mean = beta_mean, sigma = beta_var)

    # sampling gamma with metropolis-hastings

    gamma_list <- mh_gamma(
      y = y, X = X, Z = Z, beta = beta_samples[i, ],
      gamma = gamma_samples[i - 1, ], g_gamma = g_gamma,
      xi_squared = xi_samples[i - 1, ], prop_var = prop_var
    )

    gamma_samples[i, ] <- gamma_list$gamma
    acceptance_count <- acceptance_count + gamma_list$accepted

    # sampling tau
    tau_samples[i, ] <- 1 / stats::rgamma(
      n = ncol(tau_samples),
      shape = a_tau + K / 2,
      scale = b_tau + 0.5 * crossprod(beta_samples[i - 1, ])
    )

    # sampling xi
    xi_samples[i, ] <- 1 / stats::rgamma(
      n = ncol(xi_samples),
      shape = a_xi + J / 2,
      scale = b_xi + 0.5 * crossprod(gamma_samples[i - 1, ])
    )
  }

  if (mod) {
    result_gibbs_list <- list(
      coefficient_sampling_matrices = list(
        location = beta_samples,
        scale = gamma_samples,
        location_prior = tau_samples,
        scale_prior = xi_samples
      ),
      acceptance_rate = acceptance_count / num_sim
    )

    m$mcmc_ridge <- result_gibbs_list
    return(m)
  } else {
    result_gibbs_list <- list(
      sampling_matrices = list(
        beta_samples = beta_samples,
        gamma_samples = gamma_samples,
        tau_samples = tau_samples,
        xi_samples = xi_samples
      ),
      acceptance_rate = acceptance_count / num_sim
    )

    return(result_gibbs_list)
  }
}
