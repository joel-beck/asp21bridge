#' @title Gibbs sampling algorithm
#'
#' @description The `mcmc_ridge()` function is a Markov chain Monte Carlo (MCMC)
#'              algorithm allowing to simulate from parameters of a location scale
#'              regression model or from matrices and parameters that are assigned by
#'              hand. \cr
#'              The function is primarily used in a Bayesian setting, where
#'              ridge-regularization can be obtained through the variance of the
#'              corresponding prior distributions.
#'
#' @param m Model object containing a beta and gamma predictor, as well as the according
#'          model matrices. The model is used in second stage, if one of `X`, `Z`, `y`,
#'          `beta_start` or `gamma_start` is not given, \cr
#'          Default: NULL
#' @param X Matrix containing the data according to the beta_start coefficient, \cr
#'          Default: NULL
#' @param Z Matrix containing the data according to the gamma_start coefficient, \cr
#'          Default: NULL
#' @param y Response vector; `y = X * beta + Z * gamma`,
#'          Default: NULL
#' @param num_sim Number of simulations, \cr
#'                Default: 1000
#' @param beta_start Starting vector for simulation, where `beta` is a linear predictor
#'                   for the mean (= the location), that is normally distributed, \cr
#'                   Default: NULL
#' @param gamma_start Starting vector for simulation, where `gamma` is a linear predictor
#'                    for the standard deviation (= the scale), that is normally
#'                    distributed, \cr
#'                    Default: NULL
#' @param tau_start Starting value for the variance of the normal distribution of the beta
#'                  parameter, where `tau_start` is inverse gamma (IG) distributed.
#'                  Regularization parameter in a Bayesian ridge setting, \cr
#'                  Default: 1
#' @param xi_start Starting value for the variance of the normal distribution of the gamma
#'                 parameter, where `xi_start` is inverse gamma (IG) distributed.
#'                 Regularization parameter in a Bayesian ridge setting, \cr
#'                 Default: 1
#' @param a_tau Fix shape parameter of the IG distribution of `tau_start`, \cr
#'              Default: 1
#' @param b_tau Fix scale parameter of the IG distribution of `tau_start`, \cr
#'              Default: 3
#' @param a_xi Fix shape parameter of the IG distribution of `xi_start`, \cr
#'             Default: 1
#' @param b_xi Fix scale parameter of the IG distribution of `xi_start`, \cr
#'             Default: 3
#' @param prop_var_scale Variance of proposal distribution for `gamma` sampling, \cr
#'                       Default: 3
#' @param mh_location If TRUE, location parameter is sampled with
#'                    Metropolis - Hastings Algorithm. Can be used
#'                    with every kind of model input, \cr
#'                    Default: FALSE
#' @param prop_var_loc Variance of proposal distribution for `beta` sampling, \cr
#'                     Default: 200
#'
#' @return Depending on the input structure, different output types are available.
#'         In each case a a Markov Chain of samples for the parameters is generated. \cr
#'         If matrices and parameters are only assigned by hand, a list containing the
#'         parameters with their according iterations and the acceptance rate of a
#'         Metropolis-Hastings algorithm is available via the `print()` command. \cr
#'         In case a model object is used at any stage, the output extends the model by
#'         the same list, that can be obtained via `model$mcmc_ridge` or outputs the last
#'         iteration as model parameters via `summary(model, type = "mcmc_ridge)`. \cr
#'         In both cases, the package comes with the `summary_complete(model)` function,
#'         which gives additional information about the Markov Chain's results. For more
#'         see `?summary_complete()`.
#'
#' @examples
#' # Gibbs sampling with lmls model input.
#' # Output extends model by Ridge-part.
#' fit <- lmls(
#'   location = y ~ x1 + x2 + z1 + z2, scale = ~ z1 + z2,
#'   data = toy_data, light = FALSE
#' )
#' fit <- mcmc_ridge(
#'   m = fit,
#'   num_sim = 1000
#' )
#' summary(fit, type = "mcmc_ridge")
#' summary_complete(fit)
#'
#' # Sampling with by-hand-assigned-input and using Metropolis Hastings for location
#' # and scale parameter.
#' # Outputs list with each parameter and iteration as well as accuracy rates.
#' beta <- c(2, -4, 5, 1)
#' gamma <- c(0.003, 0.002)
#' fit <- mcmc_ridge(
#'   X = cbind(toy_data$x1, toy_data$x2, toy_data$z1, toy_data$z2),
#'   Z = cbind(toy_data$z1, toy_data$z2),
#'   y = toy_data$y,
#'   beta_start = beta, gamma_start = gamma,
#'   tau_start = 3, xi_start = 0.5,
#'   prop_var_scale = 2.3,
#'   mh_location = TRUE
#' )
#' print(fit)
#' summary_complete(fit)
#'
#' # Gibbs sampling with mixed input. The function uses matrices assigned by
#' # hand on first stage and takes missing parameters out of the model.
#' # Output extends model by Ridge-part, since model objects were used.
#' gamma <- c(0.003, 0.002)
#' fit <- lmls(
#'   location = y ~ x1 + x2 + z1 + z2, scale = ~ z1 + z2,
#'   data = toy_data, light = FALSE
#' )
#' fit <- mcmc_ridge(
#'   m = fit,
#'   Z = cbind(toy_data$z1, toy_data$z2),
#'   gamma_start = gamma,
#'   num_sim = 1000
#' )
#' summary(fit, type = "mcmc_ridge")
#' summary_complete(fit)
#' @export


mcmc_ridge <- function(m = NULL, X = NULL, Z = NULL, y = NULL, num_sim = 1000,
                          beta_start = NULL, gamma_start = NULL, tau_start = 1, xi_start = 1,
                          a_tau = 1, b_tau = 1, a_xi = 1, b_xi = 1,
                          prop_var_scale = 3, mh_location = FALSE, prop_var_loc = 200) {
  mod <- FALSE
  mcmc_ridge_m <- m

  # validate input ----------------------------------------------------------

  input_list <- validate_input(
    m = m, X = X, Z = Z, y = y,
    beta_start = beta_start, gamma_start = gamma_start,
    mod = mod, mcmc_ridge_m = mcmc_ridge_m
  )
  for (j in seq_along(input_list)) {
    assign(names(input_list[j]), input_list[[j]])
  }

  intercept_list <- add_intercept(
    X = X, Z = Z, beta_start = beta_start, gamma_start = gamma_start, mod = mod
  )

  for (j in seq_along(intercept_list)) {
    assign(names(intercept_list[j]), intercept_list[[j]])
  }

  # initialize variables ----------------------------------------------------

  n <- length(y)
  acc_count_scale <- 0
  if (mh_location) {
    acc_count_loc <- 0
  }

  # Berechnung der Variablen für die Full Conditionals
  K <- length(beta_start) - 1
  J <- length(gamma_start) - 1
  W <- matrix(rep(0, times = n * (K + 1)), nrow = n)
  u <- numeric(length = n)

  beta_samples <- init_sampling_matrix(
    name = "beta", nrow = num_sim, start_value = beta_start
  )
  gamma_samples <- init_sampling_matrix(
    name = "gamma", nrow = num_sim, start_value = gamma_start
  )
  tau_samples <- init_sampling_matrix(
    name = "tau^2", nrow = num_sim, start_value = tau_start
  )
  xi_samples <- init_sampling_matrix(
    name = "xi^2", nrow = num_sim, start_value = xi_start
  )

  # sampling process --------------------------------------------------------

  for (i in 2:num_sim) {
    g_gamma <- vector(mode = "numeric", length = n)
    for (k in 1:n) {
      g_gamma[k] <- exp(sum(Z[k, ] * gamma_samples[i - 1, ]))
      W[k, ] <- X[k, ] / g_gamma[k]
      u[k] <- y[k] / g_gamma[k]
    }

    if (mh_location) {
      # sampling beta with metropolis-hastings
      beta_list <- mh_beta(
        beta = beta_samples[i - 1, ], W = W, u = u,
        tau_squared = tau_samples[i - 1, ], n = n, prop_var_loc = prop_var_loc
      )
      beta_samples[i, ] <- beta_list$beta
      acc_count_loc <- acc_count_loc + beta_list$accepted
    } else {
      # sampling beta with closed form full conditional
      beta_var <- solve(crossprod(W) + (1 / tau_start^2) * diag(c(0, rep(1, times = K))))
      beta_mean <- beta_var %*% crossprod(W, u)
      beta_samples[i, ] <- mvtnorm::rmvnorm(n = 1, mean = beta_mean, sigma = beta_var)
    }

    # sampling gamma with metropolis-hastings
    gamma_list <- mh_gamma(
      y = y, X = X, Z = Z, beta = beta_samples[i, ],
      gamma = gamma_samples[i - 1, ], g_gamma = g_gamma,
      xi_squared = xi_samples[i - 1, ], prop_var_scale = prop_var_scale
    )

    gamma_samples[i, ] <- gamma_list$gamma
    acc_count_scale <- acc_count_scale + gamma_list$accepted

    # sampling tau
    tau_samples[i, ] <- invgamma::rinvgamma(
      n = ncol(tau_samples),
      shape = a_tau + K / 2,
      rate = b_tau + 0.5 * sum((beta_samples[i - 1, -1])^2)
    )

    # sampling xi
    xi_samples[i, ] <- invgamma::rinvgamma(
      n = ncol(xi_samples),
      shape = a_xi + J / 2,
      rate = b_xi + 0.5 * sum((gamma_samples[i - 1, -1])^2)
    )
  }

  # return value ------------------------------------------------------------

  result_gibbs_list <- create_output(
    mod, beta_samples, gamma_samples, tau_samples, xi_samples, num_sim,
    acc_count_scale, mh_location, acc_count_loc
  )

  if (mod) {
    m$mcmc_ridge <- result_gibbs_list
    return(m)
  } else {
    return(result_gibbs_list)
  }
}
