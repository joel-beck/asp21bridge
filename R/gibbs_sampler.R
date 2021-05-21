#' @title FUNCTION_TITLE
#'
#' @description FUNCTION_DESCRIPTION
#'
#' @param m PARAM_DESCRIPTION, Default: NULL
#' @param X PARAM_DESCRIPTION, Default: NULL
#' @param Z PARAM_DESCRIPTION, Default: NULL
#' @param y PARAM_DESCRIPTION, Default: NULL
#' @param beta_start PARAM_DESCRIPTION, Default: NULL
#' @param gamma_start PARAM_DESCRIPTION, Default: NULL
#' @param tau_start PARAM_DESCRIPTION, Default: 1
#' @param xi_start PARAM_DESCRIPTION, Default: 1
#' @param a_tau PARAM_DESCRIPTION, Default: 1
#' @param b_tau PARAM_DESCRIPTION, Default: 3
#' @param a_xi PARAM_DESCRIPTION, Default: 1
#' @param b_xi PARAM_DESCRIPTION, Default: 3
#' @param prop_var PARAM_DESCRIPTION, Default: 1
#' @param num_sim PARAM_DESCRIPTION, Default: 1000
#'
#' @return OUTPUT_DESCRIPTION
#'
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # EXAMPLE1
#' }
#' }
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
