
# Log Full Conditional für Beta
log_full_cond <- function(beta, beta_var_inv, W, u) {
  -0.5 * (sum(crossprod(beta, beta_var_inv) * beta)) - 2 * sum(beta * crossprod(W, u))
}

################ Beta Update #########################################

mh_beta <- function(beta, W, u, tau_squared, n, prop_var_loc) {
  K <- length(beta) - 1
  beta_var_inv <- crossprod(W) + (1 / tau_squared) * diag(c(0, rep(1, times = K)))
  beta_var <- solve(beta_var_inv)

  # Proposal für Beta
  Sigma <- beta_var * prop_var_loc / ((2 * n) * length(beta))
  beta_proposal <- drop(mvtnorm::rmvnorm(n = 1, mean = beta, sigma = Sigma))

  # Auswertung Full Conditional Densities
  log_full_cond_curr <- log_full_cond(
    beta = beta, beta_var_inv = beta_var_inv, W = W, u = u
  )

  log_full_cond_proposal <- log_full_cond(
    beta = beta_proposal, beta_var_inv = beta_var_inv, W = W, u = u
  )

  # Berechnung Acceptance Probability
  log_acceptance_prob <- min(0, log_full_cond_proposal - log_full_cond_curr)
  acceptance_prob <- exp(log_acceptance_prob)

  # Update
  if (stats::runif(1) < acceptance_prob) {
    new_beta <- beta_proposal
    accepted <- TRUE
  } else {
    new_beta <- beta
    accepted <- FALSE
  }

  return(
    list(beta = new_beta, accepted = accepted)
  )
}
