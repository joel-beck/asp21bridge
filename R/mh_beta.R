
# Log Full Conditional für Beta
log_full_cond <- function(beta, beta_var, W, u) {
  -0.5 * (t(beta) %*% beta_var %*% beta - 2 * t(beta) %*% t(W) %*% u)
}

################ Beta Update #########################################

mh_beta <- function(beta, W, u, tau_squared, n, prop_var_loc) {

  K <- length(beta) - 1

  beta_var <- solve(crossprod(W) + (1 / tau_squared) * diag(c(0, rep(1, times = K))))

  # Proposal für Beta
  Sigma <- beta_var * prop_var_loc / (n * length(beta))
  beta_proposal <- drop(mvtnorm::rmvnorm(n = 1, mean = beta, sigma = Sigma))

  # Auswertung Full Conditional Densities
  log_full_cond_curr <- log_full_cond(
    beta = beta, beta_var = beta_var, W = W, u = u
  )

  log_full_cond_proposal <- log_full_cond(
    beta = beta_proposal, beta_var = beta_var, W = W, u = u
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
