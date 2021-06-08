# Hilfsfunktion Current (letzter Teil der Full Conditional für Gamma)
help_function_1 <- function(y, X, Z, beta, g_gamma, n) {
  sum <- 0

  for (i in 1:n) {
    sum <- sum + ((y[i] - sum(X[i, ] * beta)) / g_gamma[i])^2
  }
  return(sum)
}

# Log Full Conditional für aktuelles Gamma
log_full_cond_1 <- function(y, X, Z, beta, gamma, g_gamma, xi_squared, n) {
  -0.5 * ((1 / xi_squared) * sum(gamma[2:length(gamma)]^2) +
    2 * sum(Z %*% drop(gamma)) +
    help_function_1(
      y = y, X = X, Z = Z, beta = beta, g_gamma = g_gamma, n = n
    ))
}

# Hilfsfunktion Proposal (letzter Teil der Full Conditional für Gamma)
help_function_2 <- function(y, X, Z, beta, gamma, n) {
  sum <- 0

  for (i in 1:n) {
    sum <- sum + ((y[i] - sum(X[i, ] * beta)) / exp(sum(Z[i, ] * gamma)))^2
  }
  return(sum)
}

# Log Full Conditional für Proposal-Gamma
log_full_cond_2 <- function(y, X, Z, beta, gamma, xi_squared, n) {
  -0.5 * ((1 / xi_squared) * sum(gamma[2:length(gamma)]^2)
    + 2 * sum(Z %*% drop(gamma)) +
    help_function_2(y = y, X = X, Z = Z, beta = beta, gamma = gamma, n = n))
}

################ Gamma Update #################################


mh_gamma <- function(y, X, Z, beta, gamma, g_gamma, xi_squared, prop_var_scale) {
  n <- length(y)

  # Proposal für Gamma
  Sigma <- diag(1 / colMeans(Z^2)) * prop_var_scale / (n * length(gamma))
  gamma_proposal <- mvtnorm::rmvnorm(n = 1, mean = gamma, sigma = Sigma)

  # Auswertung Full Conditional Densities
  log_full_cond_curr <- log_full_cond_1(
    y = y, X = X, Z = Z, beta = beta, gamma = gamma, g_gamma = g_gamma,
    xi_squared = xi_squared, n = n
  )

  log_full_cond_proposal <- log_full_cond_2(
    y = y, X = X, Z = Z, beta = beta, gamma = gamma_proposal,
    xi_squared = xi_squared, n = n
  )

  # Berechnung Acceptance Probability
  log_acceptance_prob <- min(0, log_full_cond_proposal - log_full_cond_curr)
  acceptance_prob <- exp(log_acceptance_prob)


  # Update
  if (stats::runif(1) < acceptance_prob) {
    new_gamma <- gamma_proposal
    accepted <- TRUE
  } else {
    new_gamma <- gamma
    accepted <- FALSE
  }

  return(
    list(gamma = new_gamma, accepted = accepted)
  )
}
