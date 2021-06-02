# first skeleton to work with gibbs_sampler()

mh_beta <- function(y, X, Z, beta, gamma, g_gamma, xi_squared, prop_var_loc) {
  new_beta <- beta
  accepted <- FALSE

  return(
    list(beta = new_beta, accepted = accepted)
  )
}
