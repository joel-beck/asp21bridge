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

# initializes sampling matrices with correct dimension and column names
# preallocates each row with the given starting value
init_sampling_matrix <- function(name, nrow, start_value) {
  mat <- matrix(
    data = rep(start_value, times = nrow),
    nrow = nrow, byrow = TRUE
  )
  if (length(start_value) == 1) {
    colnames(mat) <- name
  } else if (length(start_value) > 1) {
    colnames(mat) <- paste(name, 0:(length(start_value) - 1), sep = "_")
  } else {
    stop("'start_value' must be a scalar or a nonempty vector!")
  }
  return(mat)
}


# creates output list depending on logicals 'mod' and 'mh_location'
create_output <- function(mod, beta_samples, gamma_samples, tau_samples, xi_samples,
                          num_sim, acc_count_scale, mh_location, acc_count_loc) {
  if (mod) {
    result_gibbs_list <- list(
      sampling_matrices = list(
        location = beta_samples,
        scale = gamma_samples,
        location_prior = tau_samples,
        scale_prior = xi_samples
      )
    )
  } else {
    result_gibbs_list <- list(
      sampling_matrices = list(
        beta_samples = beta_samples,
        gamma_samples = gamma_samples,
        tau_samples = tau_samples,
        xi_samples = xi_samples
      )
    )
  }
  if (mh_location) {
    result_gibbs_list$acceptance_rates <- list(
      acc_rate_loc = acc_count_loc / num_sim,
      acc_rate_scale = acc_count_scale / num_sim
    )
  } else {
    result_gibbs_list$acceptance_rate <- acc_count_scale / num_sim
  }
  return(result_gibbs_list)
}


# Validate function's input conditions
validate_input <- function(m, X, Z, y,
                           beta_start, gamma_start,
                           mod, mcmc_ridge_m) {
  if (is.null(m) &
      (is.null(X) | is.null(Z) | is.null(y) | is.null(beta_start) | is.null(gamma_start))) {
    stop(paste(
      "At least either all model matrices (X, Z, y) and coefficients",
      "(beta_start, gamma_start) or a model object (m) must be given."
    ))
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
      beta_start = mcmc_ridge_m$coefficients$location,
      gamma_start = mcmc_ridge_m$coefficients$scale
    )
    for (l in 1:length(input_list)) {
      if (is.null(input_list[[l]])) {
        mod <- TRUE
        assign(names(input_list[l]), mod_name_list[[l]])
        input_list[[l]] <- mod_name_list[[l]]
      }
    }

    if ((ncol(X) != length(beta_start)) | (ncol(Z) != length(gamma_start))) {
      stop("Dimensions of design matrices do not match with length of coefficients.")
    }
    input_list <- append(input_list, mod)
    return(input_list)

  } else if ((ncol(X) != length(beta_start)) | (ncol(Z) != length(gamma_start))) {
    stop("Dimensions of design matrices do not match with length of coefficients.")
  }
}

