pacman::p_load(
  dplyr,
  purrr,
  ggplot2
)
library(asp21bridge)


### Simulate Data


### Data generation

set.seed(11)

rep <- 100

array_num_sim <- c(100,300,500,1000)

num_sim_length <- length(array_num_sim)

mean_mean <- matrix(0, nrow = 6, ncol = num_sim_length)
colnames(mean_mean) <- array_num_sim
rownames(mean_mean) <- c("beta_0","beta_1","beta_2","gamma_0",
                         "gamma_1", "gamma_2")

mean_absolute_error <- matrix(0, nrow = 6, ncol = num_sim_length)
colnames(mean_absolute_error) <- array_num_sim
rownames(mean_absolute_error) <- c("beta_0","beta_1","beta_2","gamma_0",
                                   "gamma_1", "gamma_2")

mean_squared_error <- matrix(0, nrow = 6, ncol = num_sim_length)
colnames(mean_squared_error) <- array_num_sim
rownames(mean_squared_error) <- c("beta_0","beta_1","beta_2","gamma_0",
                                  "gamma_1", "gamma_2")

mean_of_variances_within <- matrix(0, nrow = 6, ncol = num_sim_length)
colnames(mean_of_variances_within) <- array_num_sim
rownames(mean_of_variances_within) <- c("beta_0","beta_1","beta_2","gamma_0",
                                  "gamma_1", "gamma_2")

mean_array <- array(0, dim = c(6,num_sim_length,rep))
absolute_error_array <- array(0, dim = c(6,num_sim_length,rep))
squared_error_array <- array(0, dim = c(6,num_sim_length,rep))
variances_within_samples_array <- array(0, dim = c(6,num_sim_length,rep))


for(r in 1:rep){

    n <- 100 # Sample size = 100

    x0 <- rep(1,n)
    x1 <- rnorm(n, mean = 1, sd = 1)
    x2 <- rnorm(n, mean = 2, sd = 1)
    z0 <- rep(1,n)
    z1 <- rnorm(n, mean = 1, sd = 1)
    z2 <- rnorm(n, mean = 2, sd = 1)

    X <- cbind(x0, x1, x2)
    Z <- cbind(z0, z1, z2)

    beta <- c(1, -1, 4) # True: beta_0 = 1, beta_1 = -1, beta_2 = 4,
    gamma <- c(0, -0.5, 0.5) # True: gamma_0 = 2, gamma_1 = -2, gamma_2 = 1

    y <- vector(mode = "numeric", length = n)
    for (i in seq_along(y)) {
      mu <- sum(X[i,] * beta)
      sigma <- exp(sum(Z[i,] * gamma))
      y[i] <- rnorm(n = 1, mean = mu, sd = sigma)
    }

    data <- as_tibble(cbind(y,X,Z))

    model <- lmls(
      location = y ~ x1 + x2, scale = ~ z1 + z2, data = data,
      light = FALSE
    )

    fit <- mcmc_ridge(m = model)

    for (j in 1:num_sim_length){

    samples <- mcmc_ridge(m = model, X = X, Z=Z, y=y, num_sim = array_num_sim[j],
                          beta_start = fit$coefficients$location,
                          gamma_start = fit$coefficients$scale
    )

    mean_array[,j,r] <- c(colMeans(samples$sampling_matrices$beta_samples),
                          colMeans(samples$sampling_matrices$gamma_samples))

    absolute_error_array[,j,r] <- abs(mean_array[,j,r] - c(beta,gamma))
    squared_error_array[,j,r] <- (absolute_error_array[,j,r])^2
    variances_within_samples_array[,j,r] <- c(diag(var(samples$sampling_matrices$beta_samples)),
                                              diag(var(samples$sampling_matrices$gamma_samples)))


  }
}

for (i in 1:6){
  for (j in 1:num_sim_length){
    mean_mean[i,j] <- mean(mean_array[i,j,])
    mean_absolute_error[i,j] <- mean(absolute_error_array[i,j,])
    mean_squared_error[i,j] <- mean(squared_error_array[i,j,])
    mean_of_variances_within[i,j] <- mean(variances_within_samples_array[i,j,])

  }
}

matplot(t(mean_mean), x = array_num_sim, type = "b",pch=1,col = 1:6,
        main = " Mean of Posterior Means",
        xlab = "Number of MCMC Ridge Simulations", ylab = "")
legend(x = "topright",legend = rownames(mean_mean),
       col=1:6, pch=1, bty = "n")

matplot(t(mean_absolute_error), x = array_num_sim, type = "b",pch=1,col = 1:6,
        main = " MAE of Posterior Means",
        xlab = "Number of MCMC Ridge Simulations", ylab = "")
legend(x = "topright",legend = rownames(mean_absolute_error),
       col=1:6, pch=1, bty = "n")

matplot(t(mean_squared_error), x = array_num_sim, type = "b",pch=1,col = 1:6,
        main = " MSE of Posterior Means",
        xlab = "Number of MCMC Ridge Simulations", ylab = "")
legend(x = "topright",legend = rownames(mean_squared_error),
       col=1:6, pch=1, bty = "n")

matplot(t(mean_of_variances_within), x = array_num_sim, type = "b",pch=1,col = 1:6,
        main = " Mean of Variances within the Samples",
        xlab = "Number of MCMC Ridge Simulations", ylab = "")
legend(x = "topright",legend = rownames(mean_of_variances_within),
       col=1:6, pch=1, bty = "n")

#   ____________________________________________________________________________
#   Data for Second Report                                                  ####

readr::write_rds(
  x = list(
    mean_mean = mean_mean,
    mean_absolute_error = mean_absolute_error,
    mean_squared_error = mean_squared_error,
    mean_of_variances_within = mean_of_variances_within,
    array_num_sim = array_num_sim
  ),
  file = here::here("simulation-studies", "number_mcmc_ridge_simulations.rds")
)

