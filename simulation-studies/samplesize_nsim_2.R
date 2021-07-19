pacman::p_load(
  dplyr,
  purrr,
  ggplot2
)
library(asp21bridge)


### Simulate Data


### Data generation

set.seed(2)

rep <- 30

n_data <- c(30,50,100,200,300)
array_num_sim <- c(100,300,500,1000)

number_of_n <- length(n_data)

mean <- matrix(0, nrow = 4, ncol = number_of_n)
colnames(mean) <- n_data
rownames(mean) <- c("beta_0","beta_1","gamma_0",
                    "gamma_1")

abs_dev <- matrix(0, nrow = 4, ncol = number_of_n)
colnames(abs_dev) <- n_data
rownames(abs_dev) <- c("beta_0","beta_1","gamma_0",
                       "gamma_1")

variance <- matrix(0, nrow = 4, ncol = number_of_n)
colnames(variance) <- n_data
rownames(variance) <- c("beta_0","beta_1","gamma_0",
                        "gamma_1")

mean_array <- array(0, dim = c(4,number_of_n,rep))
abs_dev_array <- array(0, dim = c(4,number_of_n,rep))
variance_array <- array(0, dim = c(4,number_of_n,rep))

for(s in 1:length(array_num_sim)){

  for(r in 1:rep){

    for(j in 1:number_of_n){

      n <- n_data[j]

      x0 <- rep(1,n)
      x1 <- rnorm(n, mean = 3, sd = 1)
      z0 <- rep(1,n)
      z1 <- rnorm(n, mean = 3, sd = 1)

      X <- cbind(x0, x1)
      Z <- cbind(z0, z1)

      beta <- c(3, -1) # True: beta_0 = 3, beta_1 = -1,
      gamma <- c(2, -0.5) # True: gamma_0 = 2, gamma_1 = -0.5

      y <- vector(mode = "numeric", length = n)
      for (i in seq_along(y)) {
        mu <- sum(X[i,] * beta)
        sigma <- exp(sum(Z[i,] * gamma))
        y[i] <- rnorm(n = 1, mean = mu, sd = sigma)
      }

      data <- as_tibble(cbind(y,X,Z))

      model <- lmls(
        location = y ~ x1, scale = ~ z1, data = data,
        light = FALSE
      )

      fit <- mcmc_ridge(m = model)

      samples <- mcmc_ridge(m = model, X = X, Z=Z, y=y, num_sim = array_num_sim[s],
                            beta_start = fit$coefficients$location,
                            gamma_start = fit$coefficients$scale
      )

      mean_array[,j,r] <- c(colMeans(samples$sampling_matrices$beta_samples),
                            colMeans(samples$sampling_matrices$gamma_samples))

      abs_dev_array[,j,r] <- abs(mean_array[,j,r] - c(beta,gamma))

      variance_array[,j,r] <- c(diag(var(samples$sampling_matrices$beta_samples)),
                                diag(var(samples$sampling_matrices$gamma_samples)))

    }
  }

  for (i in 1:4){
    for (j in 1:number_of_n){
      mean[i,j] <- mean(mean_array[i,j,])
      abs_dev[i,j] <- mean(abs_dev_array[i,j,])
      variance[i,j] <- var(mean_array[i,j,])
    }
  }

  matplot(t(mean), x = n_data, type = "b",pch=1,col = 1:nrow(mean),
          main = paste("Mean of Posterior Means, num_sim = ",array_num_sim[s]),
          xlab = "Sample size", ylab = "")
  legend(x = "topright",legend = rownames(mean),
         col=1:nrow(mean), pch=1, bty = "n")

  matplot(t(abs_dev), x = n_data, type = "b",pch=1,col = 1:nrow(abs_dev),
          main = paste("MAE of Posterior Means, num_sim = ",array_num_sim[s]),
          xlab = "Sample size", ylab = "")
  legend(x = "topright",legend = rownames(abs_dev),
         col=1:nrow(abs_dev), pch=1, bty = "n")

  matplot(t(variance), x = n_data, type = "b",pch=1,col = 1:nrow(variance),
          main = paste("Variance of Posterior Means, num_sim = ",array_num_sim[s]),
          xlab = "Sample size", ylab = "")
  legend(x = "topright",legend = rownames(variance),
         col=1:nrow(variance), pch=1, bty = "n")



}
