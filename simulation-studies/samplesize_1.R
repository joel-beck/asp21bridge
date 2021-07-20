pacman::p_load(
  dplyr,
  purrr,
  ggplot2
)
library(asp21bridge)


### Simulate Data


### Data generation

set.seed(1)

rep <- 100

n_data <- c(30,50,100,200,300,500)

number_of_n <- length(n_data)

mean_mean <- matrix(0, nrow = 6, ncol = number_of_n)
colnames(mean_mean) <- n_data
rownames(mean_mean) <- c("beta_0","beta_1","beta_2","gamma_0",
                    "gamma_1", "gamma_2")

mean_absolute_error <- matrix(0, nrow = 6, ncol = number_of_n)
colnames(mean_absolute_error) <- n_data
rownames(mean_absolute_error) <- c("beta_0","beta_1","beta_2","gamma_0",
                       "gamma_1", "gamma_2")

mean_squared_error <- matrix(0, nrow = 6, ncol = number_of_n)
colnames(mean_squared_error) <- n_data
rownames(mean_squared_error) <- c("beta_0","beta_1","beta_2","gamma_0",
                        "gamma_1", "gamma_2")

mean_array <- array(0, dim = c(6,number_of_n,rep))
absolute_error_array <- array(0, dim = c(6,number_of_n,rep))
squared_error_array <- array(0, dim = c(6,number_of_n,rep))
  
  for(r in 1:rep){
    
    for(j in 1:number_of_n){
      
      n <- n_data[j]
      
      x0 <- rep(1,n)
      x1 <- rnorm(n, mean = 1, sd = 1)
      x2 <- rnorm(n, mean = 2, sd = 1)
      z0 <- rep(1,n)
      z1 <- rnorm(n, mean = 1, sd = 1)
      z2 <- rnorm(n, mean = 2, sd = 1)
      
      X <- cbind(x0, x1, x2)
      Z <- cbind(z0, z1, z2)
      
      beta <- c(1, -1, 4) # True: beta_0 = 1, beta_1 = -1, beta_2 = 4,
      gamma <- c(0, -0.5, 1) # True: gamma_0 = 2, gamma_1 = -2, gamma_2 = 1
      
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
      
      samples <- mcmc_ridge(m = model, X = X, Z=Z, y=y, num_sim = 500,
                            beta_start = fit$coefficients$location,
                            gamma_start = fit$coefficients$scale
      )
      
      mean_array[,j,r] <- c(colMeans(samples$sampling_matrices$beta_samples),
                            colMeans(samples$sampling_matrices$gamma_samples))
      
      absolute_error_array[,j,r] <- abs(mean_array[,j,r] - c(beta,gamma))
      squared_error_array[,j,r] <- (absolute_error_array[,j,r])^2
      
    }
  }
  
  for (i in 1:6){
    for (j in 1:number_of_n){
      mean_mean[i,j] <- mean(mean_array[i,j,])
      mean_absolute_error[i,j] <- mean(absolute_error_array[i,j,])
      mean_squared_error[i,j] <- mean(squared_error_array[i,j,])
    }
  }
  
  matplot(t(mean_mean), x = n_data, type = "b",pch=1,col = 1:6,
          main = " Mean of Posterior Means",
          xlab = "Sample size", ylab = "")
  legend(x = "topright",legend = rownames(mean_mean),
         col=1:6, pch=1, bty = "n")
  
  matplot(t(mean_absolute_error), x = n_data, type = "b",pch=1,col = 1:6,
          main = " MAE of Posterior Means",
          xlab = "Sample size", ylab = "")
  legend(x = "topright",legend = rownames(mean_absolute_error),
         col=1:6, pch=1, bty = "n")
  
  matplot(t(mean_squared_error), x = n_data, type = "b",pch=1,col = 1:6,
          main = " MSE of Posterior Means",
          xlab = "Sample size", ylab = "")
  legend(x = "topright",legend = rownames(mean_squared_error),
         col=1:6, pch=1, bty = "n")
