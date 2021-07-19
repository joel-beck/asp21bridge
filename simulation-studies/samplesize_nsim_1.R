pacman::p_load(
  dplyr,
  purrr,
  ggplot2
)
library(asp21bridge)


### Simulate Data


### Data generation

set.seed(10)

rep <- 30

n_data <- c(30,50,100,200,300)
array_num_sim <- c(100,300,500,1000)

number_of_n <- length(n_data)

mean <- matrix(0, nrow = 6, ncol = number_of_n)
colnames(mean) <- n_data
rownames(mean) <- c("beta_0","beta_1","beta_2","gamma_0",
                    "gamma_1", "gamma_2")

abs_dev <- matrix(0, nrow = 6, ncol = number_of_n)
colnames(abs_dev) <- n_data
rownames(abs_dev) <- c("beta_0","beta_1","beta_2","gamma_0",
                       "gamma_1", "gamma_2")

variance <- matrix(0, nrow = 6, ncol = number_of_n)
colnames(variance) <- n_data
rownames(variance) <- c("beta_0","beta_1","beta_2","gamma_0",
                        "gamma_1", "gamma_2")

coef_array <- array(0, dim = c(6,number_of_n,rep))
abs_dev_array <- array(0, dim = c(6,number_of_n,rep))
variance_array <- array(0, dim = c(6,number_of_n,rep))

for(s in 1:length(array_num_sim)){

  for(r in 1:rep){

    for(j in 1:number_of_n){

      n <- n_data[j]

      x1 <- rnorm(n, mean = 1, sd = 1)
      x2 <- rnorm(n, mean = 2, sd = 1)
      z1 <- rnorm(n, mean = 1, sd = 1)
      z2 <- rnorm(n, mean = 2, sd = 1)

      X <- cbind(x1, x2)
      Z <- cbind(z1, z2)

      beta <- c(-1, 4) # True: beta_0 = 0, beta_1 = -1, beta_2 = 4,
      gamma <- c(-0.5, 1) # True: gamma_0 = 0, gamma_1 = -0.5, gamma_2 = 1

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

      fit <- mcmc_ridge(m = model, num_sim = array_num_sim[s])

      coef_array[,j,r] <- c(fit$coefficients$location, fit$coefficients$scale)

      abs_dev_array[,j,r] <- abs(coef_array[,j,r] - c(0,beta,0,gamma))

    }
  }

  for (i in 1:6){
    for (j in 1:number_of_n){
      mean[i,j] <- mean(coef_array[i,j,])
      abs_dev[i,j] <- mean(abs_dev_array[i,j,])
      variance[i,j] <- var(coef_array[i,j,])
    }
  }

  matplot(t(mean), x = n_data, type = "b",pch=1,col = 1:nrow(mean),
          main = paste("Posterior Mean, num_sim = ",array_num_sim[s]),
          xlab = "Sample size", ylab = "")
  legend(x = "topright",legend = rownames(mean),
         col=1:nrow(mean), pch=1, bty = "n")

  matplot(t(abs_dev), x = n_data, type = "b",pch=1,col = 1:nrow(abs_dev),
          main = paste("Posterior Absolute Deviation, num_sim = ",array_num_sim[s]),
          xlab = "Sample size", ylab = "")
  legend(x = "topright",legend = rownames(abs_dev),
         col=1:nrow(abs_dev), pch=1, bty = "n")

  matplot(t(variance), x = n_data, type = "b",pch=1,col = 1:nrow(variance),
          main = paste("Posterior Variance, num_sim = ",array_num_sim[s]),
          xlab = "Sample size", ylab = "")
  legend(x = "topright",legend = rownames(variance),
         col=1:nrow(variance), pch=1, bty = "n")



}
