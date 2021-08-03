set.seed(1234)

x1 <- rnorm(n = 50, mean = 1, sd = 1)
x2 <- rnorm(n = 50, mean = 2, sd = 1)

z1 <- rnorm(n = 50, mean = 5, sd = 1)
z2 <- rnorm(n = 50, mean = 3, sd = 1)

X <- cbind(x1, x2, z1, z2)
Z <- cbind(z1, z2)

beta <- c(-2, -1, 1, 2)
gamma <- c(-1, 1)

y <- vector(mode = "numeric", length = 50)

for (i in seq_along(y)) {
  mu <- sum(X[i,] * beta)
  sigma <- exp(sum(Z[i,] * gamma))
  y[i] <- rnorm(n = 1, mean = mu, sd = sigma)
}

toy_data <- dplyr::tibble(
  y = y, x1 = x1, x2 = x2, z1 = z1, z2 = z2
)

usethis::use_data(toy_data, overwrite = TRUE)
