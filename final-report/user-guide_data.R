library(asp21bridge)
set.seed(1234)

toy_fit <- lmls(
    location = y ~ x1 + x2 + z1 + z2, scale = ~ z1 + z2,
    data = toy_data, light = FALSE
) %>%
    mcmc_ridge(num_sim = 10000)

y <- abdom$y
X <- as.matrix(abdom$x)
Z <- as.matrix(abdom$x)

abdom_fit <- mcmc_ridge(y = y, X = X, Z = Z, beta_start = 1, gamma_start = 1, num_sim = 10000)

readr::write_rds(x = list(toy_fit = toy_fit, abdom_fit = abdom_fit), file = here::here("final-report", "user-guide_data.rds"))
