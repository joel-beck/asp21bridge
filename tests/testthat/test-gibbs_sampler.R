library(lslmbRidge)
library(testthat)

test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})

fit <- lslm(
  location = y ~ x1 + x2 + z1 + z2, scale = ~ z1 + z2,
  data = toy_data, light = FALSE)
X <- cbind(toy_data$x1, toy_data$x2, toy_data$z1, toy_data$z2)
Z <- cbind(toy_data$z1, toy_data$z2)
y <- toy_data$y
beta <- c(2, -4, 5, 1)
gamma <- c(0.003, 0.002)
testmat <- matrix(rnorm(600, 0, 1), ncol = 6)    # matrix to cause errors

fit1 <- gibbs_sampler(m = fit, num_sim = 50)     # model-only input
fit2 <- gibbs_sampler(X = cbind(toy_data$x1, toy_data$x2, toy_data$z1, toy_data$z2),
                      Z = cbind(toy_data$z1, toy_data$z2),
                      y = toy_data$y,
                      beta_start = beta, gamma_start = gamma,
                      num_sim = 50)    # by-hand input only
fit3 <- gibbs_sampler(m = fit,
                      Z = cbind(rep(1, times = length( toy_data$z1)), toy_data$z1, toy_data$z2),
                      beta_start = c(1, 2, 5, 3, 2),
                      num_sim = 50)    # mixed model

test_that("Using Model Matrices extends lslm-model summary() function", {
  expect_equal(summary(fit1, type = "mcmc_ridge")$mcmc_ridge,
               fit1$mcmc_ridge)
  expect_s3_class(fit1, "lslm")
})

test_that("Output list contains each parameter, each simulation and acceptance rate", {
  expect_equal(dim(fit1$mcmc_ridge$coefficient_sampling_matrices$location), c(50, 5))
  expect_equal(dim(fit1$mcmc_ridge$coefficient_sampling_matrices$scale), c(50, 3))
  expect_length(fit1$mcmc_ridge$coefficient_sampling_matrices$location_prior, 50)
  expect_length(fit1$mcmc_ridge$coefficient_sampling_matrices$scale_prior, 50)
  expect_type(fit1$mcmc_ridge$acceptance_rate, "double")
  expect_type(fit1$mcmc_ridge, "list")
})

test_that("Using by-hand-assigned-input does not extend the model, but gives a list", {
  expect_type(fit2$mcmc_ridge, "NULL")
  expect_type(fit2, "list")
  expect_type(summary(fit2, type = "mcmc_ridge"), "NULL") # Discussion point
})

test_that("Using mixed input extends lslm-model summary() function", {
  expect_equal(summary(fit3, type = "mcmc_ridge")$mcmc_ridge,
               fit3$mcmc_ridge)
  expect_s3_class(fit3, "lslm")
})

test_that("Error Code of not matching dimensions is thrown", {
  expect_error(gibbs_sampler(m = fit, X = testmat, num_sim = 50),
               "Dimensions of design matrices do not match with length of coefficients.")
  expect_error(gibbs_sampler(X = testmat, Z = cbind(toy_data$z1, toy_data$z2),
                             y = toy_data$y,
                             beta_start = beta, gamma_start = gamma,
                             num_sim = 50),
               "Dimensions of design matrices do not match with length of coefficients.")
})

test_that("Error Code of missing input data is thrown", {
  expect_error(gibbs_sampler(Z = cbind(toy_data$z1, toy_data$z2),
               y = toy_data$y,
               beta_start = beta, gamma_start = gamma,
               num_sim = 50),
               "At least either all model matrices (X, Z, y) and coefficients (beta_start, gamma_start)
               or a model object (m) must be given.")   # only problem: formatting
})
