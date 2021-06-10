#   ____________________________________________________________________________
#   Setup                                                                   ####

fit <- lslm(
  location = y ~ x1 + x2 + z1 + z2, scale = ~ z1 + z2,
  data = toy_data, light = FALSE
)
X <- cbind(toy_data$x1, toy_data$x2, toy_data$z1, toy_data$z2)
Z <- cbind(toy_data$z1, toy_data$z2)
y <- toy_data$y
beta <- c(2, -4, 5, 1)
gamma <- c(0.003, 0.002)

# matrix to cause errors
testmat <- matrix(rnorm(600, 0, 1), ncol = 6)

# model-only input
fit1 <- gibbs_sampler(m = fit, num_sim = 50)
invisible(capture.output(summary_fit1 <- summary(fit1, type = "mcmc_ridge")))

# by-hand input only
fit2 <- gibbs_sampler(
  X = cbind(toy_data$x1, toy_data$x2, toy_data$z1, toy_data$z2),
  Z = cbind(toy_data$z1, toy_data$z2),
  y = toy_data$y,
  beta_start = beta, gamma_start = gamma,
  num_sim = 50
)

# mixed model
fit3 <- gibbs_sampler(
  m = fit,
  Z = cbind(rep(1, times = length(toy_data$z1)), toy_data$z1, toy_data$z2),
  beta_start = c(1, 2, 5, 3, 2),
  num_sim = 50
)
invisible(capture.output(summary_fit3 <- summary(fit3, type = "mcmc_ridge")))

# by hand input with mh_location = TRUE
fit4 <- gibbs_sampler(
  X = cbind(toy_data$x1, toy_data$x2, toy_data$z1, toy_data$z2),
  Z = cbind(toy_data$z1, toy_data$z2),
  y = toy_data$y,
  beta_start = beta, gamma_start = gamma,
  mh_location = TRUE,
  num_sim = 50
)



#   ____________________________________________________________________________
#   Tests                                                                   ####

test_that("Using Model Matrices extends lslm-model summary() function", {
  expect_equal(
    summary_fit1$mcmc_ridge,
    fit1$mcmc_ridge
  )
  expect_s3_class(fit1, "lslm")
})

test_that("Output list contains each parameter, each simulation and acceptance rate", {
  expect_equal(dim(fit1$mcmc_ridge$sampling_matrices$location), c(50, 5))
  expect_equal(dim(fit1$mcmc_ridge$sampling_matrices$scale), c(50, 3))
  expect_length(fit1$mcmc_ridge$sampling_matrices$location_prior, 50)
  expect_length(fit1$mcmc_ridge$sampling_matrices$scale_prior, 50)
  expect_type(fit1$mcmc_ridge$acceptance_rate, "double")
  expect_type(fit1$mcmc_ridge, "list")
})

test_that("Using by-hand-assigned-input does not extend the model, but gives a list", {
  expect_type(fit2$mcmc_ridge, "NULL")
  expect_type(fit2, "list")
})

test_that("Using mixed input extends lslm-model summary() function", {
  expect_equal(
    summary_fit3$mcmc_ridge,
    fit3$mcmc_ridge
  )
  expect_s3_class(fit3, "lslm")
})

test_that("Independent of input type, the summary_complete() function can be called", {
  expect_equal(dim(summary_complete(fit1)), c(10, 6))
  expect_equal(dim(summary_complete(fit2)), c(10, 6))
  expect_equal(dim(summary_complete(fit3)), c(10, 6))
  expect_type(summary_complete(fit1), "list")
  expect_type(summary_complete(fit2), "list")
  expect_type(summary_complete(fit3), "list")
})

test_that("Error Code of not matching dimensions is thrown", {
  expect_error(
    gibbs_sampler(m = fit, X = testmat, num_sim = 50),
    "Dimensions of design matrices do not match with length of coefficients."
  )
  expect_error(
    gibbs_sampler(
      X = testmat, Z = cbind(toy_data$z1, toy_data$z2),
      y = toy_data$y,
      beta_start = beta, gamma_start = gamma,
      num_sim = 50
    ),
    "Dimensions of design matrices do not match with length of coefficients."
  )
})

test_that("Error Code of missing input data is thrown", {
  expect_error(
    gibbs_sampler(
      Z = cbind(toy_data$z1, toy_data$z2),
      y = toy_data$y,
      beta_start = beta, gamma_start = gamma,
      num_sim = 50
    ),
    "At least either all model matrices (X, Z, y) and coefficients (beta_start, gamma_start) or a model object (m) must be given.",
    fixed = TRUE
  )
})
