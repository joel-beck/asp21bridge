library(lslmbRidge)

fit <- lslm(
  location = y ~ x1 + x2 + z1 + z2, scale = ~ z1 + z2,
  data = toy_data, light = FALSE
) %>%
  gibbs_sampler(num_sim = 50)

samples <- fit$mcmc_ridge$coefficient_sampling_matrices

test_that("burnin works for numeric vectors", {
  burned <- drop(samples$location_prior) %>% burnin(num_burn = 10)

  expect_equal(length(burned), 40)
  expect_equal(drop(samples$location_prior)[11], burned[1])
  expect_type(burned, "double")
})

test_that("burnin does not convert single column matrices to vectors", {
  burned <- samples$location_prior %>% burnin(num_burn = 10)

  expect_equal(ncol(burned), 1)
  expect_false(is.null(colnames(burned)))
})

test_that("burnin works for matrices", {
  burned <- samples$location %>% burnin(num_burn = 10)

  expect_equal(nrow(burned), 40)
  expect_equal(samples$location[11,], burned[1,])
  expect_equal(ncol(burned), 5)
})

test_that("burnin keeps colnames", {
  burned_beta <- samples$location %>% burnin(num_burn = 10)
  burned_tau <- samples$location_prior %>% burnin(num_burn = 10)

  expect_equal(colnames(burned_beta), paste("beta", 0:4, sep = "_"))
  expect_equal(colnames(burned_tau), "tau^2")
})

test_that("burnin works for lists", {
  burned <- samples %>% burnin(num_burn = 10)

  expect_equal(length(burned), 4)
  expect_equal(nrow(burned[[1]]), 40)
  expect_equal(samples[[2]][11,], burned[[2]][1,])
  expect_type(burned, "list")
})

