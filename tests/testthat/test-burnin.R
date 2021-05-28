fit <- lslm(
  location = y ~ x1 + x2 + z1 + z2, scale = ~ z1 + z2,
  data = toy_data, light = FALSE
) %>%
  gibbs_sampler(num_sim = 50)

samples <- fit$mcmc_ridge$sampling_matrices

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

test_that("burnin handles edge cases correctly", {
  burned_1 <- samples$location %>% burnin(num_burn = 0)
  burned_2 <- samples$location %>% burnin(num_burn = 49)

  expect_equal(samples$location, burned_1)
  expect_equal(samples$location[50,, drop = FALSE], burned_2)
})

test_that("burnin gives helpful error messages", {
  expect_error(samples$location %>% burnin(num_burn = -1),
               regexp = "'num_burn' must be a nonnegative integer!")
  expect_error(samples$location %>% burnin(num_burn = 0.5),
               regexp = "'num_burn' must be a nonnegative integer!")
  expect_error(samples$location %>% burnin(num_burn = 50),
               regexp = "'num_burn' must be less than the number of samples!"
  )
})

