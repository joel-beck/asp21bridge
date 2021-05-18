library(lslmbRidge)

fit <- lslm(
  location = y ~ x1 + x2 + z1 + z2, scale = ~ z1 + z2,
  data = toy_data, light = FALSE
) %>%
  gibbs_sampler(num_sim = 50)

samples <- fit$mcmc_ridge$coefficient_sampling_matrices

test_that("thinning works for numeric vectors", {
  thinned <- drop(samples$location_prior) %>% thinning(freq = 5)

  expect_equal(length(thinned), 10)
  expect_equal(drop(samples$location_prior)[11], thinned[3])
  expect_type(thinned, "double")
})

test_that("thinning does not convert single column matrices to vectors", {
  thinned <- samples$location_prior %>% thinning(freq = 5)

  expect_equal(ncol(thinned), 1)
  expect_false(is.null(colnames(thinned)))
})

test_that("thinning works for matrices", {
  thinned <- samples$location %>% thinning(freq = 5)

  expect_equal(nrow(thinned), 10)
  expect_equal(samples$location[11, ], thinned[3, ])
  expect_equal(ncol(thinned), 5)
})

test_that("thinning keeps colnames", {
  thinned_beta <- samples$location %>% thinning(freq = 5)
  thinned_tau <- samples$location_prior %>% thinning(freq = 5)

  expect_equal(colnames(thinned_beta), paste("beta", 0:4, sep = "_"))
  expect_equal(colnames(thinned_tau), "tau^2")
})

test_that("thinning works for lists", {
  thinned <- samples %>% thinning(freq = 5)

  expect_equal(length(thinned), 4)
  expect_equal(nrow(thinned[[1]]), 10)
  expect_equal(samples[[2]][11, ], thinned[[2]][3, ])
  expect_type(thinned, "list")
})

test_that("thinning handles edge cases correctly", {
  thinned_1 <- samples$location %>% thinning(freq = 1)
  thinned_2 <- samples$location %>% thinning(freq = 50)

  expect_equal(samples$location, thinned_1)
  expect_equal(samples$location[1,, drop = FALSE], thinned_2)
})

test_that("thinning gives helpful error messages", {
  expect_error(samples$location %>% thinning(freq = -1),
               regexp = "'freq' must be a positive integer!")
  expect_error(samples$location %>% thinning(freq = 0.5),
               regexp = "'freq' must be a positive integer!")
  expect_error(samples$location %>% thinning(freq = 51),
    regexp = "'freq' must not be larger than the number of samples!"
  )
})
