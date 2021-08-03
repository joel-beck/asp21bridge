fit <- lmls(
  location = y ~ x1 + x2 + z1 + z2, scale = ~ z1 + z2,
  data = toy_data, light = FALSE
) %>%
  mcmc_ridge(num_sim = 1000)

samples_list1 <- fit$mcmc_ridge
samples_list2 <- fit$mcmc_ridge$sampling_matrices
samples_matrix1 <- samples_list2$scale
samples_matrix2 <- samples_list2$location[, 2, drop = FALSE]
samples_vector <- samples_list2$location[, 2]


### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### Tests for valid input                                                   ####

test_that("model input works", {
  expect_equal(dim(summary_complete(samples = fit)), c(10, 6))
  expect_equal(
    summary_complete(samples = fit)$Parameter,
    c(
      "beta_0", "beta_1", "beta_2", "beta_3", "beta_4",
      "gamma_0", "gamma_1", "gamma_2", "tau^2", "xi^2"
    )
  )
  expect_type(summary_complete(samples = fit)[[1]], "character")
  expect_type(summary_complete(samples = fit)[[2]], "double")
})

test_that("list input works", {
  expect_equal(dim(summary_complete(samples = samples_list1)), c(10, 6))
  expect_equal(
    summary_complete(samples = samples_list1)$Parameter,
    c(
      "beta_0", "beta_1", "beta_2", "beta_3", "beta_4",
      "gamma_0", "gamma_1", "gamma_2", "tau^2", "xi^2"
    )
  )
  expect_type(summary_complete(samples = samples_list1)[[1]], "character")
  expect_type(summary_complete(samples = samples_list1)[[2]], "double")
  expect_equal(dim(summary_complete(samples = samples_list2)), c(10, 6))
  expect_equal(
    summary_complete(samples = samples_list2)$Parameter,
    c(
      "beta_0", "beta_1", "beta_2", "beta_3", "beta_4",
      "gamma_0", "gamma_1", "gamma_2", "tau^2", "xi^2"
    )
  )
  expect_type(summary_complete(samples = samples_list2)[[1]], "character")
  expect_type(summary_complete(samples = samples_list2)[[2]], "double")
})

test_that("matrix input works", {
  expect_equal(dim(summary_complete(samples = samples_matrix1)), c(3, 6))
  expect_equal(
    summary_complete(samples = samples_matrix1)$Parameter,
    c("gamma_0", "gamma_1", "gamma_2")
  )
  expect_type(summary_complete(samples = samples_matrix1)[[1]], "character")
  expect_type(summary_complete(samples = samples_matrix1)[[2]], "double")
  expect_equal(dim(summary_complete(samples = samples_matrix2)), c(1, 6))
  expect_equal(summary_complete(samples = samples_matrix2)$Parameter, "beta_1")
  expect_type(summary_complete(samples = samples_matrix2)[[1]], "character")
  expect_type(summary_complete(samples = samples_matrix2)[[2]], "double")
})

test_that("vector input works", {
  expect_equal(dim(summary_complete(samples = samples_vector)), c(1, 5))
  expect_false(
    any(names(summary_complete(samples = samples_vector)) == "Parameter")
  )
  expect_type(summary_complete(samples = samples_vector)[[1]], "double")
})


### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### Tests for include_plot argument                                         ####

test_that("include_plot = TRUE works", {
  expect_equal(
    dim(summary_complete(samples = fit, include_plot = TRUE)), c(10, 7)
  )
  expect_equal(
    names(summary_complete(samples = fit, include_plot = TRUE))[7], "Plot"
  )
  expect_type(
    summary_complete(samples = fit, include_plot = TRUE)$Plot, "list"
  )
})
