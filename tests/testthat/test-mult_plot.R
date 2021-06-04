fit <- lslm(
  location = y ~ x1 + x2 + z1 + z2, scale = ~ z1 + z2,
  data = toy_data, light = FALSE
) %>%
  gibbs_sampler(num_sim = 1000)

samples_list1 <- fit$mcmc_ridge
samples_list2 <- fit$mcmc_ridge$sampling_matrices
samples_matrix1 <- samples_list2$location
samples_matrix2 <- samples_list2$scale_prior

### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### Tests for valid input                                                   ####

test_that("input is model object", {
  expect_error(mult_plot(samples = fit), NA)
  expect_warning(mult_plot(samples = fit), NA)
  expect_message(mult_plot(samples = fit), NA)
})

test_that("input is model object with arguments", {
  expect_error(mult_plot(
    samples = fit, type = "density", robust = TRUE, latex = TRUE
  ), NA)
  expect_warning(mult_plot(
    samples = fit, type = "time", free_scale = TRUE, latex = TRUE
  ), NA)
  expect_message(mult_plot(
    samples = fit, type = "both", robust = TRUE, free_scale = TRUE,
    latex = TRUE
  ), NA)
})

test_that("input is list, but not a model object", {
  expect_error(mult_plot(samples = samples_list1), NA)
  expect_warning(mult_plot(samples = samples_list1), NA)
  expect_message(mult_plot(samples = samples_list1), NA)
  expect_error(mult_plot(samples = samples_list2), NA)
  expect_warning(mult_plot(samples = samples_list2), NA)
  expect_message(mult_plot(samples = samples_list2), NA)
})

test_that("input is list with arguments, but not a model object", {
  expect_error(mult_plot(
    samples = samples_list1, type = "density", robust = TRUE, latex = TRUE
  ), NA)
  expect_warning(mult_plot(
    samples = samples_list1, type = "time", free_scale = TRUE, latex = TRUE
  ), NA)
  expect_message(mult_plot(
    samples = samples_list1, type = "both", robust = TRUE, free_scale = TRUE,
    latex = TRUE
  ), NA)
  expect_error(mult_plot(
    samples = samples_list2, type = "density", robust = TRUE, latex = TRUE
  ), NA)
  expect_warning(mult_plot(
    samples = samples_list2, type = "time", free_scale = TRUE, latex = TRUE
  ), NA)
  expect_message(mult_plot(
    samples = samples_list2, type = "both", robust = TRUE, free_scale = TRUE,
    latex = TRUE
  ), NA)
})

test_that("input is matrix", {
  expect_error(mult_plot(samples = samples_matrix1), NA)
  expect_warning(mult_plot(samples = samples_matrix2), NA)
  expect_message(mult_plot(samples = samples_matrix1), NA)
})

test_that("input is matrix with arguments", {
  expect_error(mult_plot(
    samples = samples_matrix1, type = "density", robust = TRUE, latex = TRUE
  ), NA)
  expect_warning(mult_plot(
    samples = samples_matrix1, type = "time", free_scale = TRUE, latex = TRUE
  ), NA)
  expect_message(mult_plot(
    samples = samples_matrix2, type = "both", log = TRUE, robust = TRUE,
    free_scale = TRUE, latex = TRUE
  ), NA)
})


### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### Tests for invalid input                                                 ####

input_message <- paste(
  "Input 'samples' must be a model object of the class 'lslm',",
  "a numeric matrix or a list of vectors and matrices!"
)
logical_message <- paste(
  "Inputs 'log', 'robust', 'free_scale' and 'latex'",
  "must be either 'TRUE' or 'FALSE'!"
)

test_that("correct error messages are produced", {
  expect_error(mult_plot(samples = samples_matrix1[, 1]), input_message)
  expect_error(mult_plot(samples = samples_matrix2, robust = 1), logical_message)
  expect_error(mult_plot(samples = samples_matrix2, latex = 0), logical_message)
})


### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### Tests for latex labels                                                  ####

p1 <- mult_plot(samples = samples_list1, type = "time", latex = TRUE)

p2 <- mult_plot(samples = samples_list2, type = "density", latex = TRUE)

p3 <- mult_plot(samples = samples_matrix1, type = "both", latex = TRUE)

p4 <- mult_plot(samples = samples_matrix2, type = "both", log = TRUE, latex = TRUE)

test_that("latex labels are recognized and transformed", {
  expect_equal(
    levels(p1$data$Parameter), c(
      "beta[0]", "beta[1]", "beta[2]", "beta[3]",
      "beta[4]", "gamma[0]", "gamma[1]", "gamma[2]",
      "tau^2", "xi^2"
    )
  )
  expect_equal(
    levels(p2$data$Parameter), c(
      "beta[0]", "beta[1]", "beta[2]", "beta[3]",
      "beta[4]", "gamma[0]", "gamma[1]", "gamma[2]",
      "tau^2", "xi^2"
    )
  )
  expect_equal(
    levels(p3$data$Parameter), c(
      "beta[0]", "beta[1]", "beta[2]", "beta[3]", "beta[4]"
    )
  )
  expect_equal(levels(p4$data$Parameter), "xi^2")
})
