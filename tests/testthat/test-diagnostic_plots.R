#   ____________________________________________________________________________
#   Unit Tests for Single Chain Plotting Functions                          ####

fit <- lmls(
  location = y ~ x1 + x2 + z1 + z2, scale = ~ z1 + z2,
  data = toy_data, light = FALSE
) %>%
  mcmc_ridge(num_sim = 1000)

samples <- fit$mcmc_ridge$sampling_matrices

samples_list <- list(samples$scale_prior)
samples_matrix <- samples$location_prior
samples_vector <- samples$location[, 2]


### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### No Errors                                                               ####

test_that("no errors for list input", {
  expect_error(time_plot(samples = samples_list), NA)
  expect_error(density_plot(samples = samples_list), NA)
  expect_error(acl_plot(samples = samples_list), NA)
  expect_error(diagnostic_plots(samples = samples_list), NA)
})

test_that("no errors for matrix input", {
  expect_error(time_plot(samples = samples_matrix), NA)
  expect_error(density_plot(samples = samples_matrix), NA)
  expect_error(acl_plot(samples = samples_matrix), NA)
  expect_error(diagnostic_plots(samples = samples_matrix), NA)
})

test_that("no errors for vector input", {
  expect_error(time_plot(samples = samples_vector), NA)
  expect_error(density_plot(samples = samples_vector), NA)
  expect_error(acl_plot(samples = samples_vector), NA)
  expect_error(diagnostic_plots(samples = samples_vector), NA)
})

test_that("no errors for list input with arguments", {
  expect_error(time_plot(
    samples = samples_list, log = TRUE, robust = TRUE, latex = TRUE
  ), NA)
  expect_error(density_plot(
    samples = samples_list, log = TRUE, robust = TRUE, latex = TRUE
  ), NA)
  expect_error(acl_plot(
    samples = samples_list, lag_max = 100, latex = TRUE
  ), NA)
  expect_error(diagnostic_plots(
    samples = samples_list, lag_max = 100, log = TRUE, robust = TRUE, latex = TRUE
  ), NA)
})

test_that("no errors for matrix input with arguments", {
  expect_error(time_plot(
    samples = samples_matrix, log = TRUE, robust = TRUE, latex = TRUE
  ), NA)
  expect_error(density_plot(
    samples = samples_matrix, log = TRUE, robust = TRUE, latex = TRUE
  ), NA)
  expect_error(acl_plot(
    samples = samples_matrix, lag_max = 100, latex = TRUE
  ), NA)
  expect_error(diagnostic_plots(
    samples = samples_matrix, lag_max = 100, log = TRUE, robust = TRUE, latex = TRUE
  ), NA)
})

test_that("no errors for vector input with arguments", {
  expect_error(time_plot(
    samples = samples_vector, log = TRUE, robust = TRUE
  ), NA)
  expect_error(density_plot(
    samples = samples_vector, log = TRUE, robust = TRUE
  ), NA)
  expect_error(acl_plot(
    samples = samples_vector, lag_max = 100
  ), NA)
  expect_error(diagnostic_plots(
    samples = samples_vector, lag_max = 100, log = TRUE, robust = TRUE
  ), NA)
})



### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### No Warnings                                                             ####

test_that("no warnings for list input", {
  expect_warning(time_plot(samples = samples_list), NA)
  expect_warning(density_plot(samples = samples_list), NA)
  expect_warning(acl_plot(samples = samples_list), NA)
  expect_warning(diagnostic_plots(samples = samples_list), NA)
})

test_that("no warnings for matrix input", {
  expect_warning(time_plot(samples = samples_matrix), NA)
  expect_warning(density_plot(samples = samples_matrix), NA)
  expect_warning(acl_plot(samples = samples_matrix), NA)
  expect_warning(diagnostic_plots(samples = samples_matrix), NA)
})

test_that("no warnings for vector input", {
  expect_warning(time_plot(samples = samples_vector), NA)
  expect_warning(density_plot(samples = samples_vector), NA)
  expect_warning(acl_plot(samples = samples_vector), NA)
  expect_warning(diagnostic_plots(samples = samples_vector), NA)
})

test_that("no warnings for list input with arguments", {
  expect_warning(time_plot(
    samples = samples_list, log = TRUE, robust = TRUE, latex = TRUE
  ), NA)
  expect_warning(density_plot(
    samples = samples_list, log = TRUE, robust = TRUE, latex = TRUE
  ), NA)
  expect_warning(acl_plot(
    samples = samples_list, lag_max = 100, latex = TRUE
  ), NA)
  expect_warning(diagnostic_plots(
    samples = samples_list, lag_max = 100, log = TRUE, robust = TRUE, latex = TRUE
  ), NA)
})

test_that("no warnings for matrix input with arguments", {
  expect_warning(time_plot(
    samples = samples_matrix, log = TRUE, robust = TRUE, latex = TRUE
  ), NA)
  expect_warning(density_plot(
    samples = samples_matrix, log = TRUE, robust = TRUE, latex = TRUE
  ), NA)
  expect_warning(acl_plot(
    samples = samples_matrix, lag_max = 100, latex = TRUE
  ), NA)
  expect_warning(diagnostic_plots(
    samples = samples_matrix, lag_max = 100, log = TRUE, robust = TRUE, latex = TRUE
  ), NA)
})

test_that("no warnings for vector input with arguments", {
  expect_warning(time_plot(
    samples = samples_vector, log = TRUE, robust = TRUE
  ), NA)
  expect_warning(density_plot(
    samples = samples_vector, log = TRUE, robust = TRUE
  ), NA)
  expect_warning(acl_plot(
    samples = samples_vector, lag_max = 100
  ), NA)
  expect_warning(diagnostic_plots(
    samples = samples_vector, lag_max = 100, log = TRUE, robust = TRUE
  ), NA)
})



### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### Specific Error Messages                                                 ####

input_message <- paste(
  "Input 'samples' must be a numeric vector, a matrix with one column",
  "or a list with a single element of these types!"
)

test_that("correct error message for invalid input", {
  expect_error(time_plot(samples = fit$mcmc_ridge), input_message)
  expect_error(density_plot(samples = samples), input_message)
  expect_error(acl_plot(samples = samples$location), input_message)
  expect_error(diagnostic_plots(samples = samples$location[, 2:3], input_message))
})


unnamed_matrix <- samples$location_prior
colnames(unnamed_matrix) <- NULL
latex_message <- "Set 'latex = TRUE' only if column names are provided."

test_that("latex = FALSE for missing column names", {
  expect_error(time_plot(
    samples = samples_vector, latex = TRUE
  ), latex_message)
  expect_error(density_plot(
    samples = list(samples_vector), latex = TRUE
  ), latex_message)
  expect_error(acl_plot(
    samples = unnamed_matrix, latex = TRUE
  ), latex_message)
  expect_error(diagnostic_plots(
    samples = list(unnamed_matrix), latex = TRUE
  ), latex_message)
})



### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### Plot Title                                                              ####

t_plot <- time_plot(samples = samples_matrix, latex = TRUE)
d_plot <- density_plot(samples = samples_list, latex = TRUE)
a_plot <- acl_plot(samples = samples$location[, 2, drop = FALSE], latex = TRUE)
diag_plot <- diagnostic_plots(samples = samples_matrix, latex = TRUE)

test_that("latex symbols are included in plot title", {
  expect_type(t_plot$labels$title, type = "expression")
  expect_type(d_plot$labels$title, type = "expression")
  expect_type(a_plot$labels$title, type = "expression")
})
