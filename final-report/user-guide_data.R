library(asp21bridge)
set.seed(1234)

toy_fit <- lmls(
  location = y ~ x1 + x2 + z1 + z2, scale = ~ z1 + z2,
  data = toy_data, light = FALSE
) %>%
  mcmc_ridge(num_sim = 10000)

standardize <- function(x) (x - mean(x)) / sd(x)
y <- abdom$y
X <- as.matrix(standardize(abdom$x))
Z <- X

abdom_fit <- mcmc_ridge(y = y, X = X, Z = Z, beta_start = 1, gamma_start = 1, num_sim = 10000)

toy_trace <- mult_plot(samples = toy_fit, type = "trace", free_scale = TRUE, latex = TRUE)

variance_samples <- cbind(
  abdom_fit$sampling_matrices$tau_samples,
  abdom_fit$sampling_matrices$xi_samples
)

abdom_mult <- mult_plot(
  samples = variance_samples, type = "both", free_scale = TRUE,
  log = TRUE, latex = TRUE
)

abdom_diag_raw <- diagnostic_plots(
  samples = abdom_fit$sampling_matrices$beta_samples[, "beta_1", drop = FALSE],
  lag_max = 100, latex = TRUE
)

abdom_diag_thinned <- abdom_fit$sampling_matrices$beta_samples[, "beta_1", drop = FALSE] %>%
  burnin(num_burn = 2500) %>%
  thinning(freq = 10) %>%
  diagnostic_plots(lag_max = 100, latex = TRUE)

readr::write_rds(
  x = list(
    toy_fit = toy_fit, abdom_fit = abdom_fit, toy_trace = toy_trace,
    abdom_mult = abdom_mult, abdom_diag_raw = abdom_diag_raw,
    abdom_diag_thinned = abdom_diag_thinned
  ),
  file = here::here("final-report", "user-guide_data.rds")
)
