### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### Setup                                                                   ####

pacman::p_load(
  dplyr,
  purrr,
  ggplot2,
  gridExtra,
  scales
)
library(asp21bridge)


### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### Simulate Data                                                           ####
set.seed(142)
x1 <- rnorm(n = 50, mean = 1, sd = 1)
x2 <- rnorm(n = 50, mean = 2, sd = 1)
z1 <- rnorm(n = 50, mean = 5, sd = 1)
z2 <- rnorm(n = 50, mean = 3, sd = 1)

X <- cbind(x1, x2)
Z <- cbind(z1, z2)

beta <- c(-1, 4) # True: beta_0 = 0, beta_1 = -1, beta_2 = 4,
gamma <- c(-2, 1) # True: gamma_0 = 0, gamma_1 = -2, gamma_2 = 1

y <- vector(mode = "numeric", length = 50)
for (i in seq_along(y)) {
  mu <- sum(X[i,] * beta)
  sigma <- exp(sum(Z[i,] * gamma))
  y[i] <- rnorm(n = 1, mean = mu, sd = sigma)
}

# Values for hyper parameters
hyppar_val <- c(-1, 0, 0.5, 1, 2, 10, 50, 100, 200)


### Create lmls model
model <- lmls(
  location = y ~ x1 + x2, scale = ~ z1 + z2,
  light = FALSE
)

### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### Create initial Tibbles                                                  ####
# One tibble per hyperpar-vector
# Default for other hyper-pars = 1
# Choosed values to have specific ratios

a_tau_data <- tibble(hyppar_val = hyppar_val) %>%
  mutate(hyppar = "a_tau") %>%
  mutate(samples = map(
    .x = hyppar_val, .f = ~ mcmc_ridge(m = model, a_tau = .x, num_sim = 100)
  )) %>%
  mutate(acc_rate = map_dbl(
    .x = samples, .f = ~ .x$mcmc_ridge$acceptance_rate
  )) %>%
  mutate(results = map(
    .x = samples,
    .f = ~ summary_complete(.x) %>% select(Parameter, `Posterior Mean`)
  ))


b_tau_data <- tibble(hyppar_val = hyppar_val) %>%
  mutate(hyppar = "b_tau") %>%
  mutate(samples = map(
    .x = hyppar_val, .f = ~ mcmc_ridge(m = model, b_tau = .x, num_sim = 100)
  )) %>%
  mutate(acc_rate = map_dbl(
    .x = samples, .f = ~ .x$mcmc_ridge$acceptance_rate
  )) %>%
  mutate(results = map(
    .x = samples,
    .f = ~ summary_complete(.x) %>% select(Parameter, `Posterior Mean`)
  ))

a_xi_data <- tibble(hyppar_val = hyppar_val) %>%
  mutate(hyppar = "a_xi") %>%
  mutate(samples = map(
    .x = hyppar_val, .f = ~ mcmc_ridge(m = model, a_xi = .x, num_sim = 100)
  )) %>%
  mutate(acc_rate = map_dbl(
    .x = samples, .f = ~ .x$mcmc_ridge$acceptance_rate
  )) %>%
  mutate(results = map(
    .x = samples,
    .f = ~ summary_complete(.x) %>% select(Parameter, `Posterior Mean`)
  ))

b_xi_data <- tibble(hyppar_val = hyppar_val) %>%
  mutate(hyppar = "b_xi") %>%
  mutate(samples = map(
    .x = hyppar_val, .f = ~ mcmc_ridge(m = model, b_xi = .x, num_sim = 100)
  )) %>%
  mutate(acc_rate = map_dbl(
    .x = samples, .f = ~ .x$mcmc_ridge$acceptance_rate
  )) %>%
  mutate(results = map(
    .x = samples,
    .f = ~ summary_complete(.x) %>% select(Parameter, `Posterior Mean`)
  ))


### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### Analysis                                                                ####

# Creating results tables
results_a_tau_data <- a_tau_data %>%
  select(hyppar, hyppar_val, results) %>%
  tidyr::unnest(results)

results_b_tau_data <- b_tau_data %>%
  select(hyppar, hyppar_val, results) %>%
  tidyr::unnest(results)

results_a_xi_data <- a_xi_data %>%
  select(hyppar, hyppar_val, results) %>%
  tidyr::unnest(results)

results_b_xi_data <- b_xi_data %>%
  select(hyppar, hyppar_val, results) %>%
  tidyr::unnest(results)

# True value of -1, no variance of posterior means
results_a_tau_data %>%
  filter(Parameter == "beta_1") %>%
  print(n = Inf)

results_b_tau_data %>%
  filter(Parameter == "beta_1") %>%
  print(n = Inf)

results_b_xi_data %>%
  filter(Parameter == "beta_1") %>%
  print(n = Inf)

# True value of -2, small variance but steady small overestimation
results_a_tau_data %>%
  filter(Parameter == "gamma_1") %>%
  print(n = Inf)

results_b_xi_data %>%
  filter(Parameter == "gamma_1") %>%
  print(n = Inf)

results_b_xi_data %>%
  filter(Parameter == "gamma_1") %>%
  print(n = Inf)

# Adding absolute deviation (estimated - true value)
deviation_a_tau_data <- results_a_tau_data %>%
  filter(stringr::str_detect(string = Parameter, pattern = "beta|gamma")) %>%
  mutate(truth = rep(c(0, -1, 4, 0, -2, 1), times = n_distinct(hyppar_val))) %>%
  mutate(deviation = abs(`Posterior Mean` - truth))

deviation_b_tau_data <- results_b_tau_data %>%
  filter(stringr::str_detect(string = Parameter, pattern = "beta|gamma")) %>%
  mutate(truth = rep(c(0, -1, 4, 0, -2, 1), times = n_distinct(hyppar_val))) %>%
  mutate(deviation = abs(`Posterior Mean` - truth))

deviation_a_xi_data <- results_a_xi_data %>%
  filter(stringr::str_detect(string = Parameter, pattern = "beta|gamma")) %>%
  mutate(truth = rep(c(0, -1, 4, 0, -2, 1), times = n_distinct(hyppar_val))) %>%
  mutate(deviation = abs(`Posterior Mean` - truth))

deviation_b_xi_data <- results_b_xi_data %>%
  filter(stringr::str_detect(string = Parameter, pattern = "beta|gamma")) %>%
  mutate(truth = rep(c(0, -1, 4, 0, -2, 1), times = n_distinct(hyppar_val))) %>%
  mutate(deviation = abs(`Posterior Mean` - truth))

# two biggest and two smallest deviations for each parameter and a_tau / b_xi
# Evenly spread values of a_tau, only 100 is always at small, never at biggest
# No pattern of distinct hyperpar values obtainable at the betas
# At the gammas, values 0 < a_tau <= 1 are always at smallest
# The gammas have a higher deviation, the betas slightly none
bind_rows(
  biggest = deviation_a_tau_data %>%
    group_by(Parameter) %>%
    slice_max(order_by = deviation, n = 2) %>%
    ungroup(),
  smallest = deviation_a_tau_data %>%
    group_by(Parameter) %>%
    slice_min(order_by = deviation, n = 2) %>%
    ungroup(),
  .id = "diff"
) %>%
  select(hyppar, Parameter, diff, hyppar_val, deviation) %>%
  arrange(Parameter) %>%
  print(n = Inf)

# for beta, the values of b_xi are evenly spread
# for gamma, always b_xi = 100 are at smallest
# The gammas have a higher deviation, the betas slightly none
bind_rows(
  biggest = deviation_b_xi_data %>%
    group_by(Parameter) %>%
    slice_max(order_by = deviation, n = 2) %>%
    ungroup(),
  smallest = deviation_b_xi_data %>%
    group_by(Parameter) %>%
    slice_min(order_by = deviation, n = 2) %>%
    ungroup(),
  .id = "diff"
) %>%
  select(hyppar, Parameter, diff, hyppar_val, deviation) %>%
  arrange(Parameter) %>%
  print(n = Inf)

# how often is each hyper parameter among the four biggest or four smallest
# deviations from the true value across all parameters and a_tau / b _xi?
# Here, I moreover SPLIT b/w beta and gamma as parameters.

# small values have smaller deviation more often, while large values have larger deviation
# no remarkable differences b/w beta and gamma

bind_rows(
  biggest_beta = deviation_a_tau_data %>%
    filter(stringr::str_detect(string = Parameter, pattern = "beta")) %>%
    group_by(Parameter) %>%
    slice_max(order_by = deviation, n = 4) %>%
    ungroup(),
  smallest_beta = deviation_a_tau_data %>%
    filter(stringr::str_detect(string = Parameter, pattern = "beta")) %>%
    group_by(Parameter) %>%
    slice_min(order_by = deviation, n = 4) %>%
    ungroup(),
  biggest_gamma = deviation_a_tau_data %>%
    filter(stringr::str_detect(string = Parameter, pattern = "gamma")) %>%
    group_by(Parameter) %>%
    slice_max(order_by = deviation, n = 4) %>%
    ungroup(),
  smallest_gamma = deviation_a_tau_data %>%
    filter(stringr::str_detect(string = Parameter, pattern = "gamma")) %>%
    group_by(Parameter) %>%
    slice_min(order_by = deviation, n = 4) %>%
    ungroup(),
  .id = "diff"
) %>%
  select(hyppar, diff, hyppar_val, Parameter, deviation) %>%
  count(diff, hyppar_val) %>%
  tidyr::pivot_wider(names_from = diff, values_from = n) %>%
  arrange(hyppar_val) %>%
  print(n = Inf)


# beta's and gamma's deviation both is large for small values of hyperpar
# and small for large values
bind_rows(
  biggest_beta = deviation_b_xi_data %>%
    filter(stringr::str_detect(string = Parameter, pattern = "beta")) %>%
    group_by(Parameter) %>%
    slice_max(order_by = deviation, n = 4) %>%
    ungroup(),
  smallest_beta = deviation_b_xi_data %>%
    filter(stringr::str_detect(string = Parameter, pattern = "beta")) %>%
    group_by(Parameter) %>%
    slice_min(order_by = deviation, n = 4) %>%
    ungroup(),
  biggest_gamma = deviation_b_xi_data %>%
    filter(stringr::str_detect(string = Parameter, pattern = "gamma")) %>%
    group_by(Parameter) %>%
    slice_max(order_by = deviation, n = 4) %>%
    ungroup(),
  smallest_gamma = deviation_b_xi_data %>%
    filter(stringr::str_detect(string = Parameter, pattern = "gamma")) %>%
    group_by(Parameter) %>%
    slice_min(order_by = deviation, n = 4) %>%
    ungroup(),
  .id = "diff"
) %>%
  select(hyppar, diff, hyppar_val, Parameter, deviation) %>%
  count(diff, hyppar_val) %>%
  tidyr::pivot_wider(names_from = diff, values_from = n) %>%
  arrange(hyppar_val) %>%
  print(n = Inf)



# plot deviations for beta vs gamma Parameter for a_tau / b_xi
# largest influence in b_tau / gamma plot
p1 <- deviation_a_tau_data %>%
  filter(stringr::str_detect(string = Parameter, pattern = "beta")) %>%
  ggplot(mapping = aes(x = hyppar_val, y = deviation, color = Parameter)) +
  geom_line() +
  geom_point(size = 2) +
  labs(
    title = "Abs. dev. from true beta (#sim = 1000)",
    x = "Values for a_tau", y = "Absolute Deviation"
  ) +
  theme_light() +
  theme(panel.grid.minor = element_blank()) +
  geom_smooth(method = "lm", aes(colour = 'linear trend'), linetype = "dashed",
              se = FALSE) +
  scale_x_continuous(trans = pseudo_log_trans(), breaks = c(-1, 0, 1, 2, 10, 50, 100, 200))

p2 <- deviation_a_tau_data %>%
  filter(stringr::str_detect(string = Parameter, pattern = "gamma")) %>%
  ggplot(mapping = aes(x = hyppar_val, y = deviation, color = Parameter)) +
  geom_line() +
  geom_point(size = 2) +
  labs(
    title = "Abs. dev. from true gamma (#sim = 1000)",
    x = "Values for a_tau", y = "Absolute Deviation"
  ) +
  theme_light() +
  theme(panel.grid.minor = element_blank()) +
  geom_smooth(method = "lm", aes(colour = 'linear trend'), linetype = "dashed",
              se = FALSE) +
  scale_x_continuous(trans = pseudo_log_trans(), breaks = c(-1, 0, 1, 2, 10, 50, 100, 200))

p3 <- deviation_b_tau_data %>%
  filter(stringr::str_detect(string = Parameter, pattern = "beta")) %>%
  ggplot(mapping = aes(x = hyppar_val, y = deviation, color = Parameter)) +
  geom_line() +
  geom_point(size = 2) +
  labs(
    title = "Abs. dev. from true beta (#sim = 1000)",
    x = "Values for b_tau", y = "Absolute Deviation"
  ) +
  theme_light() +
  theme(panel.grid.minor = element_blank()) +
  geom_smooth(method = "lm", aes(colour = 'linear trend'), linetype = "dashed",
              se = FALSE) +
  scale_x_continuous(trans = pseudo_log_trans(), breaks = c(-1, 0, 1, 2, 10, 50, 100, 200))

p4 <- deviation_b_tau_data %>%
  filter(stringr::str_detect(string = Parameter, pattern = "gamma")) %>%
  ggplot(mapping = aes(x = hyppar_val, y = deviation, color = Parameter)) +
  geom_line() +
  geom_point(size = 2) +
  labs(
    title = "Abs. dev. from true gamma (#sim = 1000)",
    x = "Values for b_tau", y = "Absolute Deviation"
  ) +
  theme_light() +
  theme(panel.grid.minor = element_blank()) +
  geom_smooth(method = "lm", aes(colour = 'linear trend'), linetype = "dashed",
              se = FALSE) +
  scale_x_continuous(trans = pseudo_log_trans(), breaks = c(-1, 0, 1, 2, 10, 50, 100, 200))

(tau_dev_plot <- grid.arrange(p1, p3, p2, p4, nrow = 2))

# plot deviations for beta vs gamma Parameter for hypperparams of xi
# largest influence in b_xi / gamma plot
p5 <- deviation_a_xi_data %>%
  filter(stringr::str_detect(string = Parameter, pattern = "beta")) %>%
  ggplot(mapping = aes(x = hyppar_val, y = deviation, color = Parameter)) +
  geom_line() +
  geom_point(size = 2) +
  labs(
    title = "Abs. dev. from true beta (#sim = 1000)",
    x = "Values for a_xi", y = "Absolute Deviation"
  ) +
  theme_light() +
  theme(panel.grid.minor = element_blank()) +
  geom_smooth(method = "lm", aes(colour = 'linear trend'), linetype = "dashed",
              se = FALSE) +
  scale_x_continuous(trans = pseudo_log_trans(), breaks = c(-1, 0, 1, 2, 10, 50, 100, 200))

p6 <- deviation_a_xi_data %>%
  filter(stringr::str_detect(string = Parameter, pattern = "gamma")) %>%
  ggplot(mapping = aes(x = hyppar_val, y = deviation, color = Parameter)) +
  geom_line() +
  geom_point(size = 2) +
  labs(
    title = "Abs. dev. from true gamma (#sim = 1000)",
    x = "Values for a_xi", y = "Absolute Deviation"
  ) +
  theme_light() +
  theme(panel.grid.minor = element_blank()) +
  geom_smooth(method = "lm", aes(colour = 'linear trend'), linetype = "dashed",
              se = FALSE) +
  scale_x_continuous(trans = pseudo_log_trans(), breaks = c(-1, 0, 1, 2, 10, 50, 100, 200))

p7 <- deviation_b_xi_data %>%
  filter(stringr::str_detect(string = Parameter, pattern = "beta")) %>%
  ggplot(mapping = aes(x = hyppar_val, y = deviation, color = Parameter)) +
  geom_line() +
  geom_point(size = 2) +
  labs(
    title = "Abs. dev. from true beta (#sim = 1000)",
    x = "Values for b_xi", y = "Absolute Deviation"
  ) +
  theme_light() +
  theme(panel.grid.minor = element_blank()) +
  geom_smooth(method = "lm", aes(colour = 'linear trend'), linetype = "dashed",
              se = FALSE) +
  scale_x_continuous(trans = pseudo_log_trans(), breaks = c(-1, 0, 1, 2, 10, 50, 100, 200))

p8 <- deviation_b_xi_data %>%
  filter(stringr::str_detect(string = Parameter, pattern = "gamma")) %>%
  ggplot(mapping = aes(x = hyppar_val, y = deviation, color = Parameter)) +
  geom_line() +
  geom_point(size = 2) +
  labs(
    title = "Abs. dev. from true gamma (#sim = 1000)",
    x = "Values for b_xi", y = "Absolute Deviation"
  ) +
  theme_light() +
  theme(panel.grid.minor = element_blank()) +
  geom_smooth(method = "lm", aes(colour = 'linear trend'), linetype = "dashed",
              se = FALSE) +
  scale_x_continuous(trans = pseudo_log_trans(), breaks = c(-1, 0, 1, 2, 10, 50, 100, 200))

(xi_dev_plot <- grid.arrange(p5, p7, p6, p8, nrow = 2))

# highest and lowest acceptance rates for a_tau
a_tau_data %>%
  select(hyppar, hyppar_val, acc_rate) %>%
  slice_max(order_by = acc_rate, n = 3)

a_tau_data %>%
  select(hyppar, hyppar_val, acc_rate) %>%
  slice_min(order_by = acc_rate, n = 3)

# highest and lowest acceptance rates for b_xi
b_xi_data %>%
  select(hyppar, hyppar_val, acc_rate) %>%
  slice_max(order_by = acc_rate, n = 3)

b_xi_data %>%
  select(hyppar, hyppar_val, acc_rate) %>%
  slice_min(order_by = acc_rate, n = 3)

# correlation of hyperparams over covariates

plot_cor <- function(covariate) {
  # Plots hyperparameters' correlations against each others
  covariate <- as.character(covariate)

  deviation_binded <- tibble(a_tau = deviation_a_tau_data %>%
                               filter(stringr::str_detect(string = Parameter, pattern = covariate)) %>%
                               select(deviation)) %>%
    mutate(b_tau = deviation_b_tau_data %>%
             filter(stringr::str_detect(string = Parameter, pattern = covariate)) %>%
             select(deviation)) %>%
    mutate(a_xi = deviation_a_xi_data %>%
             filter(stringr::str_detect(string = Parameter, pattern = covariate)) %>%
             select(deviation)) %>%
    mutate(b_xi = deviation_b_xi_data %>%
             filter(stringr::str_detect(string = Parameter, pattern = covariate)) %>%
             select(deviation))

  panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...) # helper function from pairs() help-page
  {
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- cor(x, y)
    txt <- format(c(r, 0.123456789), digits = digits)[1]
    txt <- paste0(prefix, txt)
    if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = 2)
  }

  return(pairs(as.matrix(deviation_binded),
               upper.panel = panel.cor,
               lower.panel = panel.smooth,
               main = paste("Correlation plot of hyperparameters for", covariate)))
}
(beta0_cor_plot <- plot_cor("beta_0"))
(gamma0_cor_plot <- plot_cor("gamma_0"))

#   ____________________________________________________________________________
#   Save data for Second Report                                             ####

readr::write_rds(
  x = list(
    a_tau_data = a_tau_data,
    b_tau_data = b_tau_data,
    a_xi_data = a_xi_data,
    b_xi_data = b_xi_data,
    results_a_tau_data = results_a_tau_data,
    results_b_tau_data = results_b_tau_data,
    results_a_xi_data = results_a_xi_data,
    results_b_xi_data = results_b_xi_data,
    deviation_a_tau_data = deviation_a_tau_data,
    deviation_b_tau_data = deviation_b_tau_data,
    deviation_a_xi_data = deviation_a_xi_data,
    deviation_b_xi_data = deviation_b_xi_data,
    tau_dev_plot = tau_dev_plot,
    xi_dev_plot = xi_dev_plot,
    beta0_cor_plot = beta0_cor_plot,
    gamma0_cor_plot = gamma0_cor_plot
  ),
  file = here::here("simulation-studies",
                    "hyperparameters.rds")
)
