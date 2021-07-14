#   ____________________________________________________________________________
#   Setup                                                                   ####

pacman::p_load(dplyr, purrr, ggplot2, tidyr, forcats, furrr, latex2exp)
library(asp21bridge)

create_data <- function(rho, n) {
  X <- mvtnorm::rmvnorm(
    n = 50, mean = c(-5, 2, 0),
    sigma = matrix(c(
      1, rho, rho,
      rho, 3, rho,
      rho, rho, 5
    ), nrow = 3)
  )
  beta <- c(3, -1, 1)

  z1 <- 0.8 * X[, 1] + 0.2 * X[, 2]
  z2 <- X[, 2] - 0.5 * X[, 3]
  Z <- cbind(z1, z2)

  gamma <- c(2, 0)

  y <- vector(mode = "numeric", length = n)

  for (i in seq_along(y)) {
    mu <- sum(X[i, ] * beta)
    sigma <- exp(sum(Z[i] * gamma))
    y[i] <- rnorm(n = 1, mean = mu, sd = sigma)
  }

  tibble(y = y, x1 = X[, 1], x2 = X[, 2], x3 = X[, 3], z1 = z1, z2 = z2)
}

create_fit <- function(data, mh_location, num_sim) {
  data %>%
    lmls(location = y ~ x1 + x2 + x3, scale = ~ z1 + z2, light = FALSE) %>%
    mcmc(nsim = num_sim) %>%
    mcmc_ridge(num_sim = num_sim, mh_location = mh_location)
}

extract_coefs <- function(model) {
  mcmc_ridge_res <- summary_complete(model) %>%
    filter(stringr::str_detect(string = Parameter, pattern = "beta|gamma")) %>%
    select(Parameter, `Posterior Mean`)

  mcmc_res <- summary_complete(model$mcmc) %>% pull(`Posterior Mean`)
  lmls_res <- c(coef(model)$location, coef(model)$scale)
  truth <- c(0, 3, -1, 1, 0, 2, 0)

  tibble(
    Parameter = mcmc_ridge_res %>% pull(Parameter),
    mcmc_ridge = mcmc_ridge_res %>% pull(`Posterior Mean`),
    mcmc = mcmc_res,
    lmls = lmls_res,
    truth = truth
  )
}

extract_accep_rate <- function(model) {
  acc_rate_list <- model$mcmc_ridge$acceptance_rate
  if (length(acc_rate_list) == 1) {
    return(c(acc_rate_loc = NA, acc_rate_scale = acc_rate_list))
  } else {
    return(c(acc_rate_loc = acc_rate_list[[1]], acc_rate_scale = acc_rate_list[[2]]))
  }
}



#   ____________________________________________________________________________
#   Data and Plot for Single Simulation                                     ####

set.seed(111)

# Data
full_data <- crossing(rho = c(-0.5, 0, 0.9), mh_loc = c(FALSE, TRUE)) %>%
  mutate(cor_data = map(.x = rho, .f = ~ create_data(rho = .x, n = 50))) %>%
  mutate(model = map2(
    .x = cor_data,
    .y = mh_loc,
    .f = ~ create_fit(data = .x, mh_location = .y, num_sim = 1000)
  )) %>%
  mutate(coefs = map(.x = model, .f = extract_coefs)) %>%
  mutate(accep_rate = map(.x = model, .f = extract_accep_rate)) %>%
  unnest_wider(accep_rate)

df_coefs <- full_data %>%
  select(rho, mh_loc, coefs:acc_rate_scale) %>%
  unnest(coefs)


# Plot Data
plot_data_single_sim <- df_coefs %>%
  mutate(acc_label = if_else(
    condition = is.na(acc_rate_loc),
    true = stringr::str_glue("acc_scale = {acc_rate_scale}"),
    false = stringr::str_glue(
      "acc_scale = {acc_rate_scale}\n acc_loc = {acc_rate_loc}"
    )
  )) %>%
  mutate(rho = factor(rho) %>% fct_recode(
    "rho = -0.5" = "-0.5",
    "rho = 0" = "0",
    "rho = 0.9" = "0.9"
  )) %>%
  mutate(mh_loc = factor(mh_loc) %>% fct_recode(
    "Multivariate Normal" = "FALSE",
    "Metropolis - Hastings" = "TRUE"
  )) %>%
  mutate(Parameter = stringr::str_replace(
    string = Parameter, pattern = "_(.*)", replacement = "[\\1]"
  )) %>%
  mutate(Parameter = factor(Parameter) %>% fct_rev()) %>%
  pivot_longer(cols = mcmc_ridge:lmls)


# Plot
plot_single_sim <- plot_data_single_sim %>%
  ggplot(mapping = aes(x = value, y = Parameter)) +
  geom_point(
    mapping = aes(x = truth), color = "grey80", fill = "transparent",
    size = 6, shape = 21
  ) +
  geom_point(
    mapping = aes(color = name), position = position_dodge(width = 0.5),
    size = 1
  ) +
  geom_label(
    mapping = aes(label = acc_label), x = 20, y = 1.5, color = "grey50", size = 2
  ) +
  facet_grid(rows = vars(mh_loc), cols = vars(rho)) +
  labs(
    title = "Model Performance for different Predictor Correlation Structures",
    subtitle = "True coefficient values are indicated by grey circles",
    x = "Estimate", y = NULL,
    color = NULL
  ) +
  guides(color = guide_legend(reverse = TRUE)) +
  scale_x_continuous(breaks = scales::breaks_pretty()) +
  scale_y_discrete(
    labels = parse(text = levels(plot_data_single_sim$Parameter))
  ) +
  theme_light(base_size = 9) +
  theme(
    legend.position = "top",
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )

plot_single_sim







#   ____________________________________________________________________________
#   Data and Plot for Many Simulations                                      ####

SAVE_SIMS <- FALSE

if (SAVE_SIMS) {
  single_sim <- function(x) {
    crossing(rho = c(-0.5, 0, 0.9), mh_loc = c(FALSE, TRUE)) %>%
      mutate(cor_data = map(.x = rho, .f = ~ create_data(rho = .x, n = 50))) %>%
      mutate(model = map2(
        .x = cor_data,
        .y = mh_loc,
        .f = ~ create_fit(data = .x, mh_location = .y, num_sim = 1000)
      )) %>%
      mutate(coefs = map(.x = model, .f = extract_coefs)) %>%
      select(rho, mh_loc, coefs) %>%
      unnest(coefs)
  }

  plan(multisession, workers = 8)

  many_sims <- tibble(id = 1:50) %>%
    mutate(results = future_map(
      .x = id, .f = single_sim,
      .options = furrr_options(seed = 1)
    ))

  plan(sequential)

  readr::write_rds(
    x = many_sims,
    file = here::here("simulation-studies", "regr-cor-computation.rds")
  )
}

many_sims <- readr::read_rds(
  file = here::here("simulation-studies", "regr-cor-computation.rds")
)


# Data
bias_se_results <- bind_rows(many_sims$results) %>%
  filter(!mh_loc) %>%
  select(-mh_loc) %>%
  pivot_longer(cols = mcmc_ridge:lmls) %>%
  mutate(error = value - truth) %>%
  group_by(rho, Parameter, name) %>%
  summarize(
    bias = mean(error),
    se = sd(value), .groups = "drop"
  ) %>%
  pivot_wider(names_from = name, values_from = c(bias, se))

bias_se_results %>%
  filter(stringr::str_detect(Parameter, pattern = "beta")) %>%
  arrange(Parameter)

bias_se_results %>%
  filter(stringr::str_detect(Parameter, pattern = "gamma")) %>%
  arrange(Parameter)


# Plot Data
plot_data <- bind_rows(many_sims$results) %>%
  filter(!mh_loc) %>%
  select(rho, Parameter, mcmc_ridge, truth) %>%
  mutate(error = mcmc_ridge - truth) %>%
  group_by(rho, Parameter) %>%
  summarize(
    mean_estimate = mean(mcmc_ridge),
    truth = truth[1],
    se = sd(mcmc_ridge),
    lower = quantile(mcmc_ridge, probs = 0.05),
    upper = quantile(mcmc_ridge, probs = 0.95),
    .groups = "drop"
  ) %>%
  mutate(rho = factor(rho) %>% fct_recode(
    "rho = -0.5" = "-0.5",
    "rho = 0" = "0",
    "rho = 0.9" = "0.9"
  )) %>%
  mutate(Parameter = stringr::str_replace(
    string = Parameter, pattern = "_(.*)", replacement = "[\\1]"
  )) %>%
  mutate(Parameter = factor(Parameter) %>% fct_rev())


# Plot
plot_many_sims <- plot_data %>%
  ggplot(mapping = aes(x = mean_estimate, y = Parameter)) +
  geom_point(
    mapping = aes(x = truth), color = "grey80", fill = "transparent",
    size = 6, shape = 21
  ) +
  geom_errorbar(
    mapping = aes(
      # xmin = mean_estimate - 1.96 * se, xmax = mean_estimate + 1.96 * se
      xmin = lower, xmax = upper
    ),
    color = "grey10", linetype = "dashed", width = 0.3
  ) +
  geom_point(mapping = aes(x = mean_estimate)) +
  facet_wrap(facets = vars(rho)) +
  scale_y_discrete(
    labels = parse(text = levels(plot_data$Parameter))
  ) +
  labs(
    title = "Empirical 90% Confidence Intervals for Posterior Mean Estimates",
    subtitle = "True coefficient values are marked by grey circles",
    x = "Estimate", y = NULL
  ) +
  theme_light(base_size = 9) +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )

plot_many_sims



#   ____________________________________________________________________________
#   Data for Second Report                                                  ####

readr::write_rds(
  x = list(
    data_single_sim = df_coefs,
    plot_single_sim = plot_single_sim,
    data_many_sims = bias_se_results,
    plot_many_sims = plot_many_sims
  ),
  file = here::here("simulation-studies", "regressor-correlation.rds")
)
