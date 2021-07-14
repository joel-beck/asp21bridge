#   ____________________________________________________________________________
#   Setup                                                                   #####

pacman::p_load(dplyr, ggplot2, forcats, purrr, furrr)
library(asp21bridge)

create_data <- function(outcome_dist = c("norm", "t", "unif"), n) {
  outcome_dist <- match.arg(outcome_dist)

  # used for both location and scale
  x1 <- rnorm(n = n, mean = 5, sd = 4)
  x2 <- rexp(n = n, rate = 5)

  # used only for location
  x3 <- runif(n = n, min = -2, max = 12)
  x4 <- rbinom(n = n, size = 1, prob = 0.3)

  # used only for scale
  z3 <- rt(n = n, df = 10)

  X <- cbind(x1, x2, x3, x4)
  Z <- cbind(x1, x2, z3)

  # true coefficient values
  beta <- c(-3, -1, -1, 2)
  gamma <- c(1, 2, 3)

  y <- vector(mode = "numeric", length = n)

  # simulate different observation models
  for (i in seq_along(y)) {
    mu <- sum(X[i, ] * beta)
    sigma <- exp(sum(Z[i, ] * gamma))

    if (outcome_dist == "norm") {
      y[i] <- rnorm(n = 1, mean = mu, sd = sigma)
    } else if (outcome_dist == "t") {
      df <- 5
      y[i] <- mu + (sigma * sqrt((df - 2) / df)) * rt(n = 1, df = df)
    } else if (outcome_dist == "unif") {
      y[i] <- mu + (sigma * sqrt(12)) * runif(n = 1)
    }
  }

  return(tibble(y = y, x1 = x1, x2 = x2, x3 = x3, x4 = x4, z3 = z3))
}

# fit lmls(), mcmc() and mcmc_ridge() to given data
create_fit <- function(df, num_sim) {
  lmls(
    location = y ~ x1 + x2 + x3 + x4, scale = ~ x1 + x2 + z3,
    data = df, light = FALSE
  ) %>%
    mcmc(nsim = num_sim) %>%
    mcmc_ridge(num_sim = num_sim)
}

# helper function to combine estimates from all three models into single tibble
bind_results <- function(model) {
  mcmc_ridge_res <- summary_complete(model) %>%
    filter(stringr::str_detect(string = Parameter, pattern = "beta|gamma")) %>%
    select(Parameter, `Posterior Mean`)

  mcmc_res <- summary_complete(model$mcmc) %>% pull(`Posterior Mean`)
  lmls_res <- c(coef(model)$location, coef(model)$scale)
  truth <- c(0, -3, -1, -1, 2, 0, 1, 2, 3)

  tibble(
    Parameter = mcmc_ridge_res %>% pull(Parameter),
    mcmc_ridge = mcmc_ridge_res %>% pull(`Posterior Mean`),
    mcmc = mcmc_res,
    lmls = lmls_res,
    truth = truth
  )
}

# creates data, fits models and collects results
create_results <- function(seed = NULL, n, num_sim) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  full_data <- tibble(outcome_dist = c("norm", "t", "unif")) %>%
    rowwise() %>%
    mutate(df = list(create_data(outcome_dist = outcome_dist, n = n))) %>%
    mutate(model = list(create_fit(df = df, num_sim = num_sim))) %>%
    mutate(results = list(bind_results(model = model)))

  full_results <- full_data$results %>%
    purrr::set_names(nm = c("norm", "t", "unif")) %>%
    bind_rows(.id = "outcome_dist")

  return(full_results)
}

# returns list with numerical results and corresponding plot
show_results <- function(seed = NULL, n, num_sim) {
  results <- create_results(seed = seed, n = n, num_sim = num_sim)

  numerical <- results %>%
    arrange(Parameter)

  plot_data_single_sim <- results %>%
    mutate(outcome_dist = factor(outcome_dist) %>%
      fct_recode(
        "y ~ Normal" = "norm",
        "y ~ t" = "t",
        "y ~ Uniform" = "unif"
      )) %>%
    tidyr::pivot_longer(cols = mcmc_ridge:lmls) %>%
    mutate(Parameter = stringr::str_replace(
      string = Parameter, pattern = "_(.*)", replacement = "[\\1]"
    )) %>%
    mutate(Parameter = factor(Parameter) %>% fct_rev())



  plot_single_sim <- plot_data_single_sim %>%
    ggplot(mapping = aes(x = value, y = Parameter, color = name)) +
    geom_point(
      mapping = aes(x = truth), color = "grey80", fill = "transparent",
      size = 6, shape = 21
    ) +
    geom_point(
      position = position_dodge(width = 0.5), size = 1
    ) +
    facet_wrap(facets = vars(outcome_dist), scales = "free_x") +
    labs(
      title = "Posterior Means / MLE for (misspecified) Regression Models",
      subtitle = "True coefficient values are marked by grey circles",
      x = "Coefficient Estimate", y = NULL, color = NULL
    ) +
    scale_x_continuous(breaks = scales::breaks_pretty()) +
    scale_y_discrete(
      labels = parse(text = levels(plot_data_single_sim$Parameter))
    ) +
    guides(color = guide_legend(reverse = TRUE)) +
    theme_light(base_size = 9) +
    theme(
      legend.position = "top",
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      panel.grid.minor = element_blank()
    )

  list(results = numerical, plot = plot_single_sim)
}

#   ____________________________________________________________________________
#   Data and Plot for Single Simulation                                     ####

seed_2 <- show_results(seed = 2, n = 50, num_sim = 1000)

data_single_sim <- seed_2$results
plot_single_sim <- seed_2$plot


#   ____________________________________________________________________________
#   Data and Plot for Many Simulations                                      ####

SAVE_SIMS <- FALSE

if (SAVE_SIMS) {
  plan(multisession, workers = 8)
  safe_show_results <- safely(show_results)

  full_results <- tibble(id = 1:50) %>%
    mutate(samples = future_map(
      .x = id,
      .f = ~ safe_show_results(n = 50, num_sim = 1000),
      .options = furrr_options(seed = 1)
    )) %>%
    mutate(result = map(.x = samples, .f = ~ .x$result$results))

  plan(sequential)

  readr::write_rds(
    x = full_results,
    file = here::here("simulation-studies", "out-dist-computation.rds")
  )
}

# calculate bias (average error) and standard deviation of error, stratified by
# outcome distribution, model and Parameter
many_sims <- readr::read_rds(
  file = here::here("simulation-studies", "out-dist-computation.rds")
)

data_many_sims <- bind_rows(many_sims$result) %>%
  tidyr::pivot_longer(cols = mcmc_ridge:lmls) %>%
  mutate(error = value - truth) %>%
  group_by(outcome_dist, name, Parameter) %>%
  summarise(bias = mean(error), se = sd(value), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = name, values_from = c(bias, se))

plot_data_many_sims <- bind_rows(many_sims$result) %>%
  mutate(outcome_dist = factor(outcome_dist) %>%
    fct_recode(
      "y ~ Normal" = "norm",
      "y ~ t" = "t",
      "y ~ Uniform" = "unif"
    )) %>%
  tidyr::pivot_longer(cols = mcmc_ridge:lmls) %>%
  mutate(error = value - truth) %>%
  group_by(outcome_dist, name, Parameter) %>%
  summarise(
    mean_estimate = mean(value),
    truth = truth[1],
    se = sd(value),
    lower = quantile(value, probs = 0.05),
    upper = quantile(value, probs = 0.95),
    .groups = "drop"
  ) %>%
  mutate(Parameter = stringr::str_replace(
    string = Parameter, pattern = "_(.*)", replacement = "[\\1]"
  )) %>%
  mutate(Parameter = factor(Parameter) %>% fct_rev())


plot_many_sims <- plot_data_many_sims %>%
  ggplot(mapping = aes(x = mean_estimate, y = Parameter, color = name)) +
  geom_point(
    mapping = aes(x = truth), color = "grey80", fill = "transparent",
    size = 6, shape = 21
  ) +
  geom_errorbar(
    mapping = aes(
      # xmin = mean_estimate - 1.96 * se, xmax = mean_estimate + 1.96 * se,
      xmin = lower, xmax = upper, color = name
    ),
    linetype = "dashed", width = 0.5, position = position_dodge(width = 0.5),
    show.legend = FALSE
  ) +
  # geom_boxplot(mapping = aes(x = value), data = boxplot_data_many_sims,
  #              position = position_dodge(width = 0.5)) +
  geom_point(
    mapping = aes(x = mean_estimate), size = 1,
    position = position_dodge(width = 0.5)
  ) +
  facet_wrap(facets = vars(outcome_dist), scales = "free_x") +
  scale_y_discrete(
    labels = parse(text = levels(plot_data_many_sims$Parameter))
  ) +
  guides(color = guide_legend(reverse = TRUE)) +
  labs(
    title = "Empirical 90% Confidence Intervals for Posterior Mean Estimates",
    subtitle = "True coefficient values are marked by grey circles",
    x = "Estimate", y = NULL, color = NULL
  ) +
  theme_light(base_size = 9) +
  theme(
    legend.position = "top",
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )

plot_many_sims



#   ____________________________________________________________________________
#   Data for Second Report                                                  ####

readr::write_rds(
  x = list(
    data_single_sim = data_single_sim,
    plot_single_sim = plot_single_sim,
    data_many_sims = data_many_sims,
    plot_many_sims = plot_many_sims
  ),
  file = here::here("simulation-studies", "outcome-distribution.rds")
)
