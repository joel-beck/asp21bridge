pacman::p_load(dplyr, ggplot2, purrr, tidyr)
library(asp21bridge)


set.seed(1234)
n <- 50
rho <- 0.9
pairs <- 10
num_sim <- 10000


create_pairs <- function(n, rho) {
  mvtnorm::rmvnorm(n = n, mean = c(0, 0), sigma = matrix(data = c(1, rho, rho, 1), nrow = 2))
}

# standardized covariates !!
create_design_matrix <- function(n, rho, pairs) {
  map(.x = 1:pairs, .f = ~ create_pairs(n = n, rho = rho)) %>%
    do.call(what = cbind, args = .) %>%
    magrittr::set_colnames(paste0("X", 1:(2 * pairs))) %>%
    as_tibble() %>%
    mutate(across(.cols = everything(), .fns = ~ (.x - mean(.x)) / sd(.x)))
}

create_data <- function(n, rho, pairs) {
  X <- create_design_matrix(n = n, rho = rho, pairs = pairs)
  Z <- X %>% select(X1, X3)

  beta <- rep(c(1, 0), times = pairs)
  gamma <- c(1, 1)

  y <- vector(mode = "numeric", length = n)

  for (i in seq_along(y)) {
    mu <- sum(X[i, ] * beta)
    sigma <- exp(sum(Z[i, ] * gamma))
    y[i] <- rnorm(n = 1, mean = mu, sd = sigma)
  }

  return(tibble(y = y) %>% bind_cols(X))
}

create_fit <- function(df, num_sim) {
  lmls(location = y ~ ., scale = ~ X1 + X3, data = df, light = FALSE) %>%
    mcmc(nsim = num_sim) %>%
    mcmc_ridge(num_sim = num_sim)
}

# helper function to combine estimates from all three models into single tibble
bind_results <- function(fit, pairs) {
  mcmc_ridge_res <- summary_complete(fit) %>%
    filter(stringr::str_detect(string = Parameter, pattern = "beta")) %>%
    select(Parameter, `Posterior Mean`) %>%
    mutate(sort_helper = readr::parse_number(Parameter)) %>%
    arrange(sort_helper)

  mcmc_res <- summary_complete(fit$mcmc) %>%
    slice(1:21) %>%
    mutate(sort_helper = if_else(condition = Parameter == "(Intercept)",
      true = 0, false = readr::parse_number(Parameter)
    )) %>%
    arrange(sort_helper) %>%
    pull(`Posterior Mean`)

  lmls_res <- coef(fit)$location

  tibble(
    Parameter = mcmc_ridge_res %>% pull(Parameter),
    mcmc_ridge = mcmc_ridge_res %>% pull(`Posterior Mean`),
    mcmc = mcmc_res,
    lmls = lmls_res,
    truth = c(0, rep(c(1, 0), times = pairs))
  )
}

df <- create_data(n = n, rho = rho, pairs = pairs)
fit <- create_fit(df = df, num_sim = num_sim)
complete_results <- bind_results(fit = fit, pairs = pairs)

plot_data_single_sim <- complete_results %>%
  mutate(Parameter = stringr::str_replace(
    string = Parameter, pattern = "_(.*)", replacement = "[\\1]"
  )) %>%
  mutate(Parameter = factor(Parameter) %>% forcats::fct_inorder() %>%
    forcats::fct_rev()) %>%
  mutate(truth = factor(truth)) %>%
  pivot_longer(cols = mcmc_ridge:lmls)

plot_single_sim <- plot_data_single_sim %>%
  ggplot(mapping = aes(x = value, y = Parameter, color = truth)) +
  geom_vline(xintercept = 1, color = "#2CA02CFF", linetype = "dashed") +
  geom_vline(xintercept = 0, color = "#CB2314FF", linetype = "dashed") +
  geom_point() +
  geom_segment(mapping = aes(xend = 0, yend = Parameter)) +
  facet_wrap(facets = vars(name)) +
  labs(
    x = "Coefficient Estimates", y = NULL,
    title = "Posterior Means / MLE for pairwise correlated Covariates",
    color = "True Coefficient Value:"
  ) +
  scale_x_continuous(breaks = scales::breaks_pretty()) +
  scale_y_discrete(
    labels = parse(text = levels(plot_data_single_sim$Parameter))
  ) +
  scale_color_manual(values = c("0" = "#CB2314FF", "1" = "#2CA02CFF")) +
  theme_light(base_size = 9) +
  theme(
    legend.position = "top",
    plot.title = element_text(hjust = 0.5)
  )

# sum of squared coefficient values (squared 2-norm) for each model and true coefficient value
coef_norms <- complete_results %>%
  pivot_longer(cols = mcmc_ridge:lmls, names_to = "model") %>%
  mutate(abs_value = abs(value)) %>%
  group_by(model, truth) %>%
  summarise(coefficient_norm = sum(value^2)) %>%
  mutate(truth = if_else(truth == 0, "truth = 0", "truth = 1")) %>%
  pivot_wider(names_from = truth, values_from = coefficient_norm)


readr::write_rds(
  x = list(plot_single_sim = plot_single_sim, coef_norms = coef_norms),
  file = here::here("simulation-studies", "redundant-regressors.rds")
)
