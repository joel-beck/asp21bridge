### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### Setup                                                                   ####

pacman::p_load(
  dplyr,
  purrr,
  tidyr,
  ggplot2,
  scales,
  patchwork,
  latex2exp
)
library(asp21bridge)


### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### Simulate Data                                                           ####
set.seed(123)

x1 <- rnorm(n = 50, mean = 1, sd = 1)
x2 <- rnorm(n = 50, mean = 2, sd = 1)
z1 <- rnorm(n = 50, mean = 5, sd = 1)
z2 <- rnorm(n = 50, mean = 3, sd = 1)

X <- cbind(x1, x2)
Z <- cbind(z1, z2)

X_stand <- apply(X, MARGIN = 2, function(x) (x - mean(x)) / sd(x))
Z_stand <- apply(Z, MARGIN = 2, function(x) (x - mean(x)) / sd(x))

beta <- c(8, 2)
gamma <- c(3, 3)

y <- vector(mode = "numeric", length = 50)
for (i in seq_along(y)) {
  mu <- sum(X_stand[i, ] * beta)
  sigma <- exp(sum(Z_stand[i, ] * gamma))
  y[i] <- rnorm(n = 1, mean = mu, sd = sigma)
}

SIMS <- FALSE

if (SIMS) {
  df_tau <- crossing(
    a_tau = 2^(-4:8),
    b_tau = 2^(-4:8)
  ) %>%
    rowwise() %>%
    mutate(fit = list(
      mcmc_ridge(
        X = X_stand, Z = Z_stand, y = y, beta_start = c(1, 1),
        gamma_start = c(1, 1), num_sim = 1000,
        a_tau = a_tau, b_tau = b_tau, a_xi = 1, b_xi = 1
      )
    )) %>%
    mutate(summary = list(summary_complete(fit))) %>%
    unnest(summary) %>%
    select(a_tau, b_tau, Parameter, `Posterior Mean`) %>%
    filter(Parameter %in% c("beta_1", "beta_2", "gamma_1", "gamma_2")) %>%
    mutate(beta = ifelse(stringr::str_detect(Parameter, pattern = "beta"), "beta", "gamma"))


  df_xi <- crossing(
    a_xi = 2^(-4:8),
    b_xi = 2^(-4:8)
  ) %>%
    rowwise() %>%
    mutate(fit = list(
      mcmc_ridge(
        X = X_stand, Z = Z_stand, y = y, beta_start = c(1, 1),
        gamma_start = c(1, 1), num_sim = 100,
        a_xi = a_xi, b_xi = b_xi, a_tau = 1, b_tau = 1
      )
    )) %>%
    mutate(summary = list(summary_complete(fit))) %>%
    unnest(summary) %>%
    select(a_xi, b_xi, Parameter, `Posterior Mean`) %>%
    filter(Parameter %in% c("beta_1", "beta_2", "gamma_1", "gamma_2")) %>%
    mutate(beta = ifelse(stringr::str_detect(Parameter, pattern = "beta"), "beta", "gamma"))


  readr::write_rds(
    x = list(df_tau = df_tau, df_xi = df_xi),
    file = here::here("simulation-studies", "hyper-pen-computations.rds")
  )
}

l <- readr::read_rds(
  file = here::here("simulation-studies", "hyper-pen-computations.rds")
)

df_tau <- l$df_tau
df_xi <- l$df_xi


# tau hypers on beta
p_tau_beta <- df_tau %>%
  group_by(a_tau, b_tau, beta) %>%
  summarise(norm = sum(`Posterior Mean`^2)) %>%
  filter(beta == "beta") %>%
  ggplot(mapping = aes(x = a_tau, y = b_tau, fill = norm)) +
  geom_tile() +
  scale_x_continuous(
    trans = "log2", breaks = 2^seq(-4, 8, 2),
    labels = c("1/16", "1/4", 1, 4, 16, 64, 256), expand = c(0, 0)
  ) +
  scale_y_continuous(
    trans = "log2", breaks = 2^seq(-4, 8, 2),
    labels = c("1/16", "1/4", 1, 4, 16, 64, 256), expand = c(0, 0)
  ) +
  scale_fill_distiller(palette = "RdBu") +
  labs(
    title = TeX("Squared Euclidean Norm of $\\beta$"), x = TeX("$a_\\tau$"),
    y = TeX("$b_\\tau$"), fill = NULL
  ) +
  theme_minimal(base_size = 9) +
  theme(plot.title = element_text(hjust = 0.5))


# tau hypers on gamma
p_tau_gamma <- df_tau %>%
  group_by(a_tau, b_tau, beta) %>%
  summarise(norm = sum(`Posterior Mean`^2)) %>%
  filter(beta == "gamma") %>%
  ggplot(mapping = aes(x = a_tau, y = b_tau, fill = norm)) +
  geom_tile() +
  scale_x_continuous(
    trans = "log2", breaks = 2^seq(-4, 8, 2),
    labels = c("1/16", "1/4", 1, 4, 16, 64, 256), expand = c(0, 0)
  ) +
  scale_y_continuous(
    trans = "log2", breaks = 2^seq(-4, 8, 2),
    labels = c("1/16", "1/4", 1, 4, 16, 64, 256), expand = c(0, 0)
  ) +
  scale_fill_distiller(palette = "RdBu") +
  labs(
    title = TeX("Squared Euclidean Norm of $\\gamma$"), x = TeX("$a_\\tau$"),
    y = NULL, fill = NULL
  ) +
  theme_minimal(base_size = 9) +
  theme(plot.title = element_text(hjust = 0.5))


# xi hypers on beta
p_xi_beta <- df_xi %>%
  group_by(a_xi, b_xi, beta) %>%
  summarise(norm = sum(`Posterior Mean`^2)) %>%
  filter(beta == "beta") %>%
  ggplot(mapping = aes(x = a_xi, y = b_xi, fill = norm)) +
  geom_tile() +
  scale_x_continuous(
    trans = "log2", breaks = 2^seq(-4, 8, 2),
    labels = c("1/16", "1/4", 1, 4, 16, 64, 256), expand = c(0, 0)
  ) +
  scale_y_continuous(
    trans = "log2", breaks = 2^seq(-4, 8, 2),
    labels = c("1/16", "1/4", 1, 4, 16, 64, 256), expand = c(0, 0)
  ) +
  scale_fill_distiller(palette = "RdBu") +
  labs(fill = NULL, x = TeX("$a_\\xi$"), y = TeX("$b_\\xi$")) +
  theme_minimal(base_size = 9)


# xi hypers on gamma
p_xi_gamma <- df_xi %>%
  group_by(a_xi, b_xi, beta) %>%
  summarise(norm = sum(`Posterior Mean`^2)) %>%
  filter(beta == "gamma") %>%
  ggplot(mapping = aes(x = a_xi, y = b_xi, fill = norm)) +
  geom_tile() +
  scale_x_continuous(
    trans = "log2", breaks = 2^seq(-4, 8, 2),
    labels = c("1/16", "1/4", 1, 4, 16, 64, 256), expand = c(0, 0)
  ) +
  scale_y_continuous(
    trans = "log2", breaks = 2^seq(-4, 8, 2),
    labels = c("1/16", "1/4", 1, 4, 16, 64, 256), expand = c(0, 0)
  ) +
  scale_fill_distiller(palette = "RdBu") +
  labs(y = NULL, fill = NULL, x = TeX("$a_\\xi$")) +
  theme_minimal(base_size = 9)

readr::write_rds(
  x = list(
    p_tau_beta = p_tau_beta, p_tau_gamma = p_tau_gamma,
    p_xi_beta = p_xi_beta, p_xi_gamma = p_xi_gamma
  ),
  file = here::here("simulation-studies", "hyper-penalty.rds")
)

l <- readr::read_rds(file = here::here("simulation-studies", "hyper-penalty.rds"))

(l$p_tau_beta | l$p_tau_gamma) / (l$p_xi_beta | l$p_xi_gamma)
