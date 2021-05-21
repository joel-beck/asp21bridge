pacman::p_load(dplyr, tidyr, purrr, readr, here)
library(lslmbRidge)

fit <- lslm(
  location = y ~ x1 + x2 + z1 + z2, scale = ~ z1 + z2,
  data = toy_data, light = FALSE
)

# takes a long time
acceptance_simulation <- tidyr::crossing(
  prop_var = c(seq(from = 0.1, to = 0.9, by = 0.1), seq(from = 1, to = 10, by = 1)),
  num_sim = seq(from = 200, to = 1000, by = 200)
) %>%
  mutate(acceptance_rate = purrr::map2_dbl(
    .x = prop_var, .y = num_sim,
    .f = ~ gibbs_sampler(m = fit, prop_var = .x, num_sim = .y) %>%
      purrr::pluck("mcmc_ridge", 2)
  ))

write_rds(
  file = here::here("reports", "acceptance_simulation.rds"),
  x = acceptance_simulation
)
