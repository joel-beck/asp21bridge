pacman::p_load(readr, here, dplyr)

acceptance_simulation <- readr::read_rds(
  file = here::here("reports", "acceptance_simulation.rds")
)

# two highest and two lowest acceptance rates for each value of num_sim
acceptance_simulation %>%
  group_by(num_sim) %>%
  arrange(acceptance_rate) %>%
  slice(1:2, n() - 1, n()) %>%
  relocate(num_sim) %>%
  print(n = Inf)


# acceptance rates between 0.25 and 0.4
acceptance_simulation %>%
  filter(between(acceptance_rate, 0.25, 0.4)) %>%
  arrange(num_sim, prop_var) %>%
  relocate(num_sim) %>%
  print(n = Inf)
