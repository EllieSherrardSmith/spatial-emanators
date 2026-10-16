# Data preparation: turn the raw hut-row table (Hut Trial, Hut, Treatment, Time,
# au, du, df, f) into the two data structures the joint Stan model needs:
#   1) the within-hut compositional counts (used as-is, one row per hut-night)
#   2) the deterrence totals, aggregated by (Hut Trial, Time) across huts,
#      matching treatment-arm total catch against control-arm total catch
#      on the same night within the same trial.
#
# Requires: dplyr, tidyr

library(dplyr)
library(tidyr)

# raw: one row per hut-night, columns HutTrial, Hut, Treatment, Time, au, du, df, f
# T_ref: reference time (days) used by the eps40-parametrised Hill decay models
#        (nested_hill_eps40_model.stan, dirichlet_multinomial_hill_eps40_model.stan)
#        -- the time by which evidence indicates effectively no residual effect.
#        Not used by the plain nested_hill_model.stan / dirichlet_multinomial_hill_model.stan.
prepare_stan_data <- function(raw, T_ref = 400) {
  
  raw <- raw %>% mutate(N = au + du + df + f)
  
  # ---- deterrence aggregation: total catch and hut count per (trial, day, arm) ----
  agg <- raw %>%
    group_by(HutTrial, Time, Treatment) %>%
    summarise(total_catch = sum(N), n_huts = n_distinct(Hut), .groups = "drop")
  
  det <- agg %>%
    pivot_wider(
      names_from = Treatment,
      values_from = c(total_catch, n_huts),
      names_glue = "{.value}_{Treatment}"
    ) %>%
    # both arms must be present that trial-day to form a matched comparison
    filter(!is.na(total_catch_0), !is.na(total_catch_1)) %>%
    mutate(
      log_hut_ratio = log(n_huts_1 / n_huts_0),
      k = row_number()
    )
  
  if (nrow(det) == 0) {
    stop("No trial-days have both arms present -- check Treatment coding and grouping keys.")
  }
  
  # ---- link each raw hut-row to its (trial, day) deterrence group ----
  raw <- raw %>%
    left_join(det %>% select(HutTrial, Time, k), by = c("HutTrial", "Time"))
  
  if (any(is.na(raw$k))) {
    warning(
      "Some hut-rows have no matched deterrence group (their trial-day lacks one arm) ",
      "and were dropped from the compositional data to keep det_idx valid. ",
      "Consider whether those rows should be handled differently."
    )
    raw <- raw %>% filter(!is.na(k))
  }
  
  list(
    J = nrow(raw),
    y = as.matrix(raw[, c("au", "du", "df", "f")]),
    arm = raw$Treatment,
    time = raw$Time,
    
    K = nrow(det),
    Tcatch = det$total_catch_1,
    Ccatch = det$total_catch_0,
    time_det = det$Time,
    log_hut_ratio = det$log_hut_ratio,
    
    det_idx = raw$k,
    
    # weakly informative default; recompute from your chosen mean/effective sample size
    pi_prior_a = 0.699 * 6,
    pi_prior_b = (1 - 0.699) * 6,
    
    # only consumed by the eps parametrised Hill decay models
    T_ref = T_ref
    
    
  )
}

# Example usage:
# raw <- read.csv("hut_trial_data.csv")
# stan_data <- prepare_stan_data(raw)
# library(cmdstanr)
# mod <- cmdstan_model("nested_mixture_model.stan")
# fit <- mod$sample(data = stan_data, chains = 4, iter_warmup = 1000, iter_sampling = 1000)
