# Data preparation for multinomial_hill_eps_seasonal.stan. Identical to
# prepare_stan_data() in prepare_stan_data.R, plus the inputs for the seasonal
# biting probability pi(t):
#   period      -- days per seasonal cycle (e.g. 365)
#   n_harmonics -- number of Fourier harmonics H (0 = constant pi)
#   time_grid   -- trial days on which the model reports pi(t) (pi_grid in the output)
#
# Season is indexed by hut-trial day (the Time column), so no calendar dates are
# needed: trial day 1 is simply the phase origin of the Fourier cycle.
#
# Requires: dplyr, tidyr

library(dplyr)
library(tidyr)

# Fourier design matrix, columns [cos_1, sin_1, cos_2, sin_2, ...].
# Mirrors fourier_basis() in multinomial_hill_eps_seasonal.stan -- keep them in sync.
fourier_basis <- function(t, period, n_harmonics) {
  X <- matrix(0, nrow = length(t), ncol = 2 * n_harmonics)
  for (h in seq_len(n_harmonics)) {
    X[, 2 * h - 1] <- cos(2 * pi * h * t / period)
    X[, 2 * h]     <- sin(2 * pi * h * t / period)
  }
  X
}

# pi(t) for every posterior draw: returns an n_draws x length(t) matrix.
# alpha0: vector of draws; beta_pi: n_draws x (2*n_harmonics) matrix of draws.
seasonal_pi <- function(alpha0, beta_pi, t, period, n_harmonics) {
  X <- fourier_basis(t, period, n_harmonics)
  plogis(alpha0 + as.matrix(beta_pi) %*% t(X))
}

# raw: one row per hut-night, columns HutTrial, Hut, Treatment, Time, au, du, df, f
prepare_stan_data_seasonal <- function(raw, T_ref = 40, period = 365, n_harmonics = 2,
                                       time_grid = NULL) {

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

  if (is.null(time_grid)) time_grid <- seq_len(max(raw$Time))

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

    T_ref = T_ref,

    # seasonal biting probability
    period = period,
    n_harmonics = n_harmonics,
    N_grid = length(time_grid),
    time_grid = as.numeric(time_grid)
  )
}

# Example usage:
# source("R/prepare_stan_data_seasonal.R")
# stan_data <- prepare_stan_data_seasonal(raw, T_ref = 400, period = 365, n_harmonics = 2)
# library(cmdstanr)
# mod <- cmdstan_model("stan_code/multinomial_hill_eps_seasonal.stan")
# fit <- mod$sample(data = stan_data, chains = 4, iter_warmup = 1000, iter_sampling = 1000)
