# Simulate dummy hut-trial data from the same generative process assumed by
# the joint Stan model (deterrence + latent would-be-biter/non-biter mixture
# + nested repelled / killed-before-feeding / dies-after-feeding tree).
#
# Output is a data frame in exactly the raw format the trial data was
# described in: HutTrial, Hut, Treatment, Time, au, du, df, f -- ready to
# hand straight to prepare_stan_data() from prepare_stan_data.R.
#
# This is illustrative test data for exercising the fitting pipeline end to
# end and checking that the Stan model recovers known truth -- it is not a
# calibrated simulation study of any real intervention.

library(dplyr)

simulate_hut_trial_data <- function(
    n_trials   = 1,
    n_huts     = 16,          # kept even so each arm gets an equal share of huts
    n_days     = 32,
    seed       = 1,
    rotate_huts = FALSE,      # TRUE: Latin-square-style rotation, which hut is "treatment" alternates
    #       by day. FALSE: each hut keeps the same arm for the whole trial.
    
    # ---- true parameter values (edit these to explore different scenarios) ----
    pi_true    = 0.65,       # P(would bite absent intervention)
    
    b_R = c(b0 = -1.5, b1 = 1.0, b2 = -0.02),   # repelled node, logit scale
    b_K = c(b0 = -2.0, b1 = 1.5, b2 = -0.03),   # killed-before-feeding node
    b_M = c(b0 = -1.0, b1 = 0.8, b2 = -0.01),   # dies-after-feeding node
    
    rho_R = 0.05, rho_K = 0.05, rho_M = 0.05,              # within-hut overdispersion, ICC scale
    
    g0_true   = -0.37,       # deterrence intercept, log-odds scale (~30% deterrence at day 1)
    g1_true   = 0.008,       # deterrence time trend (here: waning toward ~5% deterrence by day 40)
    rho_delta = 0.05,        # deterrence overdispersion, ICC scale
    
    baseline_catch_mean = 15 # mean nightly per-hut catch absent deterrence, before Poisson noise
) {
  set.seed(seed)
  
  logit_inv <- function(x) 1 / (1 + exp(-x))
  rbeta_mv  <- function(n, mean, phi) rbeta(n, mean * phi, (1 - mean) * phi)
  
  phi_R <- (1 - rho_R) / rho_R
  phi_K <- (1 - rho_K) / rho_K
  phi_M <- (1 - rho_M) / rho_M
  
  rows <- list()
  
  for (trial in seq_len(n_trials)) {
    
    # Fixed hut -> arm assignment, used throughout this trial when rotate_huts = FALSE:
    # half the huts control, half treatment, unchanged across every day.
    fixed_arm <- rep(0:1, length.out = n_huts)
    
    for (day in seq_len(n_days)) {
      
      arm_this_day <- if (rotate_huts) {
        # Latin-square-style rotation: which physical huts are "treatment"
        # alternates by day, so arm is not confounded with hut identity.
        (seq_len(n_huts) + day) %% 2   # 0/1, balanced when n_huts is even
      } else {
        fixed_arm
      }
      
      # Nightly mosquito-abundance multiplier shared by every hut that night
      # (weather/seasonal fluctuation) -- this is exactly the nuisance
      # quantity the matched-totals deterrence likelihood conditions away.
      night_multiplier <- rlnorm(1, meanlog = 0, sdlog = 0.3)
      
      # True deterrence at this trial-day, relative to the control reference.
      delta_day <- 1 - exp(g0_true + g1_true * day)
      
      for (hut in seq_len(n_huts)) {
        arm <- arm_this_day[hut]
        
        entry_rate <- baseline_catch_mean * night_multiplier *
          if (arm == 1) (1 - delta_day) else 1
        
        N_j <- rpois(1, lambda = max(entry_rate, 0.1))
        if (N_j == 0) next   # an empty catch just contributes no row for this hut-night
        
        # Hut-night-specific process probabilities: fixed-effects mean plus
        # genuine Beta-mixing overdispersion, matching the Stan model exactly.
        meanR <- logit_inv(b_R["b0"] + b_R["b1"] * arm + b_R["b2"] * arm * day)
        meanK <- logit_inv(b_K["b0"] + b_K["b1"] * arm + b_K["b2"] * arm * day)
        meanM <- logit_inv(b_M["b0"] + b_M["b1"] * arm + b_M["b2"] * arm * day)
        
        muR <- rbeta_mv(1, meanR, phi_R)
        muK <- rbeta_mv(1, meanK, phi_K)
        muM <- rbeta_mv(1, meanM, phi_M)
        
        kappa <- if (arm == 0) 0 else muK   # the two constraints from the methodology
        
        p <- c(
          au = (1 - pi_true) * (1 - kappa) + pi_true * muR,
          du = (1 - pi_true) * kappa       + pi_true * (1 - muR) * muK,
          df = pi_true * (1 - muR) * (1 - muK) * muM,
          f  = pi_true * (1 - muR) * (1 - muK) * (1 - muM)
        )
        p <- p / sum(p)   # guard against floating-point drift before rmultinom
        
        counts <- as.vector(rmultinom(1, N_j, p))
        
        rows[[length(rows) + 1]] <- data.frame(
          HutTrial  = trial,
          Hut       = hut,
          Treatment = arm,
          Time      = day,
          au = counts[1], du = counts[2], df = counts[3], f = counts[4]
        )
      }
    }
  }
  
  bind_rows(rows)
}

# ---- Quick-start: simulate, prepare, fit, and sanity-check recovery ----
#
raw <- simulate_hut_trial_data()
# write.csv(raw, "dummy_hut_trial_data.csv", row.names = FALSE)
#
source("R/multinomial_prep_stan_data.R")
stan_data <- prepare_stan_data(raw)
#
library(cmdstanr)
mod <- cmdstan_model("stan_code/multinomial_single.stan")   # or dirichlet_multinomial_mixture_model.stan
fit <- mod$sample(data = stan_data, chains = 4, iter_warmup = 1000, iter_sampling = 1000, seed = 1, parallel_chains = TRUE)

# fit_a <- mult_mod$sample(data = ms_data, 
#                          chains = 4, iter_warmup = 2000, iter_sampling = 3000, 
#                          seed = 1, parallel_chains = TRUE, 
#                          adapt_delta = 0.98, max_treedepth = 20)
# 
# 
# fit_b <- mult_mod$sample(data = ms_datb, 
#                          chains = 4, iter_warmup = 2000, iter_sampling = 3000, 
#                          seed = 1, parallel_chains = TRUE, 
#                          adapt_delta = 0.98, max_treedepth = 20)

## And interrogate the results
fit$summary(variables = c("pi", "b0_R", "b1_R", "b2_R", "b0_K", "b1_K", "b2_K", "b0_M", "b1_M", "b2_M", "g0", "g1"))

bayesplot::mcmc_trace(fit$draws(c("b0_R", "b1_R", "b2_R", "b0_K", "b1_K", "b2_K", "b0_M", "b1_M", "b2_M")))
bayesplot::mcmc_trace(fit$draws(c("pi", "g0", "g1")))

draws <- fit$draws(format = "df", variables = c("pi", "b0_R", "b1_R", "b2_R", "b0_K", "b1_K", "b2_K", "b0_M", "b1_M", "b2_M", "g0", "g1"))


### What do these values mean in practice?
# Generate mosquito hut outcomes

time = 1:32
trt = c(0,1)
eta_R = eta_K = eta_M = array(dim = c(nrow(draws_a),length(time),length(trt)))
eta_det = array(dim = c(nrow(draws_a),length(time)))

for(tr in 1:length(trt)){
  for(t in 1:length(time)){
    for(i in 1:nrow(draws_a)){
      eta_R[i,t,tr] <- draws_a$b0_R[i] +  draws_a$b1_R[i] *  trt[tr] + draws_a$b2_R[i] *  trt[tr] * time[t]  # vector of draws (length = posterior samples)
      eta_K[i,t,tr] <- draws_a$b0_K[i] +  draws_a$b1_K[i] *  trt[tr] + draws_a$b2_K[i] *  trt[tr] * time[t]  # vector of draws (length = posterior samples)
      eta_M[i,t,tr] <- draws_a$b0_M[i] +  draws_a$b1_M[i] *  trt[tr] + draws_a$b2_M[i] *  trt[tr] * time[t]  # vector of draws (length = posterior samples)
      eta_det[i,t] <- draws_a$g0[i] + draws_a$g1[i] * time[t]  # vector of draws (length = posterior samples)
      
    }  
  }
}

# inverse logit
p_R <- plogis(eta_R)
p_K <- plogis(eta_K)
p_M <- plogis(eta_M)
p_det <- plogis(eta_det)

prob_outcome_au = prob_outcome_du = prob_outcome_df = prob_outcome_f = array(dim = c(nrow(draws_a),length(time),length(trt)))

for(tr in 1:length(trt)){
  for(t in 1:length(time)){
    for(i in 1:nrow(draws_a)){
      
      kappa <- ifelse(trt[tr] == 0, 0, p_K[i, t, 1])
      prob_outcome_au[i, t, tr] <- (1 - draws_a$pi[i]) * (1 - kappa) + draws_a$pi[i] * p_R[i, t, tr]
      prob_outcome_du[i, t, tr] <- (1 - draws_a$pi[i]) * kappa + draws_a$pi[i] * (1 - p_R[i, t, tr]) * p_K[i, t, tr]
      prob_outcome_df[i, t, tr] <- draws_a$pi[i] * (1 - p_R[i, t, tr]) * (1 - p_K[i, t, tr]) * p_M[i, t, tr]
      prob_outcome_f[i, t, tr] <- draws_a$pi[i] * (1 - p_R[i, t, tr]) * (1 - p_K[i, t, tr]) * (1 - p_M[i, t, tr])
      
    }  
  }
}

outcome_prob_res <- bind_rows(
  colMeans(prob_outcome_au) |> data.frame(time) |> pivot_longer(cols = starts_with("X"), names_to = "arm", names_prefix = "X") |> 
    mutate(name = "au"),
  colMeans(prob_outcome_du) |> data.frame(time) |> pivot_longer(cols = starts_with("X"), names_to = "arm", names_prefix = "X") |> 
    mutate(name = "du"),
  colMeans(prob_outcome_df) |> data.frame(time) |> pivot_longer(cols = starts_with("X"), names_to = "arm", names_prefix = "X") |> 
    mutate(name = "df"),
  colMeans(prob_outcome_f) |> data.frame(time) |> pivot_longer(cols = starts_with("X"), names_to = "arm", names_prefix = "X") |> 
    mutate(name = "f")) |> 
  mutate(arm = as.integer(arm)-1)


# Hut outcome dummy data
hut_outcomes <- data.frame(stan_data$y,
                           time = stan_data$time,
                           arm = stan_data$arm,
                           n = rowSums(stan_data$y))


## Plot model fit against outcomes 
df_hut_outcomes |> 
  tidyr::pivot_longer(cols = c(au, du, df, f)) |> 
  ggplot() +
  geom_point(aes(x = time, y = value/n, color = name), size = 0.5) +
  geom_line(data = outcome_prob_res, aes(x = time, y = value, color = name)) +
  facet_grid(arm~name) +
  theme_bw() + theme(strip.background = element_rect(fill = "white"))


## !!! This is different from the probabilities within the decision tree !!!

## Plot each probability through time
## We could assume that the model probability is the difference between the control baseline and the full probability
mos_prob_res <- bind_rows(
  colMeans(p_R) |> data.frame(time) |> mutate(X3 = X2-X1) |> pivot_longer(cols = starts_with("X"), names_to = "arm", names_prefix = "X") |> 
    mutate(name = "R"),
  colMeans(p_K) |> data.frame(time) |> mutate(X3 = X2-X1) |> pivot_longer(cols = starts_with("X"), names_to = "arm", names_prefix = "X") |> 
    mutate(name = "K"),
  colMeans(p_M) |> data.frame(time) |> mutate(X3 = X2-X1) |> pivot_longer(cols = starts_with("X"), names_to = "arm", names_prefix = "X") |> 
    mutate(name = "M")) |> 
  mutate(arm = as.integer(arm)-1,
         arm = ifelse(arm == 2, "Diff", as.character(arm)))

mos_prob_res |> 
  ggplot() +
  geom_line(aes(x = time, y = value, color = name)) +
  facet_grid(arm~name) +
  theme_bw() + theme(strip.background = element_rect(fill = "white")) +
  geom_hline(yintercept = 0, linetype = "dashed")
