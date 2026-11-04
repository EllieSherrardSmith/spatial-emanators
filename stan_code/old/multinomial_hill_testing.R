# Simulate dummy hut-trial data from the generative process assumed by
# nested_hill_model.stan: background mortality split into a before-/after-
# average-feeding-time hazard (nu_u, nu_f, both arms), intervention-attributable
# hazards R (repelled before feeding), K (killed before feeding), M (killed
# after feeding) -- each a Hill/log-logistic decay in time and structurally
# zero in control -- plus a Hill-decay deterrence effect on hut entry.
#
# Output is a data frame in the same raw format as before: HutTrial, Hut,
# Treatment, Time, au, du, df, f -- ready for prepare_stan_data() from
# prepare_stan_data.R (that function is unchanged; only the underlying
# generative model and the Stan file have changed).

library(dplyr)

dtaMS = readxl::read_excel("C:/Users/esherrar/Documents/RProjects/spatial-emanators_1/data/Mosquito Shield_complied data_22 July 2026.xls", sheet = "mosquito shield_2022")
dtaMS$date2 = lubridate::ymd(dtaMS$date)
dtaMS <- dtaMS[order(dtaMS$date2), ]   # ensure sorted by date
dtaMS$time_index <- as.integer(dtaMS$date2 - min(dtaMS$date2)) + 1

for(i in 1:nrow(dtaMS)){
  dtaMS$total[i] = dtaMS$aliveunfed24[i] + dtaMS$alivefed24[i] +
    dtaMS$deadunfed24[i] + dtaMS$deadfed24[i]
}

dtaMS <- dtaMS[order(dtaMS$treatment, 
                     dtaMS$time_index, 
                     dtaMS$species), ]

dtaMS_summary <- dtaMS %>%
  group_by(treatment, time_index, species) %>%
  summarise(
    across(c("aliveunfed24",
             "alivefed24",
             "deadunfed24",
             "deadfed24",
             total), sum, na.rm = TRUE)
  )
dtaMS_summary = subset(dtaMS_summary, dtaMS_summary$species == "arabiensis")

##
##Ngufor (gambiae)
dtbMS = readxl::read_excel("C:/Users/esherrar/Documents/RProjects/spatial-emanators_1/data/Ngufor - MosquitoShield EHT/Benin_Mosquito Shield_Hut data_CN.xlsx", sheet = "An. gambiae SMC" )

dtbMS1 = subset(dtbMS, dtbMS$timepoint == "a.0 month")
dtbMS2 = subset(dtbMS, dtbMS$timepoint == "b.4 months")

dtbMS1$date2 = lubridate::ymd(dtbMS1$date)
dtbMS1 <- dtbMS1[order(dtbMS1$date2), ]   # ensure sorted by date
dtbMS1$time_index <- as.integer(dtbMS1$date2 - min(dtbMS1$date2)) + 1

dtbMS2$date2 = lubridate::ymd(dtbMS2$date)
dtbMS2 <- dtbMS2[order(dtbMS2$date2), ]   # ensure sorted by date
dtbMS2$time_index <- as.integer(dtbMS2$date2 - min(dtbMS2$date2)) + 1

for(i in 1:nrow(dtbMS1)){
  dtbMS1$total[i] = dtbMS1$unf_live[i] + dtbMS1$bf_live[i] +
    dtbMS1$unf_dead[i] + dtbMS1$bf_dead[i]
}

for(i in 1:nrow(dtbMS2)){
  dtbMS2$total[i] = dtbMS2$unf_live[i] + dtbMS2$bf_live[i] +
    dtbMS2$unf_dead[i] + dtbMS2$bf_dead[i]
}

dtbMS1 <- dtbMS1[order(dtbMS1$treatment, 
                       dtbMS1$time_index), ]
dtbMS2 <- dtbMS2[order(dtbMS2$treatment, 
                       dtbMS2$time_index), ]

dtbMS1$treatment = ifelse(dtbMS1$treatment == "a.control",0,1)
dtbMS2$treatment = ifelse(dtbMS2$treatment == "a.control",0,1)


simulate_hut_trial_data_hill <- function(
    n_trials    = 1,
    n_huts      = 6,          # kept even so each arm gets an equal share of huts
    n_days      = 40,
    seed        = 1,
    rotate_huts = TRUE,       # TRUE: Latin-square-style rotation of which hut is "treatment" by day.
    # FALSE: each hut keeps the same arm for the whole trial.
    
    # ---- true parameter values (edit these to explore different scenarios) ----
    pi_true   = 0.65,         # P(would bite absent intervention)
    nu_u_true = 0.05,         # background mortality hazard, before the average feeding time
    nu_f_true = 0.03,         # background mortality hazard, after the average feeding time
    
    # Treatment-only hazards: each is a Hill/log-logistic decay f0 * thalf^n / (thalf^n + t^n)
    R0_true = 0.50, thalf_R_true = 15, nR_true = 2.0,   # repelled before feeding
    K0_true = 0.30, thalf_K_true = 20, nK_true = 2.0,   # killed before feeding
    M0_true = 0.15, thalf_M_true = 25, nM_true = 1.5,   # killed after feeding
    
    rho_R = 0.05, rho_K = 0.05, rho_M = 0.05,           # within-hut overdispersion, ICC scale (treatment rows only)
    
    delta0_true = 0.30, thalf_d_true = 15, nd_true = 2.0,  # deterrence Hill decay
    rho_delta   = 0.05,                                    # deterrence overdispersion, ICC scale
    
    baseline_catch_mean = 15  # mean nightly per-hut catch absent deterrence, before Poisson noise
) {
  set.seed(seed)
  
  hill_decay <- function(t, f0, thalf, n) f0 * thalf^n / (thalf^n + t^n)
  rbeta_mv   <- function(n_draw, mean, phi) rbeta(n_draw, mean * phi, (1 - mean) * phi)
  
  phi_R <- (1 - rho_R) / rho_R
  phi_K <- (1 - rho_K) / rho_K
  phi_M <- (1 - rho_M) / rho_M
  
  rows <- list()
  
  for (trial in seq_len(n_trials)) {
    
    fixed_arm <- rep(0:1, length.out = n_huts)
    
    for (day in seq_len(n_days)) {
      
      arm_this_day <- if (rotate_huts) {
        (seq_len(n_huts) + day) %% 2
      } else {
        fixed_arm
      }
      
      # Nightly mosquito-abundance multiplier shared by every hut that night.
      night_multiplier <- rlnorm(1, meanlog = 0, sdlog = 0.3)
      
      # True deterrence at this trial-day (Hill decay), relative to the control reference.
      delta_day <- hill_decay(day, delta0_true, thalf_d_true, nd_true)
      
      for (hut in seq_len(n_huts)) {
        arm <- arm_this_day[hut]
        
        entry_rate <- baseline_catch_mean * night_multiplier *
          if (arm == 1) (1 - delta_day) else 1
        
        N_j <- rpois(1, lambda = max(entry_rate, 0.1))
        if (N_j == 0) next   # an empty catch just contributes no row for this hut-night
        
        if (arm == 1) {
          meanR <- hill_decay(day, R0_true, thalf_R_true, nR_true)
          meanK <- hill_decay(day, K0_true, thalf_K_true, nK_true)
          meanM <- hill_decay(day, M0_true, thalf_M_true, nM_true)
          
          R <- rbeta_mv(1, meanR, phi_R)
          K <- rbeta_mv(1, meanK, phi_K)
          M <- rbeta_mv(1, meanM, phi_M)
        } else {
          R <- 0; K <- 0; M <- 0   # structural constraint: R = K = M = 0 in control
        }
        
        A  <- 1 - (1 - nu_u_true) * (1 - K)   # combined hazard, before the average feeding time
        hF <- 1 - (1 - nu_f_true) * (1 - M)   # combined hazard, after the average feeding time
        
        p <- c(
          au = (1 - A) * (1 - hF) * ((1 - pi_true) + pi_true * R),   # alive, unfed
          du = A + (1 - A) * hF * ((1 - pi_true) + pi_true * R),      # dead,  unfed
          df = pi_true * (1 - R) * (1 - A) * hF,                     # dead,  fed
          f  = pi_true * (1 - R) * (1 - A) * (1 - hF)                # alive, fed
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
raw_hill <- simulate_hut_trial_data_hill()

real_hill = data.frame(HutTrial = rep(1, nrow(dtaMS)),
                       Hut = dtaMS$hut,
                       Treatment = dtaMS$treatment,
                       Time = dtaMS$time_index,
                       au = dtaMS$aliveunfed24,
                       du = dtaMS$deadunfed24,
                       df = dtaMS$deadfed24,
                       f = dtaMS$alivefed24)

real_hillB = data.frame(HutTrial = c(rep(1, nrow(dtbMS1)),rep(2, nrow(dtbMS2))),
                       Hut = c(dtbMS1$hut,dtbMS2$hut),
                       Treatment = c(dtbMS1$treatment,dtbMS2$treatment),
                       Time = c(dtbMS1$time_index,dtbMS2$time_index),
                       au = c(dtbMS1$unf_live,dtbMS2$unf_live),
                       du = c(dtbMS1$unf_dead,dtbMS2$unf_dead),
                       df = c(dtbMS1$bf_dead,dtbMS2$bf_dead),
                       f = c( dtbMS1$bf_live,dtbMS2$bf_live)
                       )
# write.csv(raw, "dummy_hut_trial_data_hill.csv", row.names = FALSE)
#
# source("prepare_stan_data.R")     # unchanged -- same raw format, same aggregation logic
## run multinomial_prep_stan_data()
stan_data_hill <- prepare_stan_data(real_hill)
stan_data_hill <- prepare_stan_data(real_hillB)
#
library(cmdstanr)
mod_hill <- cmdstan_model("stan_code/SE_multinomial_fits/multinomial_hill.stan")
fit_hill <- mod_hill$sample(data = stan_data_hill, chains = 4, iter_warmup = 1000, iter_sampling = 1000, seed = 1, parallel_chains = 4)

## And interrogate the results
fit_hill$summary(variables = c("pi", "nu_u", "nu_f",
                               "R0", "thalf_R", "nR",
                               "K0", "thalf_K", "nK",
                               "M0", "thalf_M", "nM",
                               "delta0", "thalf_d", "nd"))

bayesplot::mcmc_trace(fit_hill$draws(c("pi", "nu_u", "nu_f",
                                  "R0", "thalf_R", "nR",
                                  "K0", "thalf_K", "nK",
                                  "M0", "thalf_M", "nM",
                                  "delta0", "thalf_d", "nd")))

draws_hill <- fit_hill$draws(format = "df", variables = c("pi", "nu_u", "nu_f",
                                                     "R0", "thalf_R", "nR",
                                                     "K0", "thalf_K", "nK",
                                                     "M0", "thalf_M", "nM",
                                                     "delta0", "thalf_d", "nd"))


### What do these values mean in practice?
### Generate mosquito outcomes

time = 1:32
trt = c(0,1)

p_R = p_K = p_M = p_det = array(dim = c(nrow(draws_hill),length(time)))

for(t in 1:length(time)){
  for(i in 1:nrow(draws_hill)){
    p_R[i,t] <- draws_hill$R0[i] * (1 / (1 + (t/draws_hill$thalf_R[i])^draws_hill$nR[i]))
    p_K[i,t] <- draws_hill$K0[i] * (1 / (1 + (t/draws_hill$thalf_K[i])^draws_hill$nK[i]))
    p_M[i,t] <- draws_hill$M0[i] * (1 / (1 + (t/draws_hill$thalf_M[i])^draws_hill$nM[i]))
    p_det[i,t] <- draws_hill$delta0[i] * (1 / (1 + (t/draws_hill$thalf_d[i])^draws_hill$nd[i]))
  }  
}

p_A <- 1 - (1 - draws_hill$nu_u) * (1 - p_K)
p_hF <- 1 - (1 - draws_hill$nu_f) * (1 - p_M)

# Full model with SE impacts
prob_outcome_au = (1 - p_A) * (1 - p_hF) * ((1 - draws_hill$pi) + draws_hill$pi * p_R)
prob_outcome_du = p_A + (1 - p_A) * p_hF * ((1 - draws_hill$pi) + draws_hill$pi * p_R)
prob_outcome_df = draws_hill$pi * (1 - p_R) * (1 - p_A) * p_hF
prob_outcome_f = draws_hill$pi * (1 - p_R) * (1 - p_A) * (1 - p_hF)

# Simplified model without SE impacts
prob_outcome_au_c = (1 - draws_hill$nu_u) * (1 - draws_hill$nu_f) * (1 - draws_hill$pi)
prob_outcome_du_c = draws_hill$nu_u + (1 - draws_hill$nu_u) * draws_hill$nu_f * (1 - draws_hill$pi)
prob_outcome_df_c = draws_hill$pi * (1 - draws_hill$nu_u) * draws_hill$nu_f
prob_outcome_f_c = draws_hill$pi * (1 - draws_hill$nu_u) * (1 - draws_hill$nu_f)

# Treatment probabilities include time-dependency
outcome_prob_res_hill <- 
  bind_rows(
    data.frame(
      au = colMeans(prob_outcome_au),
      du = colMeans(prob_outcome_du),
      df = colMeans(prob_outcome_df),
      f = colMeans(prob_outcome_f)) |> 
      mutate(arm = 1,
             time = time)) |> 
  pivot_longer(cols = c("au", "du", "df", "f")) |> 
  mutate( 
    lower = c( 
    apply(prob_outcome_au, 2, quantile, probs = 0.025), 
    apply(prob_outcome_du, 2, quantile, probs = 0.025), 
    apply(prob_outcome_df, 2, quantile, probs = 0.025), 
    apply(prob_outcome_f, 2, quantile, probs = 0.025) 
    ), 
    upper = c( 
      apply(prob_outcome_au, 2, quantile, probs = 0.975), 
      apply(prob_outcome_du, 2, quantile, probs = 0.975), 
      apply(prob_outcome_df, 2, quantile, probs = 0.975), 
      apply(prob_outcome_f, 2, quantile, probs = 0.975)
    )
  )

# Control probabilities do not change through time
controls_hill <- data.frame(
  au = mean(prob_outcome_au_c),
  du = mean(prob_outcome_du_c),
  df = mean(prob_outcome_df_c),
  f = mean(prob_outcome_f_c)) |> 
  mutate(arm = 0) |> 
  pivot_longer(cols = c("au", "du", "df", "f")) |> 
  mutate( 
    lower = c( 
      quantile(prob_outcome_au_c, 0.025), 
      quantile(prob_outcome_du_c, 0.025), 
      quantile(prob_outcome_df_c, 0.025), 
      quantile(prob_outcome_f_c, 0.025) 
    ), 
    upper = c( 
      quantile(prob_outcome_au_c, 0.975), 
      quantile(prob_outcome_du_c, 0.975), 
      quantile(prob_outcome_df_c, 0.975), 
      quantile(prob_outcome_f_c, 0.975) 
      )
  )


# Get dummy data in plotable format
hut_outcomes_hill <- data.frame(stan_data_hill$y,
                           time = stan_data_hill$time,
                           arm = stan_data_hill$arm,
                           n = rowSums(stan_data_hill$y))

library(ggplot2)
## Plot trial outcomes and model fit
hut_outcomes_hill |> 
  tidyr::pivot_longer(cols = c(au, du, df, f)) |> 
  ggplot() +
  geom_point(aes(x = time, y = value/n, color = name), size = 0.5) +
  geom_line(data = outcome_prob_res_hill, aes(x = time, y = value, color = name)) +
  geom_hline(data = controls_hill, aes(yintercept = value, color = name)) +
  facet_grid(arm~name) +
  theme_bw() + theme(strip.background = element_rect(fill = "white"))


## Plot each probability through time
mos_prob_res_hill <- bind_rows(
  data.frame(value = colMeans(p_R), 
             lower = apply(p_R, 2, quantile, probs = 0.025), 
             upper = apply(p_R, 2, quantile, probs = 0.975), 
             time, name = "Repelled"), # R is repellence, K is killing before feeding and M is killing after feeding
  data.frame(value = colMeans(p_K), 
             lower = apply(p_K, 2, quantile, probs = 0.025), 
             upper = apply(p_K, 2, quantile, probs = 0.975), 
             time, name = "Killed before feeding"),
  data.frame(value = colMeans(p_M), 
             lower = apply(p_M, 2, quantile, probs = 0.025), 
             upper = apply(p_M, 2, quantile, probs = 0.975), 
             time, name = "Killed after feeding"),
  data.frame(value = colMeans(p_det), 
             lower = apply(p_det, 2, quantile, probs = 0.025), 
             upper = apply(p_det, 2, quantile, probs = 0.975), 
             time, name = "Deterred"))

mos_prob_res_hill |> 
  ggplot() +
  geom_line(aes(x = time, y = value, color = name)) +
  facet_grid(~name) +
  theme_bw() + theme(strip.background = element_rect(fill = "white")) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_ribbon( data = mos_prob_res_hill,
               aes( x = time,
                    ymin = lower,
                    ymax = upper,
                    fill = name ),
               alpha = 0.2, colour = NA )

