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

# simulate_hut_trial_data_hill <- function(
#     n_trials    = 1,
#     n_huts      = 6,          # kept even so each arm gets an equal share of huts
#     n_days      = 40,
#     seed        = 1,
#     rotate_huts = TRUE,       # TRUE: Latin-square-style rotation of which hut is "treatment" by day.
#     # FALSE: each hut keeps the same arm for the whole trial.
#     
#     # ---- true parameter values (edit these to explore different scenarios) ----
#     pi_true   = 0.65,         # P(would bite absent intervention)
#     nu_u_true = 0.05,         # background mortality hazard, before the average feeding time
#     nu_f_true = 0.03,         # background mortality hazard, after the average feeding time
#     
#     # Treatment-only hazards: each is a Hill/log-logistic decay f0 * thalf^n / (thalf^n + t^n)
#     R0_true = 0.50, thalf_R_true = 15, nR_true = 2.0,   # repelled before feeding
#     K0_true = 0.30, thalf_K_true = 20, nK_true = 2.0,   # killed before feeding
#     M0_true = 0.15, thalf_M_true = 25, nM_true = 1.5,   # killed after feeding
#     
#     rho_R = 0.05, rho_K = 0.05, rho_M = 0.05,           # within-hut overdispersion, ICC scale (treatment rows only)
#     
#     delta0_true = 0.30, thalf_d_true = 15, nd_true = 2.0,  # deterrence Hill decay
#     rho_delta   = 0.05,                                    # deterrence overdispersion, ICC scale
#     
#     baseline_catch_mean = 15  # mean nightly per-hut catch absent deterrence, before Poisson noise
# ) {
#   set.seed(seed)
#   
#   hill_decay <- function(t, f0, thalf, n) f0 * thalf^n / (thalf^n + t^n)
#   rbeta_mv   <- function(n_draw, mean, phi) rbeta(n_draw, mean * phi, (1 - mean) * phi)
#   
#   phi_R <- (1 - rho_R) / rho_R
#   phi_K <- (1 - rho_K) / rho_K
#   phi_M <- (1 - rho_M) / rho_M
#   
#   rows <- list()
#   
#   for (trial in seq_len(n_trials)) {
#     
#     fixed_arm <- rep(0:1, length.out = n_huts)
#     
#     for (day in seq_len(n_days)) {
#       
#       arm_this_day <- if (rotate_huts) {
#         (seq_len(n_huts) + day) %% 2
#       } else {
#         fixed_arm
#       }
#       
#       # Nightly mosquito-abundance multiplier shared by every hut that night.
#       night_multiplier <- rlnorm(1, meanlog = 0, sdlog = 0.3)
#       
#       # True deterrence at this trial-day (Hill decay), relative to the control reference.
#       delta_day <- hill_decay(day, delta0_true, thalf_d_true, nd_true)
#       
#       for (hut in seq_len(n_huts)) {
#         arm <- arm_this_day[hut]
#         
#         entry_rate <- baseline_catch_mean * night_multiplier *
#           if (arm == 1) (1 - delta_day) else 1
#         
#         N_j <- rpois(1, lambda = max(entry_rate, 0.1))
#         if (N_j == 0) next   # an empty catch just contributes no row for this hut-night
#         
#         if (arm == 1) {
#           meanR <- hill_decay(day, R0_true, thalf_R_true, nR_true)
#           meanK <- hill_decay(day, K0_true, thalf_K_true, nK_true)
#           meanM <- hill_decay(day, M0_true, thalf_M_true, nM_true)
#           
#           R <- rbeta_mv(1, meanR, phi_R)
#           K <- rbeta_mv(1, meanK, phi_K)
#           M <- rbeta_mv(1, meanM, phi_M)
#         } else {
#           R <- 0; K <- 0; M <- 0   # structural constraint: R = K = M = 0 in control
#         }
#         
#         A  <- 1 - (1 - nu_u_true) * (1 - K)   # combined hazard, before the average feeding time
#         hF <- 1 - (1 - nu_f_true) * (1 - M)   # combined hazard, after the average feeding time
#         
#         p <- c(
#           au = (1 - A) * (1 - hF) * ((1 - pi_true) + pi_true * R),   # alive, unfed
#           du = A + (1 - A) * hF * ((1 - pi_true) + pi_true * R),      # dead,  unfed
#           df = pi_true * (1 - R) * (1 - A) * hF,                     # dead,  fed
#           f  = pi_true * (1 - R) * (1 - A) * (1 - hF)                # alive, fed
#         )
#         p <- p / sum(p)   # guard against floating-point drift before rmultinom
#         
#         counts <- as.vector(rmultinom(1, N_j, p))
#         
#         rows[[length(rows) + 1]] <- data.frame(
#           HutTrial  = trial,
#           Hut       = hut,
#           Treatment = arm,
#           Time      = day,
#           au = counts[1], du = counts[2], df = counts[3], f = counts[4]
#         )
#       }
#     }
#   }
#   
#   bind_rows(rows)
# }
# 
# # ---- Quick-start: simulate, prepare, fit, and sanity-check recovery ----
# #
# raw_hill <- simulate_hut_trial_data_hill()
# write.csv(raw, "dummy_hut_trial_data_hill.csv", row.names = FALSE)

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
dtaMS = subset(dtaMS, dtaMS$species == "arabiensis")



# dtaMS_summary <- dtaMS %>%
#   group_by(treatment, time_index, species) %>%
#   summarise(
#     across(c("aliveunfed24",
#              "alivefed24",
#              "deadunfed24",
#              "deadfed24",
#              total), sum, na.rm = TRUE)
#   )
# dtaMS_summary = subset(dtaMS_summary, dtaMS_summary$species == "arabiensis")

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

real_hill = data.frame(HutTrial = rep(1, nrow(dtaMS)),
                       Hut = dtaMS$hut,
                       Treatment = dtaMS$treatment,
                       Time = dtaMS$time_index,
                       au = dtaMS$aliveunfed24,
                       du = dtaMS$deadunfed24,
                       df = dtaMS$deadfed24,
                       f = dtaMS$alivefed24)
real_hill$total = real_hill$au + real_hill$df + real_hill$du + real_hill$f
summarised_data <- real_hill %>%
  group_by(HutTrial, Hut, Treatment, Time) %>%
  summarise(
    au = sum(au, na.rm = TRUE),
    du = sum(du, na.rm = TRUE),
    df = sum(df, na.rm = TRUE),
    f = sum(f, na.rm = TRUE),
    total = sum(total, na.rm = TRUE),
    .groups = "drop"
  )

real_hillB = data.frame(HutTrial = c(rep(1, nrow(dtbMS1)),rep(2, nrow(dtbMS2))),
                        Hut = c(dtbMS1$hut,dtbMS2$hut),
                        Treatment = c(dtbMS1$treatment,dtbMS2$treatment),
                        Time = c(dtbMS1$time_index,dtbMS2$time_index),
                        au = c(dtbMS1$unf_live,dtbMS2$unf_live),
                        du = c(dtbMS1$unf_dead,dtbMS2$unf_dead),
                        df = c(dtbMS1$bf_dead,dtbMS2$bf_dead),
                        f = c( dtbMS1$bf_live,dtbMS2$bf_live)
)

# source("R/prepare_stan_data.R")     # unchanged -- same raw format, same aggregation logic
# stan_data_hill <- prepare_stan_data(raw_hill, T_ref = 35)
stan_data_hill <- prepare_stan_data(summarised_data, T_ref = 50)
stan_data_hill <- prepare_stan_data(real_hillB, T_ref = 50)
#
library(cmdstanr)
# mod_hill_eps <- cmdstan_model("stan_code/SE_multinomial_fits/multinomial_hill_eps_v2.stan") ##GOOD
mod_hill_eps <- cmdstan_model("stan_code/multinomial_hill_eps_v3.stan")
fit_hill_eps <- mod_hill_eps$sample(data = stan_data_hill, chains = 4, iter_warmup = 1000, iter_sampling = 1000, seed = 1, parallel_chains = 4)

## And interrogate the results
fit_hill_eps$summary(variables = c("pi", "nu_u", "nu_f",
                               "R0", "eps_R", "nR",
                               "K0", "eps_K", "nK",
                               "M0", "eps_M", "nM",
                               "delta0", "eps_d", "nd"))

bayesplot::mcmc_trace(fit_hill_eps$draws(c("pi", "nu_u", "nu_f",
                                             "R0", "eps_R", "nR",
                                             "K0", "eps_K", "nK",
                                             "M0", "eps_M", "nM",
                                             "delta0", "eps_d", "nd")))

draws_hill_eps <- fit_hill_eps$draws(format = "df", variables = c("pi", "nu_u", "nu_f",
                                                                "R0", "eps_R", "nR",
                                                                "K0", "eps_K", "nK",
                                                                "M0", "eps_M", "nM",
                                                                "delta0", "eps_d", "nd"))


### What do these values mean in practice?
### Generate mosquito outcomes

time = 1:40
trt = c(0,1)

hill_eps <- function(f0, T_ref, n, t, eps){
  rho <- 1/eps - 1
  return(f0 / (1 + rho * (t / T_ref)^n))
}

p_R = p_K = p_M = p_det = array(dim = c(nrow(draws_hill_eps),length(time)))

for(t in 1:length(time)){
  for(i in 1:nrow(draws_hill_eps)){
    p_R[i,t] <- hill_eps(draws_hill_eps$R0[i],stan_data_hill$T_ref[1], draws_hill_eps$nR[i], t, draws_hill_eps$eps_R[i])
    p_K[i,t] <- hill_eps(draws_hill_eps$K0[i], stan_data_hill$T_ref[1], draws_hill_eps$nK[i], t, draws_hill_eps$eps_K[i])
    p_M[i,t] <- hill_eps(draws_hill_eps$M0[i], stan_data_hill$T_ref[1], draws_hill_eps$nM[i], t, draws_hill_eps$eps_M[i])
    p_det[i,t] <- hill_eps(draws_hill_eps$delta0[i], stan_data_hill$T_ref[1], draws_hill_eps$nd[i], t, draws_hill_eps$eps_d[i])
  }  
}

p_A <- 1 - (1 - draws_hill_eps$nu_u) * (1 - p_K)
p_hF <- 1 - (1 - draws_hill_eps$nu_f) * (1 - p_M)

# Full model with SE impacts deterrence
prob_deterrence = p_det
one_minus_delta = (1 - p_det)
prop_control = 1/(1 + (1 - p_det))
prop_trt = 1 - prop_control

# Full model with SE impacts
prob_outcome_au = (1 - p_A) * (1 - p_hF) * ((1 - draws_hill_eps$pi) + draws_hill_eps$pi * p_R)
prob_outcome_du = p_A + (1 - p_A) * p_hF * ((1 - draws_hill_eps$pi) + draws_hill_eps$pi * p_R)
prob_outcome_df = draws_hill_eps$pi * (1 - p_R) * (1 - p_A) * p_hF
prob_outcome_f = draws_hill_eps$pi * (1 - p_R) * (1 - p_A) * (1 - p_hF)

# Simplified model without SE impacts
prob_outcome_au_c = (1 - draws_hill_eps$nu_u) * (1 - draws_hill_eps$nu_f) * (1 - draws_hill_eps$pi)
prob_outcome_du_c = draws_hill_eps$nu_u + (1 - draws_hill_eps$nu_u) * draws_hill_eps$nu_f * (1 - draws_hill_eps$pi)
prob_outcome_df_c = draws_hill_eps$pi * (1 - draws_hill_eps$nu_u) * draws_hill_eps$nu_f
prob_outcome_f_c = draws_hill_eps$pi * (1 - draws_hill_eps$nu_u) * (1 - draws_hill_eps$nu_f)

quantile(prob_outcome_au[,1], c(0.025,0.5,0.975))
quantile(prob_outcome_du[,1], c(0.025,0.5,0.975))
quantile(prob_outcome_df[,1], c(0.025,0.5,0.975))
quantile(prob_outcome_f[,1], c(0.025,0.5,0.975))

quantile(prob_outcome_au_c, c(0.025,0.5,0.975))
quantile(prob_outcome_du_c, c(0.025,0.5,0.975))
quantile(prob_outcome_df_c, c(0.025,0.5,0.975))
quantile(prob_outcome_f_c, c(0.025,0.5,0.975))
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

times <- sort(unique(hut_outcomes_hill$time))

posterior_hill <- purrr::map_dfr(
  c("au", "du", "df", "f"),
  \(outcome) {
    mat <- get(paste0("prob_outcome_", outcome))
    
    tibble::as_tibble(mat) |>
      dplyr::mutate(draw = dplyr::row_number()) |>
      tidyr::pivot_longer(
        cols = -draw,
        names_to = "time_idx",
        values_to = "value"
      ) |>
      dplyr::mutate(
        name = outcome,
        time = times[as.integer(gsub("V", "", time_idx))],
        arm = 1
      )
  }
)

# ## Plot trial outcomes and model fit
# library(ggplot2)
# top_panel = hut_outcomes_hill |>
#   tidyr::pivot_longer(cols = c(au, du, df, f)) |>
#   # dplyr::mutate(
#   #   name = factor(name, levels = c("du", "df", "au", "f"))
#   # ) |>
#   ggplot() +
#   
#   # posterior draws — treatment arm only
#   geom_line(
#     data = posterior_hill,
#     aes(
#       x = time,
#       y = value,
#       group = draw,
#       color = name
#     ),
#     alpha = 0.015,
#     linewidth = 0.2
#   ) +
#   
#   # observed outcomes
#   geom_point(
#     aes(x = time, y = value/n, color = name, size = n),
#     alpha = 0.5
#   ) +
#   scale_size_continuous(range = c(1, 5)) +
#   # posterior mean
#   geom_line(
#     data = outcome_prob_res_hill,
#     aes(x = time, y = value, color = name),
#     linewidth = 0.8
#   ) +
#   
#   geom_hline(
#     data = controls_hill,
#     aes(yintercept = value, color = name)
#   ) +
#     scale_y_continuous(
#       limits = c(0, 1),
#       breaks = seq(0, 1, 0.2),
#       labels = function(x) x * 100,
#       expand = expansion(mult = c(0.02, 0.05))
#     ) +
#   facet_grid(
#     arm ~ name,
#     labeller = labeller(
#       arm = c(
#         `0` = "Placebo",
#         `1` = "Emanator"
#       ),
#       name = c(
#         au = "Alive unfed",
#         du = "Dead unfed",
#         df = "Dead fed",
#         f  = "Alive fed"
#       )
#     )
#   ) +
#   
#   labs(x = "Time in days",
#        y = "Outcome (%)") +
#   
#   theme_bw() +
#   theme(legend.position = "none",
#         strip.background = element_rect(fill = "white"))
# 
# 
# ## Plot each probability through time
# mos_prob_res_hill <- bind_rows(
#   data.frame(value = colMeans(p_R), 
#              lower = apply(p_R, 2, quantile, probs = 0.025), 
#              upper = apply(p_R, 2, quantile, probs = 0.975), 
#              time, name = "Repelled"), # R is repellence, K is killing before feeding and M is killing after feeding
#   data.frame(value = colMeans(p_K), 
#              lower = apply(p_K, 2, quantile, probs = 0.025), 
#              upper = apply(p_K, 2, quantile, probs = 0.975), 
#              time, name = "Killed before feeding"),
#   data.frame(value = colMeans(p_M), 
#              lower = apply(p_M, 2, quantile, probs = 0.025), 
#              upper = apply(p_M, 2, quantile, probs = 0.975), 
#              time, name = "Killed after feeding"),
#   data.frame(value = colMeans(p_det), 
#              lower = apply(p_det, 2, quantile, probs = 0.025), 
#              upper = apply(p_det, 2, quantile, probs = 0.975), 
#              time, name = "Deterred"))
# 
# second_panel = mos_prob_res_hill |> 
#   ggplot() +
#   geom_line(aes(x = time, y = value, color = name)) +
#   facet_grid(~name) +
#   geom_hline(yintercept = 0, linetype = "dashed") +
#   scale_y_continuous(
#     limits = c(0, 1),
#     breaks = seq(0, 1, 0.2),
#     labels = function(x) x * 100,
#     expand = expansion(mult = c(0.02, 0.05))
#   ) +
#   labs(x = "Time in days",
#        y = "Outcome (%)") +
#   geom_ribbon( data = mos_prob_res_hill,
#                aes( x = time,
#                     ymin = lower,
#                     ymax = upper,
#                     fill = name ),
#                alpha = 0.2, colour = NA ) +
#   theme_bw() + theme(legend.position = "none",
#                      strip.background = element_rect(fill = "white"))
# 
# library(patchwork)
# top_panel / second_panel + 
#   plot_layout( heights = c(3, 2) ) + 
#   plot_annotation( tag_levels = "A" ) 
# 
# 
# ## Expectations
# mos_prob_res_outcomes_hill <- bind_rows(
#   data.frame(value = colMeans(p_det) + (1 - colMeans(p_det)) * colMeans(p_R), time, name = "Repelled"),
#   data.frame(value = (1 - colMeans(p_det)) * colMeans(p_K), time, name = "Dead prior to biting"),
#   data.frame(value = (1 - colMeans(p_det)) * (1 - colMeans(p_R) - colMeans(p_K)) * colMeans(p_M), time, name = "Dead after biting"),
#   data.frame(value = (1 - colMeans(p_det)) * (1 - colMeans(p_R) - colMeans(p_K)) * (1 - colMeans(p_M)), time, name = "Bites and survives"))
# 
# benin = mos_prob_res_outcomes_hill |> 
#   ggplot() +
#   geom_line(aes(x = time, y = value, color = name)) +
#   geom_hline(
#     yintercept = 0.01,
#     linetype = "dashed"
#   ) +
#   facet_grid(~name) +
#   theme_bw() + theme(strip.background = element_rect(fill = "white")) +
#   geom_hline(yintercept = 0, linetype = "solid")
# 
# library(patchwork)
# tanzania / benin + 
#   plot_layout( heights = c(2, 2) ) + 
#   plot_annotation( tag_levels = "A" ) 

library(ggplot2)
# ------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------

mechanism_cols <- c(
  "deterrence"  = "#4C78A8",
  "repellence"  = "#59A14F",
  "kill before" = "#E15759",
  "kill after"  = "#F28E2B",
  "survive"     = "#B07AA1"
)

matrix_long <- function(mat, mechanism, row_name, arm, time) {
  
  tibble::as_tibble(mat) |>
    dplyr::mutate(draw = dplyr::row_number()) |>
    tidyr::pivot_longer(
      cols = -draw,
      names_to = "time_idx",
      values_to = "value"
    ) |>
    dplyr::mutate(
      time = time[as.integer(gsub("V", "", time_idx))],
      row = row_name,
      mechanism = mechanism,
      arm = arm
    ) |>
    dplyr::select(
      draw, time, row, mechanism, arm, value
    )
}


vector_to_time_matrix <- function(x, n_time) {
  matrix(
    rep(x, times = n_time),
    nrow = length(x),
    ncol = n_time
  )
}


n_draws <- nrow(p_det)
n_time  <- ncol(p_det)


# ------------------------------------------------------------
# Derived SE probabilities
# ------------------------------------------------------------

mu_det <- p_det
mu_R   <- p_R
mu_K   <- p_K
mu_M   <- p_M

rSe = mu_det + (1 - mu_det) * (1 - mu_K) * mu_R * (1-mu_M) #: mosquitos must survive in order to be counted as repeating.
dSe = (1 - mu_det) * (mu_K + (1-mu_K) * mu_R * mu_M) #: mosquitos that are repelled can also die without feeding
bSe = 1 - rSe - dSe #

## Given a mosquito has bitten:
sSe = 1 - mu_M 

# rSE <- mu_det + (1 - mu_det) * mu_R 
# dSE <- (1 - mu_det) * mu_K 
# bSE <- 1 - rSE - dSE 

prob_kill_after <- bSe * mu_M 
prob_survive <- bSe * sSe

prob_deterred_from_entering = mu_det
prob_in_hut_repel <- (1 - mu_K) * mu_R * (1-mu_M)
prob_in_hut_die_before <- (mu_K + (1-mu_K) * mu_R * mu_M)
prob_in_hut_feed <- 1 - prob_in_hut_repel - prob_in_hut_die_before
prob_in_hut_feed_die <- prob_in_hut_feed * mu_M
prob_in_hut_feed_survive <- prob_in_hut_feed * (1 - mu_M)

summaries_of_effects_f = function(x){
  d1 = mean(x)
  d2 = as.numeric(quantile(x, c(0.025,0.975)))
  return(c(d1,d2))
}

summaries_of_effects_f(prob_deterred_from_entering[,1])
summaries_of_effects_f(prob_in_hut_die_before[,1])
summaries_of_effects_f(prob_in_hut_repel[,1])
summaries_of_effects_f(prob_in_hut_feed_survive[,1])


# ------------------------------------------------------------
# PLACEBO
# ------------------------------------------------------------

placebo_plot <- dplyr::bind_rows(
  
  matrix_long(
    prop_control,
    "deterrence",
    "Placebo",
    0,
    time
  ),
  
  matrix_long(
    vector_to_time_matrix(prob_outcome_au_c, n_time),
    "repellence",
    "Placebo",
    0,
    time
  ),
  
  matrix_long(
    vector_to_time_matrix(prob_outcome_du_c, n_time),
    "kill before",
    "Placebo",
    0,
    time
  ),
  
  matrix_long(
    vector_to_time_matrix(prob_outcome_df_c, n_time),
    "kill after",
    "Placebo",
    0,
    time
  ),
  
  matrix_long(
    vector_to_time_matrix(prob_outcome_f_c, n_time),
    "survive",
    "Placebo",
    0,
    time
  )
)


# ------------------------------------------------------------
# EMANATOR
# ------------------------------------------------------------

emanator_plot <- dplyr::bind_rows(
  
  matrix_long(
    prop_trt,
    "deterrence",
    "Emanator",
    1,
    time
  ),
  
  matrix_long(
    prob_outcome_au,
    "repellence",
    "Emanator",
    1,
    time
  ),
  
  matrix_long(
    prob_outcome_du,
    "kill before",
    "Emanator",
    1,
    time
  ),
  
  matrix_long(
    prob_outcome_df,
    "kill after",
    "Emanator",
    1,
    time
  ),
  
  matrix_long(
    prob_outcome_f,
    "survive",
    "Emanator",
    1,
    time
  )
)


# ------------------------------------------------------------
# MECHANISTIC PROBABILITIES
# ------------------------------------------------------------

probabilities_plot <- dplyr::bind_rows(
  
  matrix_long(
    mu_det,
    "deterrence",
    "Probabilities",
    NA,
    time
  ),
  
  matrix_long(
    prob_in_hut_repel,
    "repellence",
    "Probabilities",
    NA,
    time
  ),
  
  matrix_long(
    prob_in_hut_die_before,
    "kill before",
    "Probabilities",
    NA,
    time
  ),
  
  matrix_long(
    prob_in_hut_feed_die,
    "kill after",
    "Probabilities",
    NA,
    time
  ),
  
  matrix_long(
    prob_in_hut_feed_survive,
    "survive",
    "Probabilities",
    NA,
    time
  )
)

quantile(mu_det[,1], c(0.025,0.5,0.975))
mean(mu_det[,1])

quantile(mu_R[,1], c(0.025,0.5,0.975))
mean(mu_R[,1])

quantile(mu_K[,1], c(0.025,0.5,0.975))
mean(mu_K[,1])

quantile(prob_survive[,1], c(0.025,0.5,0.975))
mean(prob_survive[,1])


# ============================================================
# Combine all posterior draws
# ============================================================

plot_draws <- bind_rows(
  placebo_plot,
  emanator_plot,
  probabilities_plot
) |>
  mutate(
    row = factor(
      row,
      levels = c(
        "Placebo",
        "Emanator",
        "Probabilities"
      )
    ),
    mechanism = factor(
      mechanism,
      levels = c(
        "deterrence",
        "repellence",
        "kill before",
        "kill after",
        "survive"
      )
    )
  )


# ============================================================
# Summarise posterior draws
# ============================================================

plot_data <- plot_draws |>
  group_by(
    row,
    mechanism,
    time
  ) |>
  summarise(
    estimate = median(value, na.rm = TRUE),
    lower = quantile(
      value,
      probs = 0.025,
      na.rm = TRUE
    ),
    upper = quantile(
      value,
      probs = 0.975,
      na.rm = TRUE
    ),
    .groups = "drop"
  )


observed_outcomes <- hut_outcomes_hill |>
  group_by(time, arm) |>
  summarise(
    total = sum(n),
    repellence = sum(au),
    `kill before` = sum(du),
    `kill after` = sum(df),
    survive = sum(f),
    .groups = "drop"
  ) |>
  pivot_longer(
    cols = c(
      repellence,
      `kill before`,
      `kill after`,
      survive
    ),
    names_to = "mechanism",
    values_to = "count"
  ) |>
  mutate(
    row = if_else(
      arm == 0,
      "Placebo",
      "Emanator"
    ),
    value = count / total,
    n = total
  ) |>
  select(
    time,
    arm,
    row,
    mechanism,
    value,
    n
  )


observed_deterrence <- hut_outcomes_hill |>
  group_by(time, arm) |>
  summarise(
    n = sum(n),
    .groups = "drop"
  ) |>
  group_by(time) |>
  mutate(
    value = n / sum(n)
  ) |>
  ungroup() |>
  mutate(
    row = if_else(
      arm == 0,
      "Placebo",
      "Emanator"
    ),
    mechanism = "deterrence"
  ) |>
  select(
    time,
    arm,
    n,
    row,
    mechanism,
    value
  )


observed_plot <- bind_rows(
  observed_deterrence,
  observed_outcomes
) |>
  mutate(
    row = factor(
      row,
      levels = c(
        "Placebo",
        "Emanator"
      )
    ),
    mechanism = factor(
      mechanism,
      levels = c(
        "deterrence",
        "repellence",
        "kill before",
        "kill after",
        "survive"
      )
    )
  )


outcome_labels <- c(
  "deterrence"  = "Location",
  "repellence"  = "Alive and unfed",
  "kill before" = "Dead and unfed",
  "kill after"  = "Dead and fed",
  "survive"     = "Alive and fed"
)

p_outcomes <- ggplot(
  plot_data |>
    filter(row %in% c("Placebo", "Emanator")),
  aes(
    x = time,
    y = estimate,
    colour = mechanism,
    fill = mechanism
  )
) +
  
  geom_ribbon(
    aes(ymin = lower, ymax = upper),
    alpha = 0.18,
    colour = NA
  ) +
  
  geom_line(
    linewidth = 0.9
  ) +
  
  geom_point(
    data = observed_plot,
    aes(x = time, y = value, size = n),
    shape = 21,
    fill = "white",
    colour = "black",
    stroke = 0.7,
    inherit.aes = FALSE
  ) +
  scale_size_continuous(
    range = c(1, 5)
  )+
  
  facet_grid(
    row ~ mechanism,
    labeller = labeller(
      mechanism = outcome_labels
    )
  ) +
  
  scale_colour_manual(values = mechanism_cols) +
  scale_fill_manual(values = mechanism_cols) +
  
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.25),
    labels = scales::label_percent(accuracy = 1),
    expand = expansion(mult = c(0, 0.02))
  ) +
  
  scale_x_continuous(
    breaks = c(1, 10, 20, 30, 40)
  ) +
  
  labs(
    x = NULL,
    y = "Outcomes (%)"
  ) +
  
  theme_classic(base_size = 11) +
  
  theme(
    strip.background = element_rect(
      fill = "grey95",
      colour = NA
    ),
    strip.text.x = element_text(
      face = "bold",
      size = 9
    ),
    strip.text.y = element_text(
      face = "bold"
    ),
    panel.spacing = unit(0.7, "lines"),
    panel.grid.major.y = element_line(
      colour = "grey90",
      linewidth = 0.35
    ),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )

mechanism_labels <- c(
  "deterrence"  = "Deterrence",
  "repellence"  = "Repellence",
  "kill before" = "Dead before feeding",
  "kill after"  = "Dead after feeding",
  "survive"     = "Feeds and survives"
)

p_mechanisms <- ggplot(
  plot_data |>
    filter(row == "Probabilities"),
  aes(
    x = time,
    y = estimate,
    colour = mechanism,
    fill = mechanism
  )
) +
  
  geom_ribbon(
    aes(ymin = lower, ymax = upper),
    alpha = 0.18,
    colour = NA
  ) +
  
  geom_line(
    linewidth = 0.9
  ) +
  
  facet_grid(
    row ~ mechanism, 
    labeller = labeller( 
      row = c( "Probabilities" = "arabiensis" ), 
      mechanism = mechanism_labels 
      ) 
    ) +
  
  scale_colour_manual(values = mechanism_cols) +
  scale_fill_manual(values = mechanism_cols) +
  
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.25),
    labels = scales::label_percent(accuracy = 1),
    expand = expansion(mult = c(0, 0.02))
  ) +
  
  scale_x_continuous(
    breaks = c(1, 10, 20, 30, 40)
  ) +
  
  labs(
    x = "Time in days",
    y = "Probabilities (%)"
  ) +
  
  theme_classic(base_size = 11) +
  
  theme(
    strip.background = element_rect(
      fill = "grey90",
      colour = NA
    ),
    strip.text.x = element_text(
      face = "bold",
      size = 10
    ),
    strip.text.y = element_text(
      face = "bold"
    ),
    panel.spacing = unit(0.7, "lines"),
    panel.grid.major.y = element_line(
      colour = "grey90",
      linewidth = 0.35
    ),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )

library(patchwork)

final_plot <- p_outcomes / p_mechanisms +
  plot_layout(
    heights = c(2, 1)
  ) + 
  plot_annotation( tag_levels = "A" ) 

final_plot

summary_results <- posterior::summarise_draws(
  draws_hill_eps,
  mean,
  median,
  sd,
  ~ quantile(.x, probs = 0.025),
  ~ quantile(.x, probs = 0.975)
)
