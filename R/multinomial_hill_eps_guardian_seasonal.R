# Seasonal version of multinomial_hill_eps_guardian_v2.R: fits
# multinomial_hill_eps_seasonal.stan to the Guardian hut-trial data. The constant
# biting probability pi is replaced by a seasonal pi(t), a Fourier series (on the
# logit scale) in hut-trial day -- no calendar dates needed, trial day 1 is the
# phase origin of the cycle. Everything else (eps-parametrised Hill decays for
# repellence, killing before/after feeding and deterrence) is as in v2.
#
# See R/multinomial_hill_eps_seasonal_testing.R for a simulate-and-recover check.

library(dplyr)
library(tidyr)

source("R/prepare_stan_data_seasonal.R")   # prepare_stan_data_seasonal(), seasonal_pi()

dtaMS = readxl::read_excel("C:/Users/esherrar/Documents/RProjects/spatial-emanators_1/data/Guardian_Sarah_2022.xlsx", sheet = "hut study")
dtaMS$date2 = lubridate::ymd(dtaMS$date)
dtaMS <- dtaMS[order(dtaMS$date2), ]   # ensure sorted by date
dtaMS$time_index <- as.integer(dtaMS$date2 - min(dtaMS$date2)) + 1

for(i in 1:nrow(dtaMS)){
  dtaMS$total[i] = dtaMS$`unfed alive 24`[i] + dtaMS$`fed alive 24`[i] +
    dtaMS$`unfed dead 24`[i] + dtaMS$`fed dead 24`[i]
}

dtaMS$Treatment = ifelse(dtaMS$treatment == "Control", 0, 1)
dtaMS <- dtaMS[order(dtaMS$Treatment,
                     dtaMS$time_index,
                     dtaMS$species), ]

dtaMS = subset(dtaMS, dtaMS$species == "arabiensis")


# dtaMS_summary <- dtaMS %>%
#   group_by(Treatment, time_index, species) %>%
#   summarise(
#     across(c("unfed alive 24",
#              "fed alive 24",
#              "unfed dead 24",
#              "fed dead 24",
#              total), sum, na.rm = TRUE)
#   )
# dtaMS_summary = subset(dtaMS_summary, dtaMS_summary$species == "arabiensis")
# 


real_hill = data.frame(HutTrial = rep(1, nrow(dtaMS)),
                       Hut = dtaMS$hut,
                       Treatment = dtaMS$Treatment,
                       Time = dtaMS$time_index,
                       au = dtaMS$`unfed alive 24`,
                       du = dtaMS$`unfed dead 24`,
                       df = dtaMS$`fed dead 24`,
                       f = dtaMS$`fed alive 24`)
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



# Seasonal settings: one cycle per 365 hut-trial days, two harmonics
# (allows an asymmetric / two-peaked season; set n_harmonics = 1 for a single sinusoid)
period      <- 365
n_harmonics <- 2
#filtered_data
# stan_data_hill <- prepare_stan_data_seasonal(real_hill, T_ref = 400,
#                                              period = period, n_harmonics = n_harmonics)
stan_data_hill <- prepare_stan_data_seasonal(summarised_data, T_ref = 400,
                                             period = period, n_harmonics = n_harmonics)


#
# ?install_cmdstan() 
library(cmdstanr)
mod_hill_eps <- cmdstan_model("stan_code/multinomial_hill_eps_seasonal.stan")
fit_hill_eps <- mod_hill_eps$sample(data = stan_data_hill, chains = 4, iter_warmup = 1000, iter_sampling = 1000, seed = 1, parallel_chains = 4)

## And interrogate the results
pars <- c("alpha0", "beta_pi", "pi_mean", "nu_u", "nu_f",
          "R0", "eps_R", "nR",
          "K0", "eps_K", "nK",
          "M0", "eps_M", "nM",
          "delta0", "eps_d", "nd")

fit_hill_eps$summary(variables = pars)

bayesplot::mcmc_trace(fit_hill_eps$draws(pars))

draws_hill_eps <- fit_hill_eps$draws(format = "df", variables = pars)


### What do these values mean in practice?
### Generate mosquito outcomes

time = 1:400
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

# Seasonal biting probability: draws x time matrix
p_pi <- seasonal_pi(draws_hill_eps$alpha0,
                    as.data.frame(draws_hill_eps)[, paste0("beta_pi[", seq_len(2 * n_harmonics), "]")],
                    time, period, n_harmonics)

p_A <- 1 - (1 - draws_hill_eps$nu_u) * (1 - p_K)
p_hF <- 1 - (1 - draws_hill_eps$nu_f) * (1 - p_M)

# Full model with SE impacts deterrence
prob_deterrence = p_det
one_minus_delta = (1 - p_det)
prop_control = 1/(1 + (1 - p_det))
prop_trt = 1 - prop_control

# Full model with SE impacts
prob_outcome_au = (1 - p_A) * (1 - p_hF) * ((1 - p_pi) + p_pi * p_R)
prob_outcome_du = p_A + (1 - p_A) * p_hF * ((1 - p_pi) + p_pi * p_R)
prob_outcome_df = p_pi * (1 - p_R) * (1 - p_A) * p_hF
prob_outcome_f = p_pi * (1 - p_R) * (1 - p_A) * (1 - p_hF)

plot(p_pi[1,] ~ c(1:400), type = "l",
col = "blue", ylim = c(0,1))
rr = sample(1:4000, 1000, replace = FALSE)
for(i in 1:1000){
  lines(p_pi[rr[i],] ~ c(1:400),
        col = adegenet::transp(col = "blue",
                               alpha = 0.2)
        )}
# Simplified model without SE impacts -- now varies through time with pi(t)
prob_outcome_au_c = (1 - draws_hill_eps$nu_u) * (1 - draws_hill_eps$nu_f) * (1 - p_pi)
prob_outcome_du_c = draws_hill_eps$nu_u + (1 - draws_hill_eps$nu_u) * draws_hill_eps$nu_f * (1 - p_pi)
prob_outcome_df_c = p_pi * (1 - draws_hill_eps$nu_u) * draws_hill_eps$nu_f
prob_outcome_f_c = p_pi * (1 - draws_hill_eps$nu_u) * (1 - draws_hill_eps$nu_f)

quantile(prob_outcome_au[,1], c(0.025,0.5,0.975))
quantile(prob_outcome_du[,1], c(0.025,0.5,0.975))
quantile(prob_outcome_df[,1], c(0.025,0.5,0.975))
quantile(prob_outcome_f[,1], c(0.025,0.5,0.975))

quantile(prob_outcome_au_c[,1], c(0.025,0.5,0.975))
quantile(prob_outcome_du_c[,1], c(0.025,0.5,0.975))
quantile(prob_outcome_df_c[,1], c(0.025,0.5,0.975))
quantile(prob_outcome_f_c[,1], c(0.025,0.5,0.975))

# Get data in plotable format
hut_outcomes_hill <- data.frame(stan_data_hill$y,
                                time = stan_data_hill$time,
                                arm = stan_data_hill$arm,
                                n = rowSums(stan_data_hill$y))

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

prob_in_hut_repel <- (1 - mu_K) * mu_R * (1-mu_M)
prob_in_hut_die_before <- (mu_K + (1-mu_K) * mu_R * mu_M)
prob_in_hut_feed <- 1 - prob_in_hut_repel - prob_in_hut_die_before
prob_in_hut_feed_die <- prob_in_hut_feed * mu_M
prob_in_hut_feed_survive <- prob_in_hut_feed * (1 - mu_M)

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
    prob_outcome_au_c,
    "repellence",
    "Placebo",
    0,
    time
  ),
  
  matrix_long(
    prob_outcome_du_c,
    "kill before",
    "Placebo",
    0,
    time
  ),
  
  matrix_long(
    prob_outcome_df_c,
    "kill after",
    "Placebo",
    0,
    time
  ),
  
  matrix_long(
    prob_outcome_f_c,
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
    fill = NA,#"white",
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
    breaks = c(1, 60, 120, 180, 240, 300, 360)
  ) +
  
  labs(
    x = NULL,
    y = "Outcomes (%)"
  ) +
  
  theme_classic(base_size = 11) +
  
  theme(
    # Rotate x-axis labels
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
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
      row = c( "Probabilities" = "gambiae" ), 
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
    breaks = c(1, 60, 120, 180, 240, 300, 360)
  ) +
  
  labs(
    x = "Time in days",
    y = "Probabilities (%)"
  ) +
  
  theme_classic(base_size = 11) +
  
  theme(
    # Rotate x-axis labels
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
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

# ------------------------------------------------------------
# Seasonal biting probability pi(t)
# ------------------------------------------------------------

pi_plot_data <- data.frame(
  time,
  estimate = apply(p_pi, 2, median),
  lower = apply(p_pi, 2, quantile, probs = 0.025),
  upper = apply(p_pi, 2, quantile, probs = 0.975)
)

# Observed control-arm fed fraction per night; its expectation is pi(t) * (1 - nu_u)
observed_control_fed <- hut_outcomes_hill |>
  filter(arm == 0) |>
  group_by(time) |>
  summarise(value = sum(df + f) / sum(n), n = sum(n), .groups = "drop")

p_biting <- ggplot(pi_plot_data, aes(x = time, y = estimate)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey40", alpha = 0.2) +
  geom_line(linewidth = 0.9, colour = "grey20") +
  geom_point(
    data = observed_control_fed,
    aes(x = time, y = value, size = n),
    shape = 21, fill = NA, colour = "black", stroke = 0.7
  ) +
  scale_size_continuous(range = c(1, 5)) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.25),
    labels = scales::label_percent(accuracy = 1),
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_x_continuous(
    breaks = c(1, 60, 120, 180, 240, 300, 360)
  ) +
  labs(
    x = "Time in days",
    y = "Biting probability (%)"
  ) +
  theme_classic(base_size = 11) +
  theme(
    panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.35),
    legend.position = "none"
  )

library(patchwork)

final_plot <- p_outcomes / p_mechanisms + #/ p_biting +
  plot_layout(
    heights = c(2, 1)
  ) +
  plot_annotation( tag_levels = "A" ) 

final_plot

draws <- fit_hill_eps$draws(variables = pars)

summary_results <- posterior::summarise_draws(
  draws,
  mean,
  sd,
  ~ quantile(.x, probs = 0.025),
  ~ quantile(.x, probs = 0.975)
)
