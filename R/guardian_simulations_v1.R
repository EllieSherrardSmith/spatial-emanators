devtools::load_all()
devtools::document()

library("malariaEquilibrium")
library("malariasimulation")
library(cali)

####################################################
##
## PARAMETER SENSITIVITY ANALYSIS 
##
####################################################

## parameters taken from malariaverse
## demography
## treatment
## seasonality
KEN_df <- site::fetch_site(iso3c = "KEN", version = "malariaverse_01_2025") # Done
scene::plot_interventions(
  interventions = KEN_df$interventions, 
  population =  KEN_df$population$population_total, 
  group_var = "name_1", 
  include = c("itn_use", "tx_cov", "irs_cov", "rtss_cov", "r21_cov", "smc_cov"), 
  labels = c("ITN usage", "Treatment", "IRS", "RTSS", "R21", "SMC")
)

which(unique(KEN_df$sites$name_1) == "Busia")
EIR_list = array(dim = c(1,2))
site_filter = list()

site_filter[[1]] <- KEN_df$sites[which(unique(KEN_df$sites$name_1) == "Busia"),]

# Subset to site 
fit_site <- site::subset_site(
  site = KEN_df,
  site_filter = site_filter[[1]]
)

# Convert to malariasimulation parameters
fit_site_parameters <- site::site_parameters(
  interventions = fit_site$interventions,
  demography = fit_site$demography,
  vectors = fit_site$vectors$vector_species,
  seasonality = fit_site$seasonality$seasonality_parameters
)
  
fit_site_parameters$prevalence_rendering_min_ages = c( 1, 0,  5,  15,   0) * 365 ## Prev in 6 months to 5 years measured, 
fit_site_parameters$prevalence_rendering_max_ages = c(10, 5, 15, 100, 100) * 365
  
fit_site_parameters$incidence_rendering_min_ages = c(1, 0,  5,  15,   0) * 365 ## as prev
fit_site_parameters$incidence_rendering_max_ages = c(10,5, 15, 100, 100) * 365
  
fit_site_parameters$clinical_incidence_rendering_min_ages = c(1, 0,  5,  15,   0) * 365 ## as prev
fit_site_parameters$clinical_incidence_rendering_max_ages = c(10,5, 15, 100, 100) * 365

fit_site_parameters$severe_incidence_rendering_min_ages = c(1,  0,  5,  15,   0) * 365 ## as prev
fit_site_parameters$severe_incidence_rendering_max_ages = c(10, 5, 15, 100, 100) * 365

  ## From MAP data 31 July 2025
  ## ** Some mismatches exist in the names - correct these first
map_Data <- read.csv("C:/Users/esherrar/Documents/RProjects/spatial-emanators_1/simulating projections for Unitaid/Subnational Unit-data MAP 31 July 2025.csv", header = TRUE)
sort(unique(KEN_df$sites$name_1))
sort(unique(map_Data$Name[map_Data$ISO3 == "KEN"]))

prev_year = 2010:2022 ## 
checker_year = 365*unique(ceiling((3651:8030)/365))
prevalence_2_10yrs_MAP = map_Data$Value[map_Data$Units == "per 100 Children" &
                                            map_Data$Name == "Busia"]/100

## From population data
popn_data <- KEN_df$population$population_total
population_trend_Estimate = popn_data$pop[popn_data$name_1 == "Busia"][16:25]

##--- Setting up equilibrium

year <- 365 
years <- 25 # outputs wanted are 2015-2024, simulate from 2010? 
sim_length <- year * years

# Prepare a summary function that returns the mean PfPR2-10 from each simulation output: 
# annual_pfpr_summary <- function(x){
#   #set the years to target
#   year <- ceiling(x$timestep[3651:8030] / 365)
#   #estimate the Pfpr
#   pfpr <- x$n_detect_lm_730_3650[3651:8030]  / x$n_age_730_3650[3651:8030]
#   ##tabulise
#   tapply(pfpr, year, mean)
# }
# 
# # Establish a target PfPR2-10 value:
# target_pfpr <- prevalence_2_10yrs_MAP

# Add a parameter to the parameter list specifying the number of timesteps 
# fit_site_parameters$timesteps <- sim_length

# Run the calibrate() function:
# set.seed(123)
# out <- calibrate(
#   parameters = fit_site_parameters,
#   target = target_pfpr[2:13],
#   summary_function = annual_pfpr_summary,
#   eq_prevalence = target_pfpr[1]
# )
## out - ESTIMATED AT 34.22956

# fit_site_parameters$pev = TRUE
# fit_site_parameters$pev_epi_coverages = rep(0,25)
# fit_site_parameters$spraying = FALSE
# 
# fit_site_parameters$human_population <- 30000
# fit_site_parameters <- set_equilibrium(fit_site_parameters, init_EIR = 34.22956)
# 
# output0 <- run_simulation(timesteps = sim_length, parameters = fit_site_parameters)  
# ## Generate prevalence
# output0$pfpr <- output0$n_detect_lm_0_36500  / output0$n_age_0_36500
# 
# # store[[1]] = output0
# # params_store[[1]] = fit_site_parameters
# # fit_site_parameters$total_M
# # 
# # saveRDS(store[[1]],paste0("C:/Users/esherrar/Documents/RProjects/spatial-emanators/CF_Sims_store_Busia_KENtrial"))
# # saveRDS(params_store[[1]],paste0("C:/Users/esherrar/Documents/RProjects/spatial-emanators/Base_parameters_store_Busia_KENtrial"))
# 
# # cases_pp_all_age_cf = array(dim = c(2,4))
# # people_age_cohort_cf = array(dim = c(2,7))
# # for(i in 1:2){
# #   ## all-age
# #   cases_pp_all_age_cf[i,1] = sum(output0$n_inc_clinical_0_36500[((20+i)*365+315):((22+i)*365-1+315)])/
# #     mean(output0$n_age_0_36500[((20+i)*365+315):((22+i)*365-1+315)])
# #   people_age_cohort_cf[i,1] = mean(output0$n_age_0_36500[((20+i)*365+315):((22+i)*365-1+315)])
# #   ## U5
# #   cases_pp_all_age_cf[i,2] = sum(output0$n_inc_clinical_0_1825[((20+i)*365+315):((22+i)*365-1+315)])/
# #     mean(output0$n_age_0_1825[((20+i)*365+315):((22+i)*365-1+315)])
# #   people_age_cohort_cf[i,2] = mean(output0$n_age_0_1825[((20+i)*365+315):((22+i)*365-1+315)])
# #   ## 5-15
# #   cases_pp_all_age_cf[i,3] = sum(output0$n_inc_clinical_1825_5475[((20+i)*365+315):((22+i)*365-1+315)])/
# #     mean(output0$n_age_1825_5475[((20+i)*365+315):((22+i)*365-1+315)])
# #   people_age_cohort_cf[i,3] = mean(output0$n_age_1825_5475[((20+i)*365+315):((22+i)*365-1+315)])
# #   ## >15
# #   cases_pp_all_age_cf[i,4] = sum(output0$n_inc_clinical_5475_36500[((20+i)*365+315):((22+i)*365-1+315)])/
# #     mean(output0$n_age_5475_36500[((20+i)*365+315):((22+i)*365-1+315)])
# #   people_age_cohort_cf[i,4] = mean(output0$n_age_5475_36500[((20+i)*365+315):((22+i)*365-1+315)])
# # }
# # ## percentage of people in each age cohort
# # people_age_cohort_cf[,5] = people_age_cohort_cf[,2]/people_age_cohort_cf[,1]
# # people_age_cohort_cf[,6] = people_age_cohort_cf[,3]/people_age_cohort_cf[,1]
# # people_age_cohort_cf[,7] = people_age_cohort_cf[,4]/people_age_cohort_cf[,1]
# 
# par(mfrow = c(1,1))
# par(mar = c(4,4,1,1))
# plot(output0$pfpr ~ output0$timestep,
#      type = "l",
#      xlim = c(20*365,25*365),xaxt="n",
#      yaxt = "n",ylim = c(0,1.15),
#      xlab = "Year",ylab = "",
#      bty="n")
# mtext("PfPR all age LM (%)", side = 2, line = 3)
# mtext("Intervention use (%)", side = 2, line = 2)
# # points(prevalence_2_10yrs_MAP[2:13] ~ checker_year, pch = 19, col ="red")
# axis(1,at = seq(1,25*365,365),labels = 2000:2024)
# axis(2, las = 2, at = seq(0,0.8,0.2),
#      labels = seq(0,80,20))
# # abline(v = 21*365)
# # polygon(c(output0$timestep[(14*365):(25*365)],
# #           rev(output1$timestep[(14*365):(25*365)])),
# #         c(output1$pfpr[(14*365):(25*365)],
# #           rev(output0$pfpr[(14*365):(25*365)])),
# #         col = adegenet::transp("blue",0.2),border=NA)
# 
# # text(365*20,1.1,"Busia County", cex=1.1)
# add_itn = KEN_df$interventions$itn_use[KEN_df$interventions$name_1 == "Busia"]
# add_smc = KEN_df$interventions$smc_cov[KEN_df$interventions$name_1 == "Busia"]
# add_trt = KEN_df$interventions$tx_cov[KEN_df$interventions$name_1 == "Busia"]
# add_trt_act = add_trt*KEN_df$interventions$prop_act[KEN_df$interventions$name_1 == "Busia"]
# annual_year = 365*unique(ceiling((1:(25*365))/365))
# 
# lines(add_itn ~
#         annual_year, col = "purple",lwd=2)
# lines(add_smc ~
#         annual_year, col = "#CC38A1",lwd=2)
# lines(add_trt ~
#         annual_year, col = "aquamarine2",lwd=2)
# lines(add_trt_act ~
#         annual_year, col = "aquamarine2",lwd=2, lty = 2)
# legend("topleft", legend = c("ITN","SMC","Treatment","(with ACT)"),
#        lty = c(1,1,1,2),
#        lwd = 2,
#        col = c("purple","#CC38A1","aquamarine2","aquamarine2"),
#        bty="n")
# 
# ## Now figure out how to adjust the parameter list to include spatial emanators,
# ## and remove anything / change anything to reflect trial better... 

##############################
##
## read in the extracted parameters from the statistical model
library(dplyr)


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
    pi_prior_a = 6,# 0.699 * 6,
    pi_prior_b = 4,# (1 - 0.699) * 6,
    
    # only consumed by the eps parametrised Hill decay models
    T_ref = T_ref
    
    
  )
}

# source("R/prepare_stan_data.R")     # unchanged -- same raw format, same aggregation logic
# stan_data_hill <- prepare_stan_data(raw_hill, T_ref = 35)
# stan_data_hilla <- prepare_stan_data(summarised_data, T_ref = 400)

#
library(cmdstanr)
# mod_hill_eps <- cmdstan_model("C:/Users/esherrar/Documents/RProjects/spatial-emanators/stan_code/SE_multinomial_fits/multinomial_hill_eps_v2.stan")
# fit_hill_epsAa <- mod_hill_eps$sample(data = stan_data_hilla, chains = 4, iter_warmup = 1000, iter_sampling = 1000, seed = 1, parallel_chains = 4)

# ## And interrogate the results
# fit_hill_epsAa$summary(variables = c("pi", "nu_u", "nu_f",
#                                    "R0", "eps_R", "nR",
#                                    "K0", "eps_K", "nK",
#                                    "M0", "eps_M", "nM",
#                                    "delta0", "eps_d", "nd"))
# 
# bayesplot::mcmc_trace(fit_hill_epsAa$draws(c("pi", "nu_u", "nu_f",
#                                            "R0", "eps_R", "nR",
#                                            "K0", "eps_K", "nK",
#                                            "M0", "eps_M", "nM",
#                                            "delta0", "eps_d", "nd")))
# 
# 
# draws_hill_eps_aa <- fit_hill_epsAa$draws(format = "df", variables = c("pi", "nu_u", "nu_f",
#                                                                   "R0", "eps_R", "nR",
#                                                                   "K0", "eps_K", "nK",
#                                                                   "M0", "eps_M", "nM",
#                                                                   "delta0", "eps_d", "nd"))
# 
# 
##############################
stan_data_hilla <- prepare_stan_data_seasonal(summarised_data, T_ref = 400,
                                             period = period, n_harmonics = n_harmonics)


#
# ?install_cmdstan() 
library(cmdstanr)
mod_hill_eps <- cmdstan_model("stan_code/multinomial_hill_eps_seasonal.stan")
fit_hill_epsAa <- mod_hill_eps$sample(data = stan_data_hilla, chains = 4, iter_warmup = 1000, iter_sampling = 1000, seed = 1, parallel_chains = 4)

## And interrogate the results
pars <- c("alpha0", "beta_pi", "pi_mean", "nu_u", "nu_f",
          "R0", "eps_R", "nR",
          "K0", "eps_K", "nK",
          "M0", "eps_M", "nM",
          "delta0", "eps_d", "nd")

fit_hill_epsAa$summary(variables = pars)

bayesplot::mcmc_trace(fit_hill_epsAa$draws(pars))

draws_hill_eps_aa <- fit_hill_epsAa$draws(format = "df", variables = pars)



dtsens = data.frame(eir = rep(34.22956,3),
                 Q0 = rep(0.7,3),
                 phi_indoors = rep(0.87,3),
                 phi_bed = rep(0.81,3))

fit_site_parameters$model_seasonality
fit_site_parameters$g0
fit_site_parameters$g
fit_site_parameters$h



year <- 365
month <- 30
sim_length <- 10 * year ## 2015 - 2024
human_population <- 10000
starting_EIR <- 34.22956 ## to be estimated later

simparams <- get_parameters(
  list(
    human_population = human_population,
    # 6 months to younger than 10 years trial cohort
    prevalence_rendering_min_ages = c(0.5, 0.5,    5,    0) * 365, ## Prev in 6 months to 5 years measured, 
    prevalence_rendering_max_ages = c(5,    10,   15,   100) * 365, 
    
    clinical_incidence_rendering_min_ages = c(0.5,0.5, 5,  0) * 365, ## as prev
    clinical_incidence_rendering_max_ages = c(5,  10,  15, 100) * 365,
    
    model_seasonality = TRUE, ## Seasonality 
    
    g0 = fit_site_parameters$g0,
    g = c(fit_site_parameters$g),
    h = c(fit_site_parameters$h),
    
    individual_mosquitoes = FALSE, ## Update next
    carrying_capacity = FALSE
  )
)
# https://pmc.ncbi.nlm.nih.gov/articles/PMC7204121/
# gambiae ss 78.9%, arabiensis 20.9%

#https://link.springer.com/article/10.1186/s13071-021-04748-9
# high proportions of secondary vectors

# Paper from trial not out yet 
# set species
arab_params['Q0'] <- 0.7 # human blood index: assumption
arab_params['mum'] <- 0.112 # daily mortality rate
arab_params['phi_indoors'] <- 0.86 # indoor feeding: from PNAS
arab_params['phi_bednets'] <- 0.8 # bed feeding: from PNAS

gamb_params['Q0'] <- 0.8 # human blood index: assumption
gamb_params['mum'] <- 0.132 # daily mortality rate
gamb_params['phi_indoors'] <- 0.9 # indoor feeding: from PNAS
gamb_params['phi_bednets'] <- 0.85 # bed feeding: from PNAS

simparams <- set_species(simparams,
                         species=list(arab_params, gamb_params, fun_params),
                         proportions=c(0.20,0.60,0.20) ## crude from PMI funestus dominant
)

drug_params <- set_drugs(simparams, list(AL_params, DHA_PQP_params))

treatment_params <- set_clinical_treatment(
  parameters = drug_params,
  drug = 1,
  timesteps =  c(1), # Treatment coverage changes on day 300 and day 600
  coverages =  c(0.4)) # The initial treatment coverage (0%) is the default 
# and does not need to be set

# Set treatment program for DHA-PQP (drug index = 2)
treatment_params <- set_clinical_treatment(
  parameters = treatment_params,
  drug =  2,
  timesteps = c(1),
  coverages =  c(0.3)) 

# Use set_equilibrium to update the parameter set for a given initial EIR
treatment_params <- set_equilibrium(treatment_params, starting_EIR)

simparams <- set_equilibrium(treatment_params, starting_EIR)

## Run the baseline simulation
output0 <- run_simulation(sim_length, simparams)
output0$prpf = output0$p_detect_lm_182.5_3650 / output0$n_age_182.5_3650
output0$inc_6m_10y = 1000 * output0$p_inc_clinical_182.5_3650 / output0$n_age_182.5_3650

## a what if scenario to remove the PBO ITN distribution
# Set ITNs
bednetstimesteps <- c(195, 195+1095)   # The bed nets will be distributed at the end of the first and the 4th year. 

bednetparamscomp <- set_bednets(
  simparams,
  timesteps = bednetstimesteps,
  coverages = c(.6, .8),  # Paper states 98% coverage in trial zone - this is our compromise.
  retention = 2.1 * year, # Nets are kept on average 2.1 years from Andrew Glover 
  dn0 = matrix(c(0.3207, 0.285), nrow = 2, ncol = 3), # Matrix of death probabilities for each mosquito species over time
  rn = matrix(c(0.6530, 0.6797), nrow = 2, ncol = 3), # Matrix of repelling probabilities for each mosquito species over time
  rnm = matrix(c(.24, .24), nrow = 2, ncol = 3), # Matrix of minimum repelling probabilities for each mosquito species over time
  gamman = c(3.492 * 365, 2.909 * 365) # Vector of bed net half-lives for each distribution timestep
)

## Run the ITN only simulation
output1comp <- run_simulation(sim_length, bednetparamscomp)
output1comp$prpf = output1comp$p_detect_lm_182.5_3650 / output1comp$n_age_182.5_3650
output1comp$inc_6m_10y = 1000 * output1comp$p_inc_clinical_182.5_3650 / output1comp$n_age_182.5_3650


# resistant to multiple pyrethroid insecticides
# INTERVENTION DISTRIBUTIONS SEPTEMBER 2021 (Ochomo et al 2025)
#  Busia County received pyrethroid-only ITNs (Permanet 2.0) in 2018. 
#  A mass campaign with pyrethroid plus PBO ITNs (Olyset Plus) was done in July–August, 2021 achieving coverage of up to 98% in the study area
# 2015 pyrethroid only: assumed resistance 0.2	0.3207	0.653026400621824	3.49217977638829
# 2018 pyrethroid only: assumed resistance 0.5	0.285	0.679742902810301	2.90912222304644
# 2021 PBO: assumed resistance 0.8	0.36675	0.615695222697275	3.37981553949983
  
# Set ITNs
bednetstimesteps <- c(195, 195+1095, 195+1095+1095)   # The bed nets will be distributed at the end of the first and the 4th year. 

bednetparams <- set_bednets(
  simparams,
  timesteps = bednetstimesteps,
  coverages = c(0.6, 0.8, 0.9),  # Paper states 98% coverage in trial zone - this is our compromise.
  retention = 2.1 * year, # Nets are kept on average 2.1 years from Andrew Glover 
  dn0 = matrix(c(0.3207, 0.285, 0.36675), nrow = 3, ncol = 3), # Matrix of death probabilities for each mosquito species over time
  rn = matrix(c(0.6530, 0.6797, 0.6157), nrow = 3, ncol = 3), # Matrix of repelling probabilities for each mosquito species over time
  rnm = matrix(c(.24, .24, .24), nrow = 3, ncol = 3), # Matrix of minimum repelling probabilities for each mosquito species over time
  gamman = c(3.492 * 365, 2.909 * 365, 3.3798 * 365) # Vector of bed net half-lives for each distribution timestep
)
  
## Run the ITN only simulation
output1 <- run_simulation(sim_length, bednetparams)
output1$prpf = output1$p_detect_lm_182.5_3650 / output1$n_age_182.5_3650

dtsens = data.frame(reps = 1:50,
                    R0_aa = rep(NA,50),
                    K0_aa = rep(NA,50),
                    M0_aa = rep(NA,50),
                    delta0_aa = rep(NA,50),
                    nR_aa = rep(NA,50),
                    nK_aa = rep(NA,50),
                    nM_aa = rep(NA,50),
                    nd_aa = rep(NA,50),
                    epsR_aa = rep(NA,50),
                    epsK_aa = rep(NA,50),
                    epsM_aa = rep(NA,50),
                    epsd_aa = rep(NA,50),
                    
                    T_ref = 400
)
for(i in 1:50){
  rr = sample(1:4000,size = 1,replace = FALSE)
  
  dtsens$R0_aa[i] = draws_hill_eps_aa$R0[rr]
  dtsens$K0_aa[i] = draws_hill_eps_aa$K0[rr]
  dtsens$M0_aa[i] = draws_hill_eps_aa$M0[rr]
  dtsens$delta0_aa[i] = draws_hill_eps_aa$delta0[rr]
  dtsens$nR_aa[i] = draws_hill_eps_aa$nR[rr]
  dtsens$nK_aa[i] = draws_hill_eps_aa$nK[rr]
  dtsens$nM_aa[i] = draws_hill_eps_aa$nM[rr]
  dtsens$nd_aa[i] = draws_hill_eps_aa$nd[rr]
  dtsens$epsR_aa[i] = draws_hill_eps_aa$eps_R[rr]
  dtsens$epsK_aa[i] = draws_hill_eps_aa$eps_K[rr]
  dtsens$epsM_aa[i] = draws_hill_eps_aa$eps_M[rr]
  dtsens$epsd_aa[i] = draws_hill_eps_aa$eps_d[rr]
  }


# Set Spatial emanators indoors
## start and end time time c(3, 6) * year + 1 
dt_row = 1
spatial_emanator_timesteps = round(seq(310+1095+1095, 310+1095+1095+730, by = 365),0) # Each round of Mosquito Shield spatial emanator is implemented in the 1st and second year.
spatial_emanator_timesteps = 310+1095+1095
spatial_emanator_only_params <- set_spatial_emanator( 
  bednetparams,
  timesteps = spatial_emanator_timesteps,
  coverages = rep(0.8, length(spatial_emanator_timesteps)), # # The first round covers 30% of the population and the second covers 80%. 
  R0 = matrix(rep(c(dtsens$R0_aa[dt_row],dtsens$R0_aa[dt_row],dtsens$R0_aa[dt_row]),each = length(spatial_emanator_timesteps)), ncol=3), # Matrix of impact parameters indoors
  K0 = matrix(rep(c(dtsens$K0_aa[dt_row],dtsens$K0_aa[dt_row],dtsens$K0_aa[dt_row]),each = length(spatial_emanator_timesteps)), ncol=3), # Matrix of impact parameters indoors
  M0 = matrix(rep(c(dtsens$M0_aa[dt_row],dtsens$M0_aa[dt_row],dtsens$M0_aa[dt_row]),each = length(spatial_emanator_timesteps)), ncol=3), # Matrix of impact parameters indoors
  delta0 = matrix(rep(c(dtsens$delta0_aa[dt_row],dtsens$delta0_aa[dt_row],dtsens$delta0_aa[dt_row]),each = length(spatial_emanator_timesteps)), ncol=3), # Matrix of impact parameters indoors
  
  nR = matrix(rep(c(dtsens$nR_aa[dt_row],dtsens$nR_aa[dt_row],dtsens$nR_aa[dt_row]),each = length(spatial_emanator_timesteps)), ncol=3), # Matrix of impact parameters indoors
  nK = matrix(rep(c(dtsens$nK_aa[dt_row],dtsens$nK_aa[dt_row],dtsens$nK_aa[dt_row]),each = length(spatial_emanator_timesteps)), ncol=3), # Matrix of impact parameters indoors
  nM = matrix(rep(c(dtsens$nM_aa[dt_row],dtsens$nM_aa[dt_row],dtsens$nM_aa[dt_row]),each = length(spatial_emanator_timesteps)), ncol=3), # Matrix of impact parameters indoors
  nd = matrix(rep(c(dtsens$nd_aa[dt_row],dtsens$nd_aa[dt_row],dtsens$nd_aa[dt_row]),each = length(spatial_emanator_timesteps)), ncol=3), # Matrix of impact parameters indoors
  
  eps_R = matrix(rep(c(dtsens$epsR_aa[dt_row],dtsens$epsR_aa[dt_row],dtsens$epsR_aa[dt_row]),each = length(spatial_emanator_timesteps)), ncol=3), # Matrix of impact parameters indoors
  eps_K = matrix(rep(c(dtsens$epsK_aa[dt_row],dtsens$epsK_aa[dt_row],dtsens$epsK_aa[dt_row]),each = length(spatial_emanator_timesteps)), ncol=3), # Matrix of impact parameters indoors
  eps_M = matrix(rep(c(dtsens$epsM_aa[dt_row],dtsens$epsM_aa[dt_row],dtsens$epsM_aa[dt_row]),each = length(spatial_emanator_timesteps)), ncol=3), # Matrix of impact parameters indoors
  eps_d = matrix(rep(c(dtsens$epsd_aa[dt_row],dtsens$epsd_aa[dt_row],dtsens$epsd_aa[dt_row]),each = length(spatial_emanator_timesteps)), ncol=3), # Matrix of impact parameters indoors
  
  T_ref = matrix(rep(c(dtsens$T_ref[dt_row],dtsens$T_ref[dt_row],dtsens$T_ref[dt_row]),each = length(spatial_emanator_timesteps)), ncol=3) # Matrix of impact parameters indoors
  
)


## Run the SE simulation
output2 <- run_simulation(sim_length, spatial_emanator_only_params)
output2$prpf1 = output2$p_detect_lm_182.5_3650 / output2$n_age_182.5_3650


## prev
plot(output1comp$prpf ~ output1comp$timestep, type = "l", 
     col = "grey80",lwd = 2,
     xlim = c(2000,3400), ylim = c(0,0.6),
     ylab = "Prevalence in 0.5 - 10 years (%)",
     xlab = "Date",xaxt = "n", yaxt = "n")
axis(1, at = seq(1,10*365,by = 365),labels = seq(2015,2024,by = 1))
axis(2, las = 2, at = seq(0,1,by = 0.2),labels = seq(0,100,by = 20))
lines(output1$prpf ~ output1$timestep, col = "blue", lwd = 2)
lines(output2$prpf ~ output2$timestep, col = "darkred", lwd = 2)

abline(v = spatial_emanator_timesteps, lty = 2, col = "grey80")
abline(v = bednetstimesteps, lty = 2, lwd = 2, col = "grey30")

## Repeat for incidence
# output1comp$inc_6m_10y = 1000*output1comp$p_inc_clinical_182.5_3650/output1comp$n_age_182.5_3650
# output1$inc_6m_10y = 1000*output1$p_inc_clinical_182.5_3650/output1$n_age_182.5_3650
# output2$inc_6m_10y = 1000*output2$p_inc_clinical_182.5_3650/output2$n_age_182.5_3650
# plot(output1comp$inc_6m_10y ~ output1comp$timestep, type = "l", 
#      col = "grey80",lwd = 2,
#      xlim = c(2000,3400),#ylim = c(0,1),
#      ylab = "Clinical incidence per 1,000 children (0.5 - 10 years)",
#      xlab = "Date",xaxt = "n", yaxt = "n")
# axis(1, at = seq(1,10*365,by = 365),labels = seq(2015,2024,by = 1))
# axis(2, las = 2, at = seq(0,10,by = 2))
# lines(output1$inc_6m_10y ~ output1$timestep, col = "blue", lwd = 2)
# lines(output2$inc_6m_10y ~ output2$timestep, col = "darkred", lwd = 2)


## TRACKING UNCERTAINTY THROUGH
store1 <- data.frame(time = output0$timestep)
store_inc1 <- data.frame(time = output0$timestep)
store_ara1 <- data.frame(time = output0$timestep)
store_gam1 <- data.frame(time = output0$timestep)
store_fun1 <- data.frame(time = output0$timestep)
store_eir1gam <- data.frame(time = output0$timestep)
store_eir1ara <- data.frame(time = output0$timestep)
store_eir1fun <- data.frame(time = output0$timestep)

store2 <- data.frame(time = output0$timestep)
store_inc2 <- data.frame(time = output0$timestep)
store_ara2 <- data.frame(time = output0$timestep)
store_gam2 <- data.frame(time = output0$timestep)
store_fun2 <- data.frame(time = output0$timestep)
store_eir2gam <- data.frame(time = output0$timestep)
store_eir2ara <- data.frame(time = output0$timestep)
store_eir2fun <- data.frame(time = output0$timestep)

pyr_only_res_20 = read.csv("C:/Users/esherrar/Documents/RProjects/spatial-emanators/data/pyr_only_res_20.csv")
pyr_only_res_50 = read.csv("C:/Users/esherrar/Documents/RProjects/spatial-emanators/data/pyr_only_res_50.csv")
pyr_pbo_res_80 = read.csv("C:/Users/esherrar/Documents/RProjects/spatial-emanators/data/pyr_PBO_res_80.csv")

for(i in 1:50){
  
  
  # Set ITNs
  bednetstimesteps <- c(195, 195+1095, 195+1095+1095)   # The bed nets will be distributed at the end of the first and the 4th year. 
  
  bednetparams <- set_bednets(
    simparams,
    timesteps = bednetstimesteps,
    coverages = c(.6, .8, .9),  # Paper states 98% coverage in trial zone - this is our compromise.
    retention = 2.1 * year, # Nets are kept on average 2.1 years from Andrew Glover 
    dn0 = matrix(c(pyr_only_res_20[i,2], pyr_only_res_50[i,2], pyr_pbo_res_80[i,2]), nrow = 3, ncol = 3), # Matrix of death probabilities for each mosquito species over time
    rn = matrix(c(pyr_only_res_20[i,3], pyr_only_res_50[i,3], pyr_pbo_res_80[i,3]), nrow = 3, ncol = 3), # Matrix of repelling probabilities for each mosquito species over time
    rnm = matrix(c(.24, .24, .24), nrow = 3, ncol = 3), # Matrix of minimum repelling probabilities for each mosquito species over time
    gamman = c(pyr_only_res_20[i,4] * 365, pyr_only_res_50[i,4] * 365, pyr_pbo_res_80[i,4] * 365) # Vector of bed net half-lives for each distribution timestep
  )
  
  ## Run the ITN only simulation
  output1 <- run_simulation(sim_length, bednetparams)
  output1$prpf = output1$p_detect_lm_182.5_3650 / output1$n_age_182.5_3650
  output1$inc = 1000*output1$p_inc_clinical_182.5_3650 / output1$n_age_182.5_3650
  
  store1[[paste0("draw_", i)]] <- output1$prpf
  store_inc1[[paste0("draw_", i)]] <- output1$inc
  store_ara1[[paste0("draw_",i)]] <- output1$total_M_arab
  store_gam1[[paste0("draw_",i)]] <- output1$total_M_gamb
  store_fun1[[paste0("draw_",i)]] <- output1$total_M_fun
  store_eir1gam[[paste0("draw_",i)]] <- output1$EIR_gamb
  store_eir1ara[[paste0("draw_",i)]] <- output1$EIR_arab
  store_eir1fun[[paste0("draw_",i)]] <- output1$EIR_fun
  
  # i_draw <- sample(x = 1:4000, size = 1, replace = FALSE)  # or sample a random draw if you like
  
  # spatial_emanator_timesteps = round(seq(310+1095+1095, 310+1095+1095+730, by = 30),0) # Each round of Mosquito Shield spatial emanator is implemented in the 1st and second year.
  spatial_emanator_timesteps = c(310+1095+1095,310+1095+1095+365)
  spatial_emanator_only_params <- set_spatial_emanator( 
    bednetparams,
    timesteps = spatial_emanator_timesteps,
    coverages = rep(0.8, length(spatial_emanator_timesteps)), # # The first round covers 30% of the population and the second covers 80%. 
    R0 = matrix(rep(c(dtsens$R0_aa[i],dtsens$R0_aa[i],dtsens$R0_aa[i]),each = length(spatial_emanator_timesteps)), ncol=3), # Matrix of impact parameters indoors
    K0 = matrix(rep(c(dtsens$K0_aa[i],dtsens$K0_aa[i],dtsens$K0_aa[i]),each = length(spatial_emanator_timesteps)), ncol=3), # Matrix of impact parameters indoors
    M0 = matrix(rep(c(dtsens$M0_aa[i],dtsens$M0_aa[i],dtsens$M0_aa[i]),each = length(spatial_emanator_timesteps)), ncol=3), # Matrix of impact parameters indoors
    delta0 = matrix(rep(c(dtsens$delta0_aa[i],dtsens$delta0_aa[i],dtsens$delta0_aa[i]),each = length(spatial_emanator_timesteps)), ncol=3), # Matrix of impact parameters indoors
    
    nR = matrix(rep(c(dtsens$nR_aa[i],dtsens$nR_aa[i],dtsens$nR_aa[i]),each = length(spatial_emanator_timesteps)), ncol=3), # Matrix of impact parameters indoors
    nK = matrix(rep(c(dtsens$nK_aa[i],dtsens$nK_aa[i],dtsens$nK_aa[i]),each = length(spatial_emanator_timesteps)), ncol=3), # Matrix of impact parameters indoors
    nM = matrix(rep(c(dtsens$nM_aa[i],dtsens$nM_aa[i],dtsens$nM_aa[i]),each = length(spatial_emanator_timesteps)), ncol=3), # Matrix of impact parameters indoors
    nd = matrix(rep(c(dtsens$nd_aa[i],dtsens$nd_aa[i],dtsens$nd_aa[i]),each = length(spatial_emanator_timesteps)), ncol=3), # Matrix of impact parameters indoors
    
    eps_R = matrix(rep(c(dtsens$epsR_aa[i],dtsens$epsR_aa[i],dtsens$epsR_aa[i]),each = length(spatial_emanator_timesteps)), ncol=3), # Matrix of impact parameters indoors
    eps_K = matrix(rep(c(dtsens$epsK_aa[i],dtsens$epsK_aa[i],dtsens$epsK_aa[i]),each = length(spatial_emanator_timesteps)), ncol=3), # Matrix of impact parameters indoors
    eps_M = matrix(rep(c(dtsens$epsM_aa[i],dtsens$epsM_aa[i],dtsens$epsM_aa[i]),each = length(spatial_emanator_timesteps)), ncol=3), # Matrix of impact parameters indoors
    eps_d = matrix(rep(c(dtsens$epsd_aa[i],dtsens$epsd_aa[i],dtsens$epsd_aa[i]),each = length(spatial_emanator_timesteps)), ncol=3), # Matrix of impact parameters indoors
    
    T_ref = matrix(rep(c(dtsens$T_ref[i],dtsens$T_ref[i],dtsens$T_ref[i]),each = length(spatial_emanator_timesteps)), ncol=3) # Matrix of impact parameters indoors
    
  )
  
  
  output2 <- run_simulation(sim_length, spatial_emanator_only_params)
  output2$prpf = output2$p_detect_lm_182.5_3650 / output2$n_age_182.5_3650
  output2$inc = 1000*output2$p_inc_clinical_182.5_3650 / output2$n_age_182.5_3650
  
  print(i)
  
  store2[[paste0("draw_", i)]] <- output2$prpf
  store_inc2[[paste0("draw_", i)]] <- output2$inc
  store_ara2[[paste0("draw_",i)]] <- output2$total_M_arab
  store_gam2[[paste0("draw_",i)]] <- output2$total_M_gamb
  store_fun2[[paste0("draw_",i)]] <- output2$total_M_fun
  store_eir2gam[[paste0("draw_",i)]] <- output2$EIR_gamb
  store_eir2ara[[paste0("draw_",i)]] <- output2$EIR_arab
  store_eir2fun[[paste0("draw_",i)]] <- output2$EIR_fun
  
}

# saveRDS(store1,"C:/Users/esherrar/Documents/RProjects/spatial-emanators/simulations/Guardian_store1")
# saveRDS(store_inc1,"C:/Users/esherrar/Documents/RProjects/spatial-emanators/simulations/Guardian_store_inc1")
# saveRDS(store_ara1,"C:/Users/esherrar/Documents/RProjects/spatial-emanators/simulations/Guardian_store_ara1")
# saveRDS(store_gam1,"C:/Users/esherrar/Documents/RProjects/spatial-emanators/simulations/Guardian_store_gam1")
# saveRDS(store_fun1,"C:/Users/esherrar/Documents/RProjects/spatial-emanators/simulations/Guardian_store_fun1")
# saveRDS(store_eir1gam,"C:/Users/esherrar/Documents/RProjects/spatial-emanators/simulations/Guardian_store_eir1gam")
# saveRDS(store_eir1ara,"C:/Users/esherrar/Documents/RProjects/spatial-emanators/simulations/Guardian_store_eir1ara")
# saveRDS(store_eir1fun,"C:/Users/esherrar/Documents/RProjects/spatial-emanators/simulations/Guardian_store_eir1fun")
# saveRDS(store2,"C:/Users/esherrar/Documents/RProjects/spatial-emanators/simulations/Guardian_store2")
# saveRDS(store_inc2,"C:/Users/esherrar/Documents/RProjects/spatial-emanators/simulations/Guardian_store_inc2")
# saveRDS(store_ara2,"C:/Users/esherrar/Documents/RProjects/spatial-emanators/simulations/Guardian_store_ara2")
# saveRDS(store_gam2,"C:/Users/esherrar/Documents/RProjects/spatial-emanators/simulations/Guardian_store_gam2")
# saveRDS(store_fun2,"C:/Users/esherrar/Documents/RProjects/spatial-emanators/simulations/Guardian_store_fun2")
# saveRDS(store_eir2gam,"C:/Users/esherrar/Documents/RProjects/spatial-emanators/simulations/Guardian_store_eir2gam")
# saveRDS(store_eir2ara,"C:/Users/esherrar/Documents/RProjects/spatial-emanators/simulations/Guardian_store_eir2ara")
# saveRDS(store_eir2fun,"C:/Users/esherrar/Documents/RProjects/spatial-emanators/simulations/Guardian_store_eir2fun")


plot(store_ara1[,2] ~ store_ara1[,1], type = "l", 
     col = "white",lwd = 2,
     xlim = c(2000,3400), ylim = c(0,400000),
     ylab = "Mosquito densities",
     xlab = "Date",xaxt = "n", yaxt = "n")
axis(1, at = seq(1,10*365,by = 365),labels = seq(2015,2024,by = 1))
# axis(2, las = 2, at = seq(0,400000,by = 100000),labels = seq(0,400,by = 100))
# lines(output1$prpf ~ output1$timestep, col = "blue", lwd = 2)

for(i in 1:50){
  
  lines(store_ara1[,i+1] ~ store_ara1[,1], col = adegenet::transp("darkblue",0.2), lwd = 2)
  lines(store_ara2[,i+1] ~ store_ara1[,1], col = adegenet::transp("darkred",0.2), lwd = 2)
  lines(store_gam1[,i+1] ~ store_ara1[,1], col = adegenet::transp("blue",0.2), lwd = 1)
  lines(store_gam2[,i+1] ~ store_ara1[,1], col = adegenet::transp("red",0.2), lwd = 1)
  lines(store_fun1[,i+1] ~ store_ara1[,1], col = adegenet::transp("lightblue",0.2), lwd = 1)
  lines(store_fun2[,i+1] ~ store_ara1[,1], col = adegenet::transp("pink",0.2), lwd = 1)
  
}

spatial_emanator_timesteps = 310+1095+1095

abline(v = spatial_emanator_timesteps[1], lty = 2, col = "grey80")
abline(v = spatial_emanator_timesteps[1]+365, lty = 2, col = "grey80")
# abline(v = spatial_emanator_timesteps[1]+365+365, lty = 2, col = "grey80")
abline(v = bednetstimesteps, lty = 2, col = "grey30")




store1 = readRDS("C:/Users/esherrar/Documents/RProjects/spatial-emanators/simulations/Guardian_store1")
store_inc1 = readRDS("C:/Users/esherrar/Documents/RProjects/spatial-emanators/simulations/Guardian_store_inc1")
store2 = readRDS("C:/Users/esherrar/Documents/RProjects/spatial-emanators/simulations/Guardian_store2")
store_inc2 = readRDS("C:/Users/esherrar/Documents/RProjects/spatial-emanators/simulations/Guardian_store_inc2")

store_ara1 = readRDS("C:/Users/esherrar/Documents/RProjects/spatial-emanators/simulations/Guardian_store_ara1")
store_gam1 = readRDS("C:/Users/esherrar/Documents/RProjects/spatial-emanators/simulations/Guardian_store_gam1")
store_fun1 = readRDS("C:/Users/esherrar/Documents/RProjects/spatial-emanators/simulations/Guardian_store_fun1")
store_ara2 = readRDS("C:/Users/esherrar/Documents/RProjects/spatial-emanators/simulations/Guardian_store_ara2")
store_gam2 = readRDS("C:/Users/esherrar/Documents/RProjects/spatial-emanators/simulations/Guardian_store_gam2")
store_fun2 = readRDS("C:/Users/esherrar/Documents/RProjects/spatial-emanators/simulations/Guardian_store_fun2")
spatial_emanator_timesteps = 310+1095+1095

## Prepare figures 
## Final analysis for MS paper... 
(mean(as.numeric(store1[2752,2:51])) - mean(as.numeric(c(store2[2752,2:51]))) ) / 
  mean(as.numeric(store1[2752,2:51]))
(mean(as.numeric(store1[2752,2:51])) - quantile(as.numeric(c(store2[2752,2:51])),0.025) ) / 
  mean(as.numeric(store1[2752,2:51]))
(mean(as.numeric(store1[2752,2:51])) - quantile(as.numeric(c(store2[2752,2:51])),0.975) ) / 
  mean(as.numeric(store1[2752,2:51]))

(mean(as.numeric(store1[3117,2:51])) - mean(as.numeric(c(store2[3117,2:51]))) ) / 
  mean(as.numeric(store1[3117,2:51]))
(mean(as.numeric(store1[3117,2:51])) - quantile(as.numeric(c(store2[3117,2:51])),0.025) ) / 
  mean(as.numeric(store1[3117,2:51]))
(mean(as.numeric(store1[3117,2:51])) - quantile(as.numeric(c(store2[3117,2:51])),0.975) ) / 
  mean(as.numeric(store1[3117,2:51]))

# mean reduction in incidence
fc1 = ft1 = array(dim = c(365,3))
for(t in 2501:2865){
  fc1[t-2500,] = quantile(as.numeric(store_inc1[t,2:51]),c(0.025,0.5,0.975), na.rm = TRUE) 
  ft1[t-2500,] = quantile(as.numeric(store_inc2[t,2:51]),c(0.025,0.5,0.975), na.rm = TRUE)
}
(mean(fc1[,2]) - mean(ft1[,2]) )/mean(fc1[,2])
(mean(fc1[,2]) - mean(ft1[,1]) )/mean(fc1[,2])
(mean(fc1[,2]) - mean(ft1[,3]) )/mean(fc1[,2])

fc1 = ft1 = array(dim = c(365,3))
for(t in 2866:3230){
  fc1[t-2865,] = quantile(as.numeric(store_inc1[t,2:51]),c(0.025,0.5,0.975), na.rm = TRUE) 
  ft1[t-2865,] = quantile(as.numeric(store_inc2[t,2:51]),c(0.025,0.5,0.975), na.rm = TRUE)
}
(mean(fc1[,2]) - mean(ft1[,2]) )/mean(fc1[,2])
(mean(fc1[,2]) - mean(ft1[,1]) )/mean(fc1[,2])
(mean(fc1[,2]) - mean(ft1[,3]) )/mean(fc1[,2])

fc1 = ft1 = array(dim = c(730,3))
for(t in 2501:3231){
  fc1[t-2501,] = quantile(as.numeric(store_inc1[t,2:51]),c(0.025,0.5,0.975), na.rm = TRUE) 
  ft1[t-2501,] = quantile(as.numeric(store_inc2[t,2:51]),c(0.025,0.5,0.975), na.rm = TRUE)
}
(mean(fc1[,2]) - mean(ft1[,2]) )/mean(fc1[,2])
(mean(fc1[,2]) - mean(ft1[,1]) )/mean(fc1[,2])
(mean(fc1[,2]) - mean(ft1[,3]) )/mean(fc1[,2])



## for A plot the relative reduction in mosquitoes
summary_fun = function(store_objectA,
                       store_objectB){
  fc1 = ft1 = array(dim = c(365,3))
  
  for(t in 2501:2865){
    fc1[t-2500,] = quantile(as.numeric(store_objectA[t,2:51]),c(0.025,0.5,0.975), na.rm = TRUE) 
    ft1[t-2500,] = quantile(as.numeric(store_objectB[t,2:51]),c(0.025,0.5,0.975), na.rm = TRUE)
  }
  rel_low = (fc1[1:365,2] - ft1[1:365,1])/fc1[1:365,2]
  rel_mean = (fc1[1:365,2] - ft1[1:365,2])/fc1[1:365,2]
  rel_upp = (fc1[1:365,2] - ft1[1:365,3])/fc1[1:365,2]
  
  return(c(
    mean(rel_mean, na.rm = TRUE),
    mean(rel_upp, na.rm = TRUE),
    mean(rel_low, na.rm = TRUE)
  ))
}
summary_fun = function(store_objectA,
                       store_objectB){
  fc1 = ft1 = array(dim = c(365,3))
  
  for(t in 2501:2865){
    fc1[t-2500,] = quantile(as.numeric(store_objectA[t,2:51]),c(0.025,0.5,0.975), na.rm = TRUE) 
    ft1[t-2500,] = quantile(as.numeric(store_objectB[t,2:51]),c(0.025,0.5,0.975), na.rm = TRUE)
  }
  rel_low = (fc1[1:365,2] - ft1[1:365,1])/fc1[1:365,2]
  rel_mean = (fc1[1:365,2] - ft1[1:365,2])/fc1[1:365,2]
  rel_upp = (fc1[1:365,2] - ft1[1:365,3])/fc1[1:365,2]
  
  return(c(
    mean(rel_mean, na.rm = TRUE),
    mean(rel_upp, na.rm = TRUE),
    mean(rel_low, na.rm = TRUE)
  ))
}
summary_fun(store1,store2)
summary_fun(store_inc1,store_inc2)
summary_fun(store_ara1,store_ara2)
summary_fun(store_gam1,store_gam2)
summary_fun(store_fun1,store_fun2)
summary_fun(store_eir1gam,store_eir2gam)
summary_fun(store_eir1ara,store_eir2ara)
summary_fun(store_eir1fun,store_eir2fun)

summary_fun_y2 = function(store_objectA,
                          store_objectB){
  fc1 = ft1 = array(dim = c(365,3))
  
  for(t in 2866:3230){
    fc1[t-2865,] = quantile(as.numeric(store_objectA[t,2:51]),c(0.025,0.5,0.975), na.rm = TRUE) 
    ft1[t-2865,] = quantile(as.numeric(store_objectB[t,2:51]),c(0.025,0.5,0.975), na.rm = TRUE)
  }
  rel_low = (fc1[1:365,2] - ft1[1:365,1])/fc1[1:365,2]
  rel_mean = (fc1[1:365,2] - ft1[1:365,2])/fc1[1:365,2]
  rel_upp = (fc1[1:365,2] - ft1[1:365,3])/fc1[1:365,2]
  
  return(c(
    mean(rel_mean, na.rm = TRUE),
    mean(rel_upp, na.rm = TRUE),
    mean(rel_low, na.rm = TRUE)
  ))
}
summary_fun_y2(store1,store2)
summary_fun_y2(store_inc1,store_inc2)
summary_fun_y2(store_ara1,store_ara2)
summary_fun_y2(store_gam1,store_gam2)
summary_fun_y2(store_fun1,store_fun2)
summary_fun_y2(store_eir1gam,store_eir2gam)
summary_fun_y2(store_eir1ara,store_eir2ara)
summary_fun_y2(store_eir1fun,store_eir2fun)

summary_fun_both = function(store_objectA,
                            store_objectB){
  fc1 = ft1 = array(dim = c(730,3))
  
  for(t in 2501:3230){
    fc1[t-2500,] = quantile(as.numeric(store_objectA[t,2:51]),c(0.025,0.5,0.975), na.rm = TRUE) 
    ft1[t-2500,] = quantile(as.numeric(store_objectB[t,2:51]),c(0.025,0.5,0.975), na.rm = TRUE)
  }
  rel_low = (fc1[1:730,2] - ft1[1:730,1])/fc1[1:730,2]
  rel_mean = (fc1[1:730,2] - ft1[1:730,2])/fc1[1:730,2]
  rel_upp = (fc1[1:730,2] - ft1[1:730,3])/fc1[1:730,2]
  
  return(c(
    mean(rel_mean, na.rm = TRUE),
    mean(rel_upp, na.rm = TRUE),
    mean(rel_low, na.rm = TRUE)
  ))
}
store_objectA = store_inc1
store_objectB = store_inc2
fc1 = ft1 = array(dim = c(730,3))
for(t in 2501:3230){
  fc1[t-2500,] = quantile(as.numeric(store_inc1[t,2:51]),c(0.025,0.5,0.975), na.rm = TRUE) 
  ft1[t-2500,] = quantile(as.numeric(store_inc2[t,2:51]),c(0.025,0.5,0.975), na.rm = TRUE)
}

mean(fc1[1:730,2] - ft1[1:730,1])/mean(fc1[1:730,2])
mean(fc1[1:730,2] - ft1[1:730,2])/mean(fc1[1:730,2])
mean(fc1[1:730,2] - ft1[1:730,3])/mean(fc1[1:730,2])

summary_fun_both(store1,store2)
summary_fun_both(store_inc1,store_inc2)
summary_fun_both(store_gam1,store_gam2)
summary_fun_both(store_ara1,store_ara2)

species = c("gambiae", "arabiensis", "funestus")
reductions = c(summary_fun(store_gam1,store_gam2)[1], summary_fun(store_ara1,store_ara2)[1], summary_fun(store_fun1,store_fun2)[1], 
               summary_fun_y2(store_gam1,store_gam2)[1], summary_fun_y2(store_ara1,store_ara2)[1], summary_fun_y2(store_fun1,store_fun2)[1])
lower =      c(summary_fun(store_gam1,store_gam2)[2], summary_fun(store_ara1,store_ara2)[2], summary_fun(store_fun1,store_fun2)[2], 
               summary_fun_y2(store_gam1,store_gam2)[2], summary_fun_y2(store_ara1,store_ara2)[2], summary_fun_y2(store_fun1,store_fun2)[2])
upper =      c(summary_fun(store_gam1,store_gam2)[3], summary_fun(store_ara1,store_ara2)[3], summary_fun(store_fun1,store_fun2)[3], 
               summary_fun_y2(store_gam1,store_gam2)[3], summary_fun_y2(store_ara1,store_ara2)[3], summary_fun_y2(store_fun1,store_fun2)[3])#y2

results = c("did_HLC","CDC 1st year","CDC 2nd year")
MEAN_EMPIRICL_results = c(0.274, 0.3,  0.44)
low_EMPIRICL_results = c(-0.108, 0.11, 0.05)
upp_EMPIRICL_results =  c(0.519, 0.45, 0.64)


baseC = as.numeric(store1[2400,2:51])
baseT = as.numeric(store2[2400,2:51])
FOL_C1 = as.numeric(store1[2765,2:51])
FOL_T1 = as.numeric(store2[2765,2:51])
FOL_C2 = as.numeric(store1[3130,2:51])
FOL_T2 = as.numeric(store2[3130,2:51])


## prevalence reductions
quantile(unlist(store1[2501:2865,]),c(0.025,0.5,0.975))
quantile(unlist(store2[2501:2865,]),c(0.025,0.5,0.975))

(quantile(unlist(store1[2501:2865,]),0.5)-
    quantile(unlist(store2[2501:2865,]),c(0.025,0.5,0.975)))/
  quantile(unlist(store1[2501:2865,]),0.5)

(quantile(unlist(store1[2866:3230,]),0.5)-
    quantile(unlist(store2[2866:3230,]),c(0.025,0.5,0.975)))/
  quantile(unlist(store1[2866:3230,]),0.5)

(quantile(unlist(store1[2501:3230,]),0.5)-
    quantile(unlist(store2[2501:3230,]),c(0.025,0.5,0.975)))/
  quantile(unlist(store1[2501:3230,]),0.5)



#############################
##
## PLOT Figure 4


library(ggplot2)
library(dplyr)

df <- data.frame( 
  period = c( "Overall", 
              "Year 1", "Year 1", "Year 1", "Year 1", "Year 1", 
              "Year 2", "Year 2", "Year 2", "Year 2" ), 
  label = c( "Core region overall", "HLC", "Light traps", 
             "gambiae", "arabiensis", "funestus", 
             "Light traps", 
             "gambiae", "arabiensis", "funestus" ), 
  estimate = c(0.115, 
               MEAN_EMPIRICL_results[1], MEAN_EMPIRICL_results[2], 
               reductions[1], reductions[2], reductions[3], 
               MEAN_EMPIRICL_results[3], 
               reductions[4], reductions[5], reductions[6] ), 
  lower = c( -0.416, 
             low_EMPIRICL_results[1], low_EMPIRICL_results[2], 
             lower[1], lower[2], lower[3], 
             low_EMPIRICL_results[3], 
             lower[4], lower[5], lower[6] ), 
  upper = c( 1.129, 
             upp_EMPIRICL_results[1], upp_EMPIRICL_results[2], 
             upper[1], upper[2], upper[3], 
             upp_EMPIRICL_results[3], 
             upper[4], upper[5], upper[6] ), 
  source = c( "Overall", "Trial", "Trial", 
              "Simulated", "Simulated", "Simulated", 
              "Trial", 
              "Simulated", "Simulated", "Simulated" ), 
  species = c( "Combined", 
               "Combined", "Combined",
               "gambiae", "arabiensis", "funestus", 
               "Combined", 
               "gambiae", "arabiensis", "funestus" )
)

df$period <- factor( df$period, 
                     levels = c("Overall", "Year 1", "Year 2") ) 
df$label <- factor( df$label, 
                    levels = c( "Core region overall", 
                                "HLC", "Light traps", 
                                "gambiae", "arabiensis", "funestus" ) ) 
df$source <- factor( df$source, 
                     levels = c("Overall", "Trial", "Simulated") ) 
df$species <- factor( df$species, levels = c("Combined",
                                             "gambiae", "arabiensis", "funestus") )


# panA = ggplot(df, 
#               aes(x = label, y = estimate)) + 
#   # Confidence intervals 
#   geom_errorbar( 
#     aes( ymin = lower, ymax = upper, colour = source ), 
#     width = 0, linewidth = 0.8, linetype = "dashed" ) + 
#   # Point estimates 
#   geom_point( aes( colour = source, shape = species ), size = 3.5 ) + 
#   # Zero/reference line 
#   geom_hline( yintercept = 0, colour = "grey75", linewidth = 0.5 ) + 
#   # Facets 
#   facet_wrap( ~ period, nrow = 1, scales = "free_x" ) + 
#   # Colours 
#   scale_colour_manual( values = c( "Overall" = "grey30", 
#                                    "Trial" = "#2878B5", 
#                                    "Simulated" = "#A23E48" ), 
#                        drop = FALSE ) +
#   # Shapes 
#   scale_shape_manual( values = c( "Combined" = 19, 
#                                   "gambiae" = 16, "arabiensis" = 17, "funestus" = 18), 
#                       drop = FALSE ) + 
#   # Y axis 
#   scale_y_continuous( limits = c(-0.5, 1.2), 
#                       breaks = seq(-0.4, 1.2, by = 0.2), 
#                       labels = function(x) paste0(x * 100, "%"), 
#                       expand = expansion(mult = c(0.02, 0.05)) ) + 
#   # Labels 
#   labs( x = NULL, y = "Reduction in mosquito density (%)", 
#         colour = NULL, shape = NULL ) + 
#   # Clean theme 
#   guides( shape = "none" ) +
#   theme_classic(base_size = 13) + theme( 
#     # Facet labels 
#     strip.background = element_blank(), 
#     strip.text = element_text( size = 14, face = "bold", colour = "grey20" ), 
#     # X labels 
#     axis.text.x = element_text( angle = 45, 
#                                 hjust = 1, vjust = 1, 
#                                 colour = "grey20" ), 
#     # Y labels 
#     axis.text.y = element_text( colour = "grey20" ), 
#     axis.title.y = element_text( face = "bold", 
#                                  margin = margin(r = 12) ), 
#     # Legend 
#     legend.position = "none", 
#     # legend.box = "vertical", 
#     # Spacing 
#     panel.spacing = unit(2, "lines"), 
#     plot.margin = margin( 10, 15, 10, 10 ) )



############
##
## Panel 2 can be the prevalence trends with point estimates overlaid as boxplots
library(ggplot2)
library(dplyr)
library(tidyr)


# =========================================================
# 1. Create data frame for the simulated trajectories
# =========================================================

# Time
time <- store_ara1[, 1]

# Store 1: blue trajectories
traj1 <- as.data.frame(store1[, -1])
names(traj1) <- paste0("sim", 1:ncol(traj1))

traj1$time <- time

traj1_long <- traj1 %>%
  pivot_longer(
    cols = starts_with("sim"),
    names_to = "simulation",
    values_to = "prevalence"
  ) %>%
  mutate(
    model = "Model 1"
  )


# Store 2: dark-red trajectories
traj2 <- as.data.frame(store2[, -1])
names(traj2) <- paste0("sim", 1:ncol(traj2))

traj2$time <- time

traj2_long <- traj2 %>%
  pivot_longer(
    cols = starts_with("sim"),
    names_to = "simulation",
    values_to = "prevalence"
  ) %>%
  mutate(
    model = "Model 2"
  )


# Combine the two sets of trajectories
trajectories <- bind_rows(
  traj1_long,
  traj2_long
)


# =========================================================
# 2. Create data frame for prevalence estimates
# =========================================================

# ---------------------------------------------------------
# Replace these with your actual prevalence estimates.
#
# You want one vector for each time point:
#
#   baseline
#   year 1
#   year 2
#
# Each vector should contain the individual prevalence
# estimates that you want represented by the boxplot.
# ---------------------------------------------------------
baseC = as.numeric(store1[2400,2:51])
baseT = as.numeric(store2[2400,2:51])
FOL_C1 = as.numeric(store1[2765,2:51])
FOL_T1 = as.numeric(store2[2765,2:51])
FOL_C2 = as.numeric(store1[3130,2:51])
FOL_T2 = as.numeric(store2[3130,2:51])

prevalence_estimates <- data.frame(
  
  period = c(
    rep("Baseline treatment", length(baseT)),
    rep("Year 1 (SE)",    length(FOL_T1)),
    rep("Year 2 (SE)",    length(FOL_T2)),
    rep("Baseline placebo", length(baseC)),
    rep("Year 1 (control)",    length(FOL_C1)),
    rep("Year 2 (control)",    length(FOL_C2))
    
  ),
  
  prevalence = c(baseT, FOL_T1, FOL_T2, baseC, FOL_C1, FOL_C2 )
)


# # Set order of the three periods
# prevalence_estimates$period <- factor(
#   prevalence_estimates$period,
#   levels = c(
#     "Baseline treatment","Baseline placebo",
#     "Year 1 (SE)", "Year 1 (Placebo)",
#     "Year 2 (SE)", "Year 2 (Placebo)"
#   )
# )


# =========================================================
# 3. Identify the x positions for the boxplots
# =========================================================

# These should correspond to the actual intervention /
# observation dates in your simulation.

boxplot_times <- data.frame(
  
  period = factor(
    c("Baseline treatment", "Baseline placebo",
      "Year 1 (SE)", "Year 1 (control)",
      "Year 2 (SE)", "Year 2 (control)"),
    levels = c("Baseline treatment", "Baseline placebo",
               "Year 1 (SE)", "Year 1 (control)",
               "Year 2 (SE)", "Year 2 (control)")
  ),
  
  time = c(
    spatial_emanator_timesteps[1]-90,spatial_emanator_timesteps[1]-90,
    spatial_emanator_timesteps[1]-90 + 365,spatial_emanator_timesteps[1]-90 + 365,
    spatial_emanator_timesteps[1]-90 + 2 * 365,spatial_emanator_timesteps[1]-90 + 2 * 365
  )
)


# Add the actual time to each prevalence observation
prevalence_estimates <- prevalence_estimates %>%
  left_join(
    boxplot_times,
    by = "period"
  )

prevalence_estimates <- prevalence_estimates %>%
  mutate(
    group = ifelse(
      grepl("treatment|SE", period),
      "Treatment",
      "Control"
    )
  )

# ==============================================
panB <- ggplot() + 
  geom_line( data = trajectories %>% 
               filter(model == "Model 1"), 
             aes( x = time, y = prevalence, group = simulation ), 
             colour = scales::alpha("blue", 0.20), linewidth = 0.7 ) + 
  
  geom_line( data = trajectories %>% 
               filter(model == "Model 2"), 
             aes( x = time, y = prevalence, group = simulation ), 
             colour = scales::alpha("darkred", 0.20), linewidth = 0.7 ) + 
  
  geom_vline( xintercept = c(spatial_emanator_timesteps,spatial_emanator_timesteps+365), linetype = "dashed", 
              colour = "grey80", linewidth = 0.7 ) + 
  
  geom_vline( xintercept = bednetstimesteps, 
              linetype = "dashed", colour = "grey30", linewidth = 0.7 ) + 
  
  geom_boxplot(
    data = prevalence_estimates,
    aes(
      x = time,
      y = prevalence,
      group = period,
      fill = group
    ),
    width = 35,
    position = position_dodge(width = 45),
    colour = "black",
    alpha = 0.8,
    linewidth = 0.7,
    outlier.shape = 21,
    outlier.fill = "white",
    outlier.colour = "black"
  ) +
  
  scale_fill_manual(
    values = c(
      "Treatment" = "darkred",
      "Control" = "blue"
    )
  ) +
  
  guides(
    fill = guide_legend(title = NULL)
  ) +
  
  # geom_jitter(
  #   data = prevalence_estimates,
  #   aes(
  #     x = time,
  #     y = prevalence,
  #     group = period
  #   ),
  #   position = position_jitterdodge(
  #     jitter.width = 8,
  #     dodge.width = 45
  #   ),
  #   size = 1.5,
  #   colour = "black",
  #   alpha = 0.4
  # ) +
  scale_x_continuous( 
    limits = c(2000, 3400), 
    breaks = seq( 1, 10 * 365, by = 365 ), 
    labels = seq( 2015, 2024, by = 1 ), 
    expand = expansion(mult = c(0, 0.01)) ) + 
  scale_y_continuous( limits = c(0, 0.6), 
                      breaks = seq( 0, 0.6, by = 0.2 ), 
                      labels = function(x) { paste0(x * 100, "%") }, 
                      expand = expansion(mult = c(0, 0.03)) ) + 
  
  labs( x = NULL, 
        y = "Model simulated prevalence \nin 6 months to 10 years (%)" ) + 
  
  theme_classic(base_size = 11) +
  theme(
    axis.text = element_text(
      colour = "grey20",
      size = 9
    ),
    axis.title.y = element_text(
      face = "bold",
      size = 10,
      margin = margin(r = 8)
    ),
    axis.text.x = element_text(
      colour = "grey20",
      size = 9
    ),
    legend.position = "none",
    plot.margin = margin(10, 15, 10, 10)
  )
  
  # theme_classic( base_size = 13 ) + 
  # theme( axis.text = element_text( colour = "grey20" ), 
  #        axis.title.y = element_text( face = "bold", margin = margin(r = 12) ), 
  #        axis.text.x = element_text( colour = "grey20" ), 
  #        legend.position = "none", 
  #        plot.margin = margin( 10, 15, 10, 10 ) )

#######################################
##
## And for clinical incidence
## Panel 3

# =========================================================
# 1. Create data frame for the simulated trajectories
# =========================================================

# Time
time <- store_ara1[, 1]

# Store 1: blue trajectories
traj1 <- as.data.frame(store1[, -1])
names(traj1) <- paste0("sim", 1:ncol(traj1))

traj1$time <- time

traj1_long <- traj1 %>%
  pivot_longer(
    cols = starts_with("sim"),
    names_to = "simulation",
    values_to = "prevalence"
  ) %>%
  mutate(
    model = "Model 1"
  )


# Store 2: dark-red trajectories
traj2 <- as.data.frame(store2[, -1])
names(traj2) <- paste0("sim", 1:ncol(traj2))

traj2$time <- time

traj2_long <- traj2 %>%
  pivot_longer(
    cols = starts_with("sim"),
    names_to = "simulation",
    values_to = "prevalence"
  ) %>%
  mutate(
    model = "Model 2"
  )


# Combine the two sets of trajectories
trajectories <- bind_rows(
  traj1_long,
  traj2_long
)


# =========================================================
# 2. Create data frame for prevalence estimates
# =========================================================

# ---------------------------------------------------------
# Replace these with your actual prevalence estimates.
#
# You want one vector for each time point:
#
#   baseline
#   year 1
#   year 2
#
# Each vector should contain the individual prevalence
# estimates that you want represented by the boxplot.
# ---------------------------------------------------------

# =========================================================
# Incidence boxplot data
# =========================================================

# Testing times: every 30 days for 2 years
testing_times <- round(seq(310+1095+1095, 310+1095+1095+730, by = 30),0) # Each round of Mosquito Shield spatial emanator is implemented in the 1st and second year.

# Treatment and control incidence at every testing time
incT <- as.numeric(store_inc2[testing_times[1], 2:51])
incC <- as.numeric(store_inc1[testing_times[1], 2:51])


# ---------------------------------------------------------
# Build long data frame
# ---------------------------------------------------------

incidence_estimates <- bind_rows(
  
  # Treatment
  lapply(testing_times, function(t) {
    
    data.frame(
      time = t,
      period = "Treatment",
      simulation = 1:50,
      incidence = as.numeric(store_inc2[t, 2:51])
    )
    
  }) %>% bind_rows(),
  
  # Control
  lapply(testing_times, function(t) {
    
    data.frame(
      time = t,
      period = "Control",
      simulation = 1:50,
      incidence = as.numeric(store_inc1[t, 2:51])
    )
    
  }) %>% bind_rows()
)


# Check
dim(incidence_estimates)

table(incidence_estimates$period)

incidence_estimates$subsets = paste(incidence_estimates$time, incidence_estimates$period)

# =========================================================
# Incidence trajectories
# =========================================================

time_inc <- store_inc1[, 1]

traj_inc_C <- as.data.frame(store_inc1[, -1])
names(traj_inc_C) <- paste0("sim", 1:ncol(traj_inc_C))

traj_inc_C$time <- time_inc

traj_inc_C_long <- traj_inc_C %>%
  pivot_longer(
    cols = starts_with("sim"),
    names_to = "simulation",
    values_to = "incidence"
  ) %>%
  mutate(
    model = "Control"
  )


traj_inc_T <- as.data.frame(store_inc2[, -1])
names(traj_inc_T) <- paste0("sim", 1:ncol(traj_inc_T))

traj_inc_T$time <- time_inc

traj_inc_T_long <- traj_inc_T %>%
  pivot_longer(
    cols = starts_with("sim"),
    names_to = "simulation",
    values_to = "incidence"
  ) %>%
  mutate(
    model = "Treatment"
  )


incidence_trajectories <- bind_rows(
  traj_inc_C_long,
  traj_inc_T_long
)

p_inc <- ggplot() +
  
  # Control trajectories
  geom_line(
    data = incidence_trajectories %>%
      filter(model == "Control"),
    aes(
      x = time,
      y = incidence,
      group = simulation
    ),
    colour = scales::alpha("blue", 0.02),
    linewidth = 0.7
  ) +
  
  # Treatment trajectories
  geom_line(
    data = incidence_trajectories %>%
      filter(model == "Treatment"),
    aes(
      x = time,
      y = incidence,
      group = simulation
    ),
    colour = scales::alpha("darkred", 0.02),
    linewidth = 0.7
  ) +
  
  # Your intervention/reference lines
  geom_vline(
    xintercept = c(spatial_emanator_timesteps,spatial_emanator_timesteps+365),
    linetype = "dashed",
    colour = "grey80",
    linewidth = 0.7
  ) +
  
  geom_vline(
    xintercept = bednetstimesteps,
    linetype = "dashed",
    colour = "grey30",
    linewidth = 0.7
  ) +
  
  # Boxplots
  geom_boxplot(
    data = incidence_estimates,
    aes(
      x = time,
      y = incidence,
      group = subsets,
      fill = period
    ),
    width = 20,
    position = position_dodge(width = 25),
    colour = "black",
    alpha = 0.8,
    linewidth = 0.6,
    outlier.shape = 21,
    outlier.fill = "white",
    outlier.colour = "black"
  ) +
  
  # # Individual simulation points
  # geom_jitter(
  #   data = incidence_estimates,
  #   aes(
  #     x = time,
  #     y = incidence,
  #     group = period
  #   ),
  #   position = position_jitterdodge(
  #     jitter.width = 5,
  #     dodge.width = 25
  #   ),
  #   size = 1.2,
  #   colour = "black",
  #   alpha = 0.25
  # ) +
  
  scale_fill_manual(
    values = c(
      "Treatment" = "darkred",
      "Control" = "blue"
    )
  ) +
  
  scale_x_continuous(
    limits = c(2000, 3400),
    breaks = seq(1, 10 * 365, by = 365),
    labels = seq(2015, 2024, by = 1),
    expand = expansion(mult = c(0, 0.01))
  ) +
  
  labs(
    x = NULL,
    y = "Model simulated incidence",
    fill = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.text = element_text(
      colour = "grey20",
      size = 9
    ),
    axis.title.y = element_text(
      face = "bold",
      size = 10,
      margin = margin(r = 8)
    ),
    axis.text.x = element_text(
      colour = "grey20",
      size = 9
    ),
    legend.position = "none",
    plot.margin = margin(10, 15, 10, 10)
  )
  # theme_classic(base_size = 13) +
  # 
  # theme(
  #   axis.text = element_text(colour = "grey20"),
  #   axis.title.y = element_text(
  #     face = "bold",
  #     margin = margin(r = 12)
  #   ),
  #   legend.position = "none",
  #   plot.margin = margin(10, 15, 10, 10)
  # )

library(patchwork)

(panB / p_inc) +
  plot_layout( heights = c(1, 1) ) +
  plot_annotation( tag_levels = "A" ) 
