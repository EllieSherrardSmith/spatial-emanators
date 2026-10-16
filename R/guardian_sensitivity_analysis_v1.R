# devtools::load_all()
# devtools::document()

# library("malariaEquilibrium")
# library("malariasimulation")
# library(cali)
# # Sensitivity Analysis for Guardian

## Vary waning use (if Richard can add!)

## phi_indoors / in bed

## ITN resistance parameters

## picking Kenya setting as example:
year <- 365
month <- 30
sim_length <- 10 * year ## 2015 - 2024
human_population <- 30000
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

## we will alternate phi_indoors and phi_bednets
gamb_params['Q0'] <- 0.8 # human blood index: assumption
gamb_params['mum'] <- 0.132 # daily mortality rate
## LOOP THROUGH PARAMTERS FOR UNCERTAINTY
gamb_params['phi_indoors'] <- 0.9 # indoor feeding: from PNAS
gamb_params['phi_bednets'] <- 0.85 # bed feeding: from PNAS

simparams <- set_species(simparams,
                         species=list(gamb_params),
                         proportions=c(1) ## crude from PMI funestus dominant
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
# Set ITNs (campaigns every 3 years)
bednetstimesteps <- c(195, 195+1095, 195+1095+1095)   # The bed nets will be distributed at the end of the first and the 4th year. 

## LOOP THROUGH RESISTANCE PARAMETERS FOR ITNS
bednetparamscomp <- set_bednets(
  simparams,
  timesteps = bednetstimesteps,
  coverages = c(.8, .8, .8),  # Paper states 98% coverage in trial zone - this is our compromise.
  retention = 2.1 * year, # Nets are kept on average 2.1 years from Andrew Glover 
  dn0 = matrix(c(0.285), nrow = 3, ncol = 1), # Matrix of death probabilities for each mosquito species over time
  rn = matrix(c(0.6797), nrow = 3, ncol = 1), # Matrix of repelling probabilities for each mosquito species over time
  rnm = matrix(c(.24, .24, .24), nrow = 3, ncol = 1), # Matrix of minimum repelling probabilities for each mosquito species over time
  gamman = c( 2.909 * 365,  2.909 * 365,  2.909 * 365) # Vector of bed net half-lives for each distribution timestep
)

output1 <- run_simulation(sim_length, bednetparamscomp)
output1$prpf = output1$p_detect_lm_182.5_3650 / output1$n_age_182.5_3650
output1$inc_6m_10y = 1000 * output1$p_inc_clinical_182.5_3650 / output1$n_age_182.5_3650


## SE PARAMETERS FOR ANNUAL USE OF GUARDIAN STARTING AT 80%
rr_tests = c(1000,500,300,150,50)

rr = rr_tests[1]
SEparamscomp <- set_spatial_emanator(bednetparamscomp,
                     timesteps = c(195+1095+1095,195+1095+1095+365,195+1095+1095+730),
                     coverages = c(0.9,0.9,0.9),
                     
                     R0 = matrix(0.286, nrow = 3, ncol = 1), 
                     nR = matrix(7.938, nrow = 3, ncol = 1), 
                     eps_R = matrix(0.029, nrow = 3, ncol = 1),
                     
                     K0 = matrix(0.075, nrow = 3, ncol = 1), 
                     nK = matrix(19.98, nrow = 3, ncol = 1), 
                     eps_K = matrix(0.103, nrow = 3, ncol = 1),
                     
                     M0 = matrix(0.065, nrow = 3, ncol = 1), 
                     nM = matrix(6.30, nrow = 3, ncol = 1), 
                     eps_M = matrix(0.059, nrow = 3, ncol = 1),
                     
                     delta0 = matrix(0.424, nrow = 3, ncol = 1), 
                     nd = matrix(10.5, nrow = 3, ncol = 1), 
                     eps_d = matrix(0.101, nrow = 3, ncol = 1),
                     
                     T_ref = matrix(400, nrow = 3, ncol = 1), 
                     retention = c(rr,rr,rr))

output2 <- run_simulation(sim_length, SEparamscomp)
output2$prpf = output2$p_detect_lm_182.5_3650 / output2$n_age_182.5_3650
output2$inc_6m_10y = 1000 * output2$p_inc_clinical_182.5_3650 / output2$n_age_182.5_3650

plot(output0$prpf ~ output0$timestep, col = "white",
     ylim = c(0,0.6), xlim = c(1095+1095,1095+1095+1395))
lines(output1$prpf ~ output1$timestep, col = "blue",lwd = 2)
# lines(output1$prpf ~ output1$timestep, col = "blue")
lines(output2$prpf ~ output2$timestep, col = "orange",lwd = 2)
abline(v = c(195+1095+1095,195+1095+1095+365,195+1095+1095+730),
       lty = 2)
abline(v = c(195+1095+1095,195+1095+1095+1095),
       lty = 1)


## LOOP THROUGH SE PARAMETERS FOR ANNUAL USE OF GUARDIAN STARTING AT 80%
rr_tests = seq(50, 500, length = 6)
phi_indoors_values = seq(0.4,1,length = 5)

sims_loop_comparator_f = function(phi_indoors_values){
  gamb_params['phi_indoors'] <- phi_indoors_values # indoor feeding: from PNAS
  gamb_params['phi_bednets'] <- phi_indoors_values - 0.05 # bed feeding: from PNAS
  
  simparams <- set_species(simparams,
                           species=list(gamb_params),
                           proportions=c(1) ## crude from PMI funestus dominant
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
  
  ### Bednet comparison only needs to run once per phi estimate
  # Set ITNs (campaigns every 3 years)
  bednetstimesteps <- c(195, 195+1095, 195+1095+1095)   # The bed nets will be distributed at the end of the first and the 4th year. 
  
  ## LOOP THROUGH RESISTANCE PARAMETERS FOR ITNS
  bednetparamscomp <- set_bednets(
    simparams,
    timesteps = bednetstimesteps,
    coverages = c(.8, .8, .8),  # Paper states 98% coverage in trial zone - this is our compromise.
    retention = 2.1 * year, # Nets are kept on average 2.1 years from Andrew Glover 
    dn0 = matrix(c(0.285), nrow = 3, ncol = 1), # Matrix of death probabilities for each mosquito species over time
    rn = matrix(c(0.6797), nrow = 3, ncol = 1), # Matrix of repelling probabilities for each mosquito species over time
    rnm = matrix(c(.24, .24, .24), nrow = 3, ncol = 1), # Matrix of minimum repelling probabilities for each mosquito species over time
    gamman = c( 2.909 * 365,  2.909 * 365,  2.909 * 365) # Vector of bed net half-lives for each distribution timestep
  )
  
  output1 <- run_simulation(sim_length, bednetparamscomp)
  output1$prpf = output1$p_detect_lm_182.5_3650 / output1$n_age_182.5_3650
  output1$inc_6m_10y = 1000 * output1$p_inc_clinical_182.5_3650 / output1$n_age_182.5_3650
  ## output a value for cases over a 3-year ITN campaign and retain this
  cases = sum(output1$inc_6m_10y[2385:3480])
  
  return(list(output1$prpf,
              cases))
}

aa = list()
for(i in 1:length(phi_indoors_values)){
  aa[[i]] = sims_loop_comparator_f(phi_indoors_values[i])
}

comparator = data.frame(timesteps = output1$timestep)
comparator$phi_indoors_1 = c(aa[[1]][[1]][1:3650])
comparator$phi_indoors_2 = c(aa[[2]][[1]][1:3650])
comparator$phi_indoors_3 = c(aa[[3]][[1]][1:3650])
comparator$phi_indoors_4 = c(aa[[4]][[1]][1:3650])
comparator$phi_indoors_5 = c(aa[[5]][[1]][1:3650])

compare_cases = numeric(5)
compare_cases[1] = aa[[1]][2]
compare_cases[2] = aa[[2]][2]
compare_cases[3] = aa[[3]][2]
compare_cases[4] = aa[[4]][2]
compare_cases[5] = aa[[5]][2]


rr_tests = seq(50, 500, length = 6)
phi_indoors_values = seq(0.4,1,length = 5)

sims_loop_SE_f = function(rr_tests, phi_indoors_values){
  gamb_params['phi_indoors'] <- phi_indoors_values # indoor feeding: from PNAS
  gamb_params['phi_bednets'] <- phi_indoors_values - 0.05 # bed feeding: from PNAS
  
  simparams <- set_species(simparams,
                           species=list(gamb_params),
                           proportions=c(1) ## crude from PMI funestus dominant
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
  
  ### Bednet comparison only needs to run once per phi estimate
  # Set ITNs (campaigns every 3 years)
  bednetstimesteps <- c(195, 195+1095, 195+1095+1095)   # The bed nets will be distributed at the end of the first and the 4th year. 
  
  ## LOOP THROUGH RESISTANCE PARAMETERS FOR ITNS
  bednetparamscomp <- set_bednets(
    simparams,
    timesteps = bednetstimesteps,
    coverages = c(.8, .8, .8),  # Paper states 98% coverage in trial zone - this is our compromise.
    retention = 2.1 * year, # Nets are kept on average 2.1 years from Andrew Glover 
    dn0 = matrix(c(0.285), nrow = 3, ncol = 1), # Matrix of death probabilities for each mosquito species over time
    rn = matrix(c(0.6797), nrow = 3, ncol = 1), # Matrix of repelling probabilities for each mosquito species over time
    rnm = matrix(c(.24, .24, .24), nrow = 3, ncol = 1), # Matrix of minimum repelling probabilities for each mosquito species over time
    gamman = c( 2.909 * 365,  2.909 * 365,  2.909 * 365) # Vector of bed net half-lives for each distribution timestep
  )
  
  SEparamscomp <- set_spatial_emanator(bednetparamscomp,
                                       timesteps = c(195+1095+1095,195+1095+1095+365,195+1095+1095+730),
                                       coverages = c(0.9,0.9,0.9),
                                       
                                       R0 = matrix(0.286, nrow = 3, ncol = 1), 
                                       nR = matrix(7.938, nrow = 3, ncol = 1), 
                                       eps_R = matrix(0.029, nrow = 3, ncol = 1),
                                       
                                       K0 = matrix(0.075, nrow = 3, ncol = 1), 
                                       nK = matrix(19.98, nrow = 3, ncol = 1), 
                                       eps_K = matrix(0.103, nrow = 3, ncol = 1),
                                       
                                       M0 = matrix(0.065, nrow = 3, ncol = 1), 
                                       nM = matrix(6.30, nrow = 3, ncol = 1), 
                                       eps_M = matrix(0.059, nrow = 3, ncol = 1),
                                       
                                       delta0 = matrix(0.424, nrow = 3, ncol = 1), 
                                       nd = matrix(10.5, nrow = 3, ncol = 1), 
                                       eps_d = matrix(0.101, nrow = 3, ncol = 1),
                                       
                                       T_ref = matrix(400, nrow = 3, ncol = 1), 
                                       retention = c(rr_tests,rr_tests,rr_tests))
  
  output2 <- run_simulation(sim_length, SEparamscomp)
  
  ## output a vector for prevalence trend and retain this
  output2$prpf = output2$p_detect_lm_182.5_3650 / output2$n_age_182.5_3650
  output2$inc_6m_10y = 1000 * output2$p_inc_clinical_182.5_3650 / output2$n_age_182.5_3650
  
  ## output a value for cases over a 3-year ITN campaign and retain this
  cases = sum(output2$inc_6m_10y[2385:3480])
  
  return(list(output2$prpf,
              cases))
}

ab1 = ab2 = ab3 = ab4 = ab5 = ab6 = list()
for(i in 1:length(phi_indoors_values)){
  ab1[[i]] = sims_loop_SE_f(rr_tests = rr_tests[1],
                            phi_indoors_values = phi_indoors_values[i])
  ab2[[i]] = sims_loop_SE_f(rr_tests = rr_tests[2],
                            phi_indoors_values = phi_indoors_values[i])
  ab3[[i]] = sims_loop_SE_f(rr_tests = rr_tests[3],
                            phi_indoors_values = phi_indoors_values[i])
  ab4[[i]] = sims_loop_SE_f(rr_tests = rr_tests[4],
                            phi_indoors_values = phi_indoors_values[i])
  ab5[[i]] = sims_loop_SE_f(rr_tests = rr_tests[5],
                            phi_indoors_values = phi_indoors_values[i])
  ab6[[i]] = sims_loop_SE_f(rr_tests = rr_tests[6],
                            phi_indoors_values = phi_indoors_values[i])
}

trt_rr1 = trt_rr2 = trt_rr3 = trt_rr4 = trt_rr5 = trt_rr6 = data.frame(timesteps = output1$timestep)
trt_rr1$phi_indoors_1_rr1 = unlist(ab1[[1]][[1]][1:3650])
trt_rr1$phi_indoors_2_rr1 = unlist(ab1[[2]][[1]][1:3650])
trt_rr1$phi_indoors_3_rr1 = unlist(ab1[[3]][[1]][1:3650])
trt_rr1$phi_indoors_4_rr1 = unlist(ab1[[4]][[1]][1:3650])
trt_rr1$phi_indoors_5_rr1 = unlist(ab1[[5]][[1]][1:3650])

trt_rr2$phi_indoors_1_rr2 = unlist(ab2[[1]][[1]][1:3650])
trt_rr2$phi_indoors_2_rr2 = unlist(ab2[[2]][[1]][1:3650])
trt_rr2$phi_indoors_3_rr2 = unlist(ab2[[3]][[1]][1:3650])
trt_rr2$phi_indoors_4_rr2 = unlist(ab2[[4]][[1]][1:3650])
trt_rr2$phi_indoors_5_rr2 = unlist(ab2[[5]][[1]][1:3650])

trt_rr3$phi_indoors_1_rr3 = unlist(ab3[[1]][[1]][1:3650])
trt_rr3$phi_indoors_2_rr3 = unlist(ab3[[2]][[1]][1:3650])
trt_rr3$phi_indoors_3_rr3 = unlist(ab3[[3]][[1]][1:3650])
trt_rr3$phi_indoors_4_rr3 = unlist(ab3[[4]][[1]][1:3650])
trt_rr3$phi_indoors_5_rr3 = unlist(ab3[[5]][[1]][1:3650])

trt_rr4$phi_indoors_1_rr4 = unlist(ab4[[1]][[1]][1:3650])
trt_rr4$phi_indoors_2_rr4 = unlist(ab4[[2]][[1]][1:3650])
trt_rr4$phi_indoors_3_rr4 = unlist(ab4[[3]][[1]][1:3650])
trt_rr4$phi_indoors_4_rr4 = unlist(ab4[[4]][[1]][1:3650])
trt_rr4$phi_indoors_5_rr4 = unlist(ab4[[5]][[1]][1:3650])

trt_rr5$phi_indoors_1_rr5 = unlist(ab5[[1]][[1]][1:3650])
trt_rr5$phi_indoors_2_rr5 = unlist(ab5[[2]][[1]][1:3650])
trt_rr5$phi_indoors_3_rr5 = unlist(ab5[[3]][[1]][1:3650])
trt_rr5$phi_indoors_4_rr5 = unlist(ab5[[4]][[1]][1:3650])
trt_rr5$phi_indoors_5_rr5 = unlist(ab5[[5]][[1]][1:3650])

trt_rr6$phi_indoors_1_rr6 = unlist(ab6[[1]][[1]][1:3650])
trt_rr6$phi_indoors_2_rr6 = unlist(ab6[[2]][[1]][1:3650])
trt_rr6$phi_indoors_3_rr6 = unlist(ab6[[3]][[1]][1:3650])
trt_rr6$phi_indoors_4_rr6 = unlist(ab6[[4]][[1]][1:3650])
trt_rr6$phi_indoors_5_rr6 = unlist(ab6[[5]][[1]][1:3650])

## make a data.frame with timesteps, retention and prevalence
prev_long = data.frame(timesteps = rep(output1$timestep,6),
                            retention = rep(rr_tests,each = length(output1$timestep)),
                            prevalence = c(trt_rr1$phi_indoors_4_rr1,
                                           trt_rr2$phi_indoors_4_rr2,
                                           trt_rr3$phi_indoors_4_rr3,
                                           trt_rr4$phi_indoors_4_rr4,
                                           trt_rr5$phi_indoors_4_rr5,
                                           trt_rr6$phi_indoors_4_rr6) )

library(ggplot2)

panelC = ggplot(prev_long,
       aes(x = timesteps,
           y = prevalence,
           colour = retention,
           group = retention)) +
  
  # Shaded region for three-year comparison
  annotate(
    "rect",
    xmin = 2385, xmax = 3480,
    ymin = 0, ymax = 0.6,
    fill = "grey80",
    alpha = 0.35
  ) +
  
  # Prevalence trajectories
  geom_line(linewidth = 0.8) +
  
  # Treatment deployment
  geom_vline(
    xintercept = c(195+1095+1095),
    linetype = "solid",
    colour = "darkblue",
    linewidth = 0.7
  ) +
  
  # Treatment deployment
  geom_vline(
    xintercept = c(195+1095+1095,195+1095+1095+365,195+1095+1095+730),
    linetype = "dashed",
    colour = "black",
    linewidth = 0.7
  ) +
  
  
  
  # Label for shaded region
  annotate(
    "text",
    x = 2932, y = 0.57,
    label = "Cases compared over three years",
    colour = "grey30",
    size = 3.5
  ) +
  
  labs(
    x = "Time (days)",
    y = "Prevalence",
    colour = "Retention time"
  ) +
  
  coord_cartesian(
    xlim = c(2000, 3600),
    ylim = c(0, 0.6)
  ) +
  
  theme_classic()


compare_cases_se = data.frame(phi = phi_indoors_values,
                              se_rr1 = 0,se_rr2 = 0,se_rr3 = 0,
                              se_rr4 = 0,se_rr5 = 0,se_rr6 = 0)
compare_cases_se$se_rr1[1] = ab1[[1]][2]
compare_cases_se$se_rr1[2] = ab1[[2]][2]
compare_cases_se$se_rr1[3] = ab1[[3]][2]
compare_cases_se$se_rr1[4] = ab1[[4]][2]
compare_cases_se$se_rr1[5] = ab1[[5]][2]

compare_cases_se$se_rr2[1] = ab2[[1]][2]
compare_cases_se$se_rr2[2] = ab2[[2]][2]
compare_cases_se$se_rr2[3] = ab2[[3]][2]
compare_cases_se$se_rr2[4] = ab2[[4]][2]
compare_cases_se$se_rr2[5] = ab2[[5]][2]

compare_cases_se$se_rr3[1] = ab3[[1]][2]
compare_cases_se$se_rr3[2] = ab3[[2]][2]
compare_cases_se$se_rr3[3] = ab3[[3]][2]
compare_cases_se$se_rr3[4] = ab3[[4]][2]
compare_cases_se$se_rr3[5] = ab3[[5]][2]

compare_cases_se$se_rr4[1] = ab4[[1]][2]
compare_cases_se$se_rr4[2] = ab4[[2]][2]
compare_cases_se$se_rr4[3] = ab4[[3]][2]
compare_cases_se$se_rr4[4] = ab4[[4]][2]
compare_cases_se$se_rr4[5] = ab4[[5]][2]

compare_cases_se$se_rr5[1] = ab5[[1]][2]
compare_cases_se$se_rr5[2] = ab5[[2]][2]
compare_cases_se$se_rr5[3] = ab5[[3]][2]
compare_cases_se$se_rr5[4] = ab5[[4]][2]
compare_cases_se$se_rr5[5] = ab5[[5]][2]

compare_cases_se$se_rr6[1] = ab6[[1]][2]
compare_cases_se$se_rr6[2] = ab6[[2]][2]
compare_cases_se$se_rr6[3] = ab6[[3]][2]
compare_cases_se$se_rr6[4] = ab6[[4]][2]
compare_cases_se$se_rr6[5] = ab6[[5]][2]


##########################################
##
## create a cases averted matrix
control_cases = c(compare_cases[[1]][1],
                  compare_cases[[2]][1],
                  compare_cases[[3]][1],
                  compare_cases[[4]][1],
                  compare_cases[[5]][1])

phi_40 = as.numeric((control_cases[1] - unlist(compare_cases_se[1,2:7])) / control_cases[1])
phi_55 = as.numeric((control_cases[2] - unlist(compare_cases_se[2,2:7])) / control_cases[2])
phi_70 = as.numeric((control_cases[3] - unlist(compare_cases_se[3,2:7])) / control_cases[3])
phi_85 = as.numeric((control_cases[4] - unlist(compare_cases_se[4,2:7])) / control_cases[4])
phi_99 = as.numeric((control_cases[5] - unlist(compare_cases_se[5,2:7])) / control_cases[5])

cases_matrix <- cbind(
  phi_40,
  phi_55,
  phi_70,
  phi_85,
  phi_99
)

cases_matrix
cases_matrix = round(cases_matrix *100,0)
rownames(cases_matrix) <- c(
  "Retention 1", "Retention 2", "Retention 3",
  "Retention 4", "Retention 5", "Retention 6"
)
library(tidyverse)

cases_df <- as.data.frame(cases_matrix) %>%
  mutate(retention = row_number()) %>%
  pivot_longer(
    cols = starts_with("phi_"),
    names_to = "phi",
    values_to = "cases_averted"
  )
retention_values <- seq(50, 500, length.out = 6)

panelD = ggplot(cases_df,
                aes(x = phi, y = factor(retention,
                                        levels = 1:6,
                                        labels = retention_values),
                    fill = cases_averted)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = round(cases_averted, 3))) +
  scale_fill_viridis_c() +
  labs(
    x = expression(phi),
    y = "Retention time (days)"
  ) +
  guides(fill = "none") +
  theme_minimal()

#### !! RUN EVERYTHING FROM PLOTS LINES 729 
#### IN guardian_simulations_v1.R

library(patchwork)

(panB + p_inc) / (panelC + panelD) +
  plot_layout( heights = c(1, 1) ) +
  plot_annotation( tag_levels = "A" ) 
