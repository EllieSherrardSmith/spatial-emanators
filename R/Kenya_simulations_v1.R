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
map_Data <- read.csv("C:/Users/esherrar/Documents/RProjects/spatial-emanators/simulating projections for Unitaid/Subnational Unit-data MAP 31 July 2025.csv", header = TRUE)
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
annual_pfpr_summary <- function(x){
  #set the years to target
  year <- ceiling(x$timestep[3651:8030] / 365)
  #estimate the Pfpr
  pfpr <- x$n_detect_lm_730_3650[3651:8030]  / x$n_age_730_3650[3651:8030]
  ##tabulise
  tapply(pfpr, year, mean)
}

# Establish a target PfPR2-10 value:
target_pfpr <- prevalence_2_10yrs_MAP

# Add a parameter to the parameter list specifying the number of timesteps 
fit_site_parameters$timesteps <- sim_length

# Run the calibrate() function:
set.seed(123)
out <- calibrate(
  parameters = fit_site_parameters,
  target = target_pfpr[2:13],
  summary_function = annual_pfpr_summary,
  eq_prevalence = target_pfpr[1]
)
## out - ESTIMATED AT 34.22956

fit_site_parameters$pev = TRUE
fit_site_parameters$pev_epi_coverages = rep(0,25)
fit_site_parameters$spraying = FALSE

fit_site_parameters$human_population <- 30000
fit_site_parameters <- set_equilibrium(fit_site_parameters, init_EIR = 34.22956)

output0 <- run_simulation(timesteps = sim_length, parameters = fit_site_parameters)  
## Generate prevalence
output0$pfpr <- output0$n_detect_lm_0_36500  / output0$n_age_0_36500

# store[[1]] = output0
# params_store[[1]] = fit_site_parameters
# fit_site_parameters$total_M
# 
# saveRDS(store[[1]],paste0("C:/Users/esherrar/Documents/RProjects/spatial-emanators/CF_Sims_store_Busia_KENtrial"))
# saveRDS(params_store[[1]],paste0("C:/Users/esherrar/Documents/RProjects/spatial-emanators/Base_parameters_store_Busia_KENtrial"))

# cases_pp_all_age_cf = array(dim = c(2,4))
# people_age_cohort_cf = array(dim = c(2,7))
# for(i in 1:2){
#   ## all-age
#   cases_pp_all_age_cf[i,1] = sum(output0$n_inc_clinical_0_36500[((20+i)*365+315):((22+i)*365-1+315)])/
#     mean(output0$n_age_0_36500[((20+i)*365+315):((22+i)*365-1+315)])
#   people_age_cohort_cf[i,1] = mean(output0$n_age_0_36500[((20+i)*365+315):((22+i)*365-1+315)])
#   ## U5
#   cases_pp_all_age_cf[i,2] = sum(output0$n_inc_clinical_0_1825[((20+i)*365+315):((22+i)*365-1+315)])/
#     mean(output0$n_age_0_1825[((20+i)*365+315):((22+i)*365-1+315)])
#   people_age_cohort_cf[i,2] = mean(output0$n_age_0_1825[((20+i)*365+315):((22+i)*365-1+315)])
#   ## 5-15
#   cases_pp_all_age_cf[i,3] = sum(output0$n_inc_clinical_1825_5475[((20+i)*365+315):((22+i)*365-1+315)])/
#     mean(output0$n_age_1825_5475[((20+i)*365+315):((22+i)*365-1+315)])
#   people_age_cohort_cf[i,3] = mean(output0$n_age_1825_5475[((20+i)*365+315):((22+i)*365-1+315)])
#   ## >15
#   cases_pp_all_age_cf[i,4] = sum(output0$n_inc_clinical_5475_36500[((20+i)*365+315):((22+i)*365-1+315)])/
#     mean(output0$n_age_5475_36500[((20+i)*365+315):((22+i)*365-1+315)])
#   people_age_cohort_cf[i,4] = mean(output0$n_age_5475_36500[((20+i)*365+315):((22+i)*365-1+315)])
# }
# ## percentage of people in each age cohort
# people_age_cohort_cf[,5] = people_age_cohort_cf[,2]/people_age_cohort_cf[,1]
# people_age_cohort_cf[,6] = people_age_cohort_cf[,3]/people_age_cohort_cf[,1]
# people_age_cohort_cf[,7] = people_age_cohort_cf[,4]/people_age_cohort_cf[,1]

par(mfrow = c(1,1))
par(mar = c(4,4,1,1))
plot(output0$pfpr ~ output0$timestep,
     type = "l",
     xlim = c(20*365,25*365),xaxt="n",
     yaxt = "n",ylim = c(0,1.15),
     xlab = "Year",ylab = "",
     bty="n")
mtext("PfPR all age LM (%)", side = 2, line = 3)
mtext("Intervention use (%)", side = 2, line = 2)
# points(prevalence_2_10yrs_MAP[2:13] ~ checker_year, pch = 19, col ="red")
axis(1,at = seq(1,25*365,365),labels = 2000:2024)
axis(2, las = 2, at = seq(0,0.8,0.2),
     labels = seq(0,80,20))
# abline(v = 21*365)
# polygon(c(output0$timestep[(14*365):(25*365)],
#           rev(output1$timestep[(14*365):(25*365)])),
#         c(output1$pfpr[(14*365):(25*365)],
#           rev(output0$pfpr[(14*365):(25*365)])),
#         col = adegenet::transp("blue",0.2),border=NA)

# text(365*20,1.1,"Busia County", cex=1.1)
add_itn = KEN_df$interventions$itn_use[KEN_df$interventions$name_1 == "Busia"]
add_smc = KEN_df$interventions$smc_cov[KEN_df$interventions$name_1 == "Busia"]
add_trt = KEN_df$interventions$tx_cov[KEN_df$interventions$name_1 == "Busia"]
add_trt_act = add_trt*KEN_df$interventions$prop_act[KEN_df$interventions$name_1 == "Busia"]
annual_year = 365*unique(ceiling((1:(25*365))/365))

lines(add_itn ~
        annual_year, col = "purple",lwd=2)
lines(add_smc ~
        annual_year, col = "#CC38A1",lwd=2)
lines(add_trt ~
        annual_year, col = "aquamarine2",lwd=2)
lines(add_trt_act ~
        annual_year, col = "aquamarine2",lwd=2, lty = 2)
legend("topleft", legend = c("ITN","SMC","Treatment","(with ACT)"),
       lty = c(1,1,1,2),
       lwd = 2,
       col = c("purple","#CC38A1","aquamarine2","aquamarine2"),
       bty="n")

## Now figure out how to adjust the parameter list to include spatial emanators,
## and remove anything / change anything to reflect trial better... 

##############################
##
## read in the extracted parameters from the statistical model
sts = readRDS("C:/Users/esherrar/Documents/RProjects/spatial-emanators/manuscript/simulations/MSb_stan_model_outputNgufor") 
stsT = readRDS("C:/Users/esherrar/Documents/RProjects/spatial-emanators/manuscript/simulations/MSbtanz_stan_model_outputTanzania")
## pick a random row from these sims
rr = sample(1:4000,size = 1,replace = FALSE)
fed_a1 = sts$beta1[rr] + sts$beta3[rr]
fed_a2 = sts$beta2[rr]
feddead_a1 = sts$alphaf1[rr] + sts$alphaf3[rr]
feddead_a2 = sts$alphaf2[rr]
unfeddead_a1 = sts$alpha1[rr] + sts$alpha3[rr]
unfeddead_a2 =  sts$alpha2[rr]
det_a1 = sts$omega1[rr]
det_a2 = sts$omega2[rr]

Tfed_a1 = stsT$beta1[rr] + stsT$beta3[rr]
Tfed_a2 = stsT$beta2[rr]
Tfeddead_a1 = stsT$alphaf1[rr] + stsT$alphaf3[rr]
Tfeddead_a2 = stsT$alphaf2[rr]
Tunfeddead_a1 = stsT$alpha1[rr] + stsT$alpha3[rr]
Tunfeddead_a2 =  stsT$alpha2[rr]
Tdet_a1 = stsT$omega1[rr]
Tdet_a2 = stsT$omega2[rr]

time = 1:365
test_mort_unfed = 1/(1 + exp(-(unfeddead_a1 + unfeddead_a2*time)))
test_mort_fed = 1/(1 + exp(-(feddead_a1 + feddead_a2*time)))
test_fed =0.69/(1 + exp(-(fed_a1 + fed_a2*time)))
test_det = 1/(1 + exp(-(det_a1 + det_a2*time)))

Ttest_mort_unfed = 1/(1 + exp(-(Tunfeddead_a1 + Tunfeddead_a2*time)))
Ttest_mort_fed = 1/(1 + exp(-(Tfeddead_a1 + Tfeddead_a2*time)))
Ttest_fed =0.69/(1 + exp(-(Tfed_a1 + Tfed_a2*time)))
Ttest_det = 1/(1 + exp(-(Tdet_a1 + Tdet_a2*time)))


plot(test_mort_unfed ~ time, type = "l", col = "darkblue", 
     lwd = 2, xlim = c(0,100),ylim = c(0,1))
lines(test_mort_fed ~ time,col = "blue", lwd = 2)
lines(test_fed ~ time,col = "darkred", lwd = 2)
lines(test_det ~ time,col = "orange", lwd = 2)

lines(Ttest_mort_unfed ~ time,col = "darkblue", lwd = 2, lty = 2)
lines(Ttest_mort_fed ~ time,col = "blue", lwd = 2, lty = 2)
lines(Ttest_fed ~ time,col = "darkred", lwd = 2, lty = 2)
lines(Ttest_det ~ time,col = "orange", lwd = 2, lty = 2)
##
##############################



dtsens = data.frame(eir = rep(34.22956,3),
                 Q0 = rep(0.7,3),
                 phi_indoors = rep(0.87,3),
                 phi_bed = rep(0.81,3))#,
                 # fed_a1 = c(-0.573, -0.84, -0.26),
                 # fed_a2 = c(0.088, 0.07, 0.10),
                 # feddead_a1 = c(-4.3, -5.26, -3.27),
                 # feddead_a2 = c(-0.01, -0.05, 0.016),
                 # unfeddead_a1 = c(0.68, 0.32, 1.03),
                 # unfeddead_a2 = c(-0.08, -0.097, -0.07),
                 # det_a1 = c(-1.31, -1.69, -0.83),
                 # det_a2 = c(-0.05, -0.07, -0.03))

## what this does...
spraying_decay <- function(t, theta, gamma) {
  1 / (1 + exp(-(theta + gamma * t)))
}
time = 1:365
lFse  = spraying_decay(time, -4.213, -0.015)
lUFse = spraying_decay(time, 1.337, -0.082)
mse   = spraying_decay(time, -1.252, -0.052)
kse   = 0.699 * spraying_decay(time, -3.527, 0.093)

jse = 1 - lUFse - lFse - kse

mse_comp <- (1 - mse)
lUFse_prime <- lUFse * mse_comp
lFse_prime <- lFse * mse_comp
kse_prime <- kse * mse_comp
jse_prime <- jse * mse_comp + mse
lse_prime <- (lUFse + lFse) * mse_comp

prob_spraying_repels <- function(ls_prime, ks_prime, js_prime, k0) {
  (1 - ks_prime / k0) * (js_prime / (ls_prime + js_prime))
}

rse_in = prob_spraying_repels(lse_prime, kse_prime, jse_prime, 0.699)
dt_row = 1

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

dtsens = data.frame(reps = 1:20,
                    fed_a1 = rep(NA,20),
                    fed_a2 = rep(NA,20),
                    feddead_a1 = rep(NA,20),
                    feddead_a2 = rep(NA,20),
                    unfeddead_a1 = rep(NA,20),
                    unfeddead_a2 = rep(NA,20),
                    det_a1 = rep(NA,20),
                    det_a2 = rep(NA,20),
                    
                    Tfed_a1 = rep(NA,20),
                    Tfed_a2 = rep(NA,20),
                    Tfeddead_a1 = rep(NA,20),
                    Tfeddead_a2 = rep(NA,20),
                    Tunfeddead_a1 = rep(NA,20),
                    Tunfeddead_a2 = rep(NA,20),
                    Tdet_a1 = rep(NA,20),
                    Tdet_a2 = rep(NA,20))
for(i in 1:20){
  rr = sample(1:4000,size = 1,replace = FALSE)
  fed_a1 = sts$beta1[rr] + sts$beta3[rr]
  fed_a2 = sts$beta2[rr]
  feddead_a1 = sts$alphaf1[rr] + sts$alphaf3[rr]
  feddead_a2 = sts$alphaf2[rr]
  unfeddead_a1 = sts$alpha1[rr] + sts$alpha3[rr]
  unfeddead_a2 =  sts$alpha2[rr]
  det_a1 = sts$omega1[rr]
  det_a2 = sts$alpha2[rr]
  
  Tfed_a1 = stsT$beta1[rr] + stsT$beta3[rr]
  Tfed_a2 = stsT$beta2[rr]
  Tfeddead_a1 = stsT$alphaf1[rr] + stsT$alphaf3[rr]
  Tfeddead_a2 = stsT$alphaf2[rr]
  Tunfeddead_a1 = stsT$alpha1[rr] + stsT$alpha3[rr]
  Tunfeddead_a2 =  stsT$alpha2[rr]
  Tdet_a1 = stsT$omega1[rr]
  Tdet_a2 = stsT$alpha2[rr]
  
  dtsens$fed_a1[i] = fed_a1 
  dtsens$fed_a2[i] = fed_a2
  dtsens$feddead_a1[i] = feddead_a1
  dtsens$feddead_a2[i] = feddead_a2
  dtsens$unfeddead_a1[i] = unfeddead_a1
  dtsens$unfeddead_a2[i] = unfeddead_a2
  dtsens$det_a1[i] = det_a1
  dtsens$det_a2[i] = det_a2
  
  dtsens$Tfed_a1[i] = Tfed_a1 
  dtsens$Tfed_a2[i] = Tfed_a2
  dtsens$Tfeddead_a1[i] = Tfeddead_a1
  dtsens$Tfeddead_a2[i] = Tfeddead_a2
  dtsens$Tunfeddead_a1[i] = Tunfeddead_a1
  dtsens$Tunfeddead_a2[i] = Tunfeddead_a2
  dtsens$Tdet_a1[i] = Tdet_a1
  dtsens$Tdet_a2[i] = Tdet_a2
  
}
time = 1:365
test_mort_unfed = 1/(1 + exp(-(dtsens$unfeddead_a1[1] + dtsens$unfeddead_a2[1]*time)))
test_mort_fed = 1/(1 + exp(-(dtsens$feddead_a1[1] + dtsens$feddead_a2[1]*time)))
test_fed =0.69/(1 + exp(-(dtsens$fed_a1[1] + dtsens$fed_a2[1]*time)))
test_det = 1/(1 + exp(-(dtsens$det_a1[1] + dtsens$det_a2[1]*time)))

Ttest_mort_unfed = 1/(1 + exp(-(dtsens$Tunfeddead_a1[1] + dtsens$Tunfeddead_a2[1]*time)))
Ttest_mort_fed = 1/(1 + exp(-(dtsens$Tfeddead_a1[1] + dtsens$Tfeddead_a2[1]*time)))
Ttest_fed =0.69/(1 + exp(-(dtsens$Tfed_a1[1] + dtsens$Tfed_a2[1]*time)))
Ttest_det = 1/(1 + exp(-(dtsens$Tdet_a1[1] + dtsens$Tdet_a2[1]*time)))

plot(test_mort_unfed ~ time, type = "l", col = "darkblue", 
     lwd = 2, xlim = c(0,100),ylim = c(0,1))
for(i in 1:20){
  test_mort_unfed = 1/(1 + exp(-(dtsens$unfeddead_a1[i] + dtsens$unfeddead_a2[i]*time)))
  test_mort_fed = 1/(1 + exp(-(dtsens$feddead_a1[i] + dtsens$feddead_a2[i]*time)))
  test_fed =0.69/(1 + exp(-(dtsens$fed_a1[i] + dtsens$fed_a2[i]*time)))
  test_det = 1/(1 + exp(-(dtsens$det_a1[i] + dtsens$det_a2[i]*time)))
  
  Ttest_mort_unfed = 1/(1 + exp(-(dtsens$Tunfeddead_a1[i] + dtsens$Tunfeddead_a2[i]*time)))
  Ttest_mort_fed = 1/(1 + exp(-(dtsens$Tfeddead_a1[i] + dtsens$Tfeddead_a2[i]*time)))
  Ttest_fed =0.69/(1 + exp(-(dtsens$Tfed_a1[i] + dtsens$Tfed_a2[i]*time)))
  Ttest_det = 1/(1 + exp(-(dtsens$Tdet_a1[i] + dtsens$Tdet_a2[i]*time)))
  
  lines(test_mort_unfed ~ time,col = "darkblue", lwd = 2)
  lines(test_mort_fed ~ time,col = "blue", lwd = 2)
  lines(test_fed ~ time,col = "darkred", lwd = 2)
  lines(test_det ~ time,col = "orange", lwd = 2)
  
  lines(Ttest_mort_unfed ~ time,col = "darkblue", lwd = 2, lty = 2)
  lines(Ttest_mort_fed ~ time,col = "blue", lwd = 2, lty = 2)
  lines(Ttest_fed ~ time,col = "darkred", lwd = 2, lty = 2)
  lines(Ttest_det ~ time,col = "orange", lwd = 2, lty = 2)
  
}


# Set Spatial emanators indoors
## start and end time time c(3, 6) * year + 1 
dt_row = 1
spatial_emanator_timesteps = round(seq(310+1095+1095, 310+1095+1095+730, by = 30),0) # Each round of Mosquito Shield spatial emanator is implemented in the 1st and second year.

spatial_emanator_only_params <- set_spatial_emanator(
  bednetparams,
  timesteps = spatial_emanator_timesteps,
  coverages = rep(0.8, length(spatial_emanator_timesteps)), # # The first round covers 30% of the population and the second covers 80%. 
  kse_in_a1 = matrix(rep(c(dtsens$Tfed_a1[dt_row],dtsens$fed_a1[dt_row],dtsens$fed_a1[dt_row]),each = length(spatial_emanator_timesteps)), ncol=3), # Matrix of impact parameters indoors
  kse_in_a2 = matrix(rep(c(dtsens$Tfed_a2[dt_row],dtsens$fed_a2[dt_row],dtsens$fed_a2[dt_row]),each = length(spatial_emanator_timesteps)), ncol=3), # Matrix of impact parameters indoors
  dfse_in_a1 = matrix(rep(c(dtsens$Tfeddead_a1[dt_row],dtsens$feddead_a1[dt_row],dtsens$feddead_a1[dt_row]),each = length(spatial_emanator_timesteps)), ncol=3), # Matrix of impact parameters indoors
  dfse_in_a2 = matrix(rep(c(dtsens$Tfeddead_a2[dt_row],dtsens$feddead_a2[dt_row],dtsens$feddead_a2[dt_row]),each = length(spatial_emanator_timesteps)), ncol=3), # Matrix of impact parameters indoors
  dufse_in_a1 = matrix(rep(c(dtsens$Tunfeddead_a1[dt_row],dtsens$unfeddead_a1[dt_row],dtsens$unfeddead_a1[dt_row]),each = length(spatial_emanator_timesteps)), ncol=3), # Matrix of impact parameters indoors
  dufse_in_a2 = matrix(rep(c(dtsens$Tunfeddead_a2[dt_row],dtsens$unfeddead_a2[dt_row],dtsens$unfeddead_a2[dt_row]),each = length(spatial_emanator_timesteps)), ncol=3), # Matrix of impact parameters indoors
  detse_in_a1 = matrix(rep(c(dtsens$Tdet_a1[dt_row],dtsens$det_a1[dt_row],dtsens$det_a1[dt_row]),each = length(spatial_emanator_timesteps)), ncol=3), # Matrix of impact parameters indoors
  detse_in_a2 = matrix(rep(c(dtsens$Tdet_a2[dt_row],dtsens$det_a2[dt_row],dtsens$det_a2[dt_row]),each = length(spatial_emanator_timesteps)), ncol=3) # Matrix of impact parameters indoors
 
)


## Run the SE simulation
output2 <- run_simulation(sim_length, spatial_emanator_only_params)
output2$prpf1 = output2$p_detect_lm_182.5_3650 / output2$n_age_182.5_3650

output2l <- run_simulation(sim_length, spatial_emanator_only_params)
output2l$prpf = output2l$p_detect_lm_182.5_3650 / output2l$n_age_182.5_3650

output2u <- run_simulation(sim_length, spatial_emanator_only_params)
output2u$prpf = output2u$p_detect_lm_182.5_3650 / output2u$n_age_182.5_3650


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
polygon(c(output2l$timestep,rev(output2u$timestep)),
  c(output2l$prpf,output2u$prpf), bty = "n", border = NA,col = adegenet::transp("darkred",0.4))

abline(v = spatial_emanator_timesteps, lty = 2, col = "grey80")
abline(v = bednetstimesteps, lty = 2, lwd = 2, col = "grey30")

## Repeat for incidence
output1comp$inc_6m_10y = 1000*output1comp$p_inc_clinical_182.5_3650/output1comp$n_age_182.5_3650
output1$inc_6m_10y = 1000*output1$p_inc_clinical_182.5_3650/output1$n_age_182.5_3650
output2$inc_6m_10y = 1000*output2$p_inc_clinical_182.5_3650/output2$n_age_182.5_3650
plot(output1comp$inc_6m_10y ~ output1comp$timestep, type = "l", 
     col = "grey80",lwd = 2,
     xlim = c(2000,3400),#ylim = c(0,1),
     ylab = "Clinical incidence per 1,000 children (0.5 - 10 years)",
     xlab = "Date",xaxt = "n", yaxt = "n")
axis(1, at = seq(1,10*365,by = 365),labels = seq(2015,2024,by = 1))
axis(2, las = 2, at = seq(0,10,by = 2))
lines(output1$inc_6m_10y ~ output1$timestep, col = "blue", lwd = 2)
lines(output2$inc_6m_10y ~ output2$timestep, col = "darkred", lwd = 2)
