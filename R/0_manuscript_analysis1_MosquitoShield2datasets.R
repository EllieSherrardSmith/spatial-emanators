############################################
##
## Cleaned analysis for Kuipou et al. 
## Sarah Moore data 2023 arabiensis and funestus (too few)
## Corine Ngufor data 2026 (two replicates for gambiae)

# --- Environment Setup ---
rm(list = ls(all.names = TRUE))  # Clear all variables from the workspace
gc()  # Run garbage collection to free up memory
if (!is.null(dev.list())) dev.off()  # Close any open graphics devices

library(tidyverse)  # For data manipulation and visualization
library(adegenet)
library(rstan)
library(readr)
library(readxl)     # For reading Excel files
library(lubridate) 

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

## organise data to be read into the model
cleaner_baseA = function(data1,time_cap){
  
  n_t = c(data1$total[data1$treatment_ID == 0],rep(100,5),
          data1$total[data1$treatment_ID == 1],rep(100,5))
  
  duf_t = c(data1$"unfed dead 24"[data1$treatment_ID == 0],rep(0,5),
            data1$"unfed dead 24"[data1$treatment_ID == 1],rep(0,5))
  
  df_t = c(data1$"fed dead 24"[data1$treatment_ID == 0],rep(0,5),
           data1$"fed dead 24"[data1$treatment_ID == 1],rep(0,5))
  
  fed_t = c(data1$"fed alive 24"[data1$treatment_ID == 0],rep(100,5),
            data1$"fed alive 24"[data1$treatment_ID == 1],rep(100,5))
  
  unfed_t = n_t - duf_t - df_t - fed_t 
  
  deterrence = ifelse(c(data1$total[data1$treatment_ID == 0] - data1$total[data1$treatment_ID == 1])<0,0,
                      c(data1$total[data1$treatment_ID == 0] - data1$total[data1$treatment_ID == 1]))

  deterrence_total = c(data1$total[data1$treatment_ID == 0] + data1$total[data1$treatment_ID == 1])

  treatments = c(data1$treatment_ID[data1$treatment_ID == 0],rep(0,5),
                 data1$treatment_ID[data1$treatment_ID == 1],rep(1,5) )
  
  time = c(data1$time_index[data1$treatment_ID == 0],rep(time_cap,5),
           data1$time_index[data1$treatment_ID == 1],rep(time_cap,5))
  
  return(list(  # if trying v3
    N = length(n_t),
    n_t = n_t,
    duf_t = duf_t,
    df_t = df_t,
    fed_t = fed_t,
    unfed_t = unfed_t,
    deterrence = c(deterrence,rep(0,5),deterrence,rep(0,5)),
    deterrence_total = c(deterrence_total,rep(100,5),deterrence_total,rep(100,5)),
    treatments = treatments,
    time = time)) 
}

cleaner_baseB = function(data1, data2, time_cap){
  
  n_t1 = c(data1$total[data1$treatment == 0],rep(100,5),
          data1$total[data1$treatment == 1],rep(100,5))
  
  duf_t1 = c(data1$unf_dead[data1$treatment == 0],rep(0,5),
            data1$unf_dead[data1$treatment == 1],rep(0,5))
  
  df_t1 = c(data1$bf_dead[data1$treatment == 0],rep(0,5),
           data1$bf_dead[data1$treatment == 1],rep(0,5))
  
  fed_t1 = c(data1$bf_live[data1$treatment == 0],rep(100,5),
            data1$bf_live[data1$treatment == 1],rep(100,5))
  
  unfed_t1 = n_t1 - duf_t1 - df_t1 - fed_t1 
  
  deterrence1 = ifelse(c(data1$total[data1$treatment == 0] - data1$total[data1$treatment == 1])<0,0,
                      c(data1$total[data1$treatment == 0] - data1$total[data1$treatment == 1]))

  deterrence_total1 = c(data1$total[data1$treatment == 0] + data1$total[data1$treatment == 1])

  # treated_count1 = data1$total[data1$treatment == 1]
  # control_count1 = ifelse(data1$total[data1$treatment == 0] == 0, 1, data1$total[data1$treatment == 0])
  # 
  treatments1 = c(data1$treatment[data1$treatment == 0],rep(0,5),
                 data1$treatment[data1$treatment == 1],rep(1,5) )
  
  time1 = c(data1$time_index[data1$treatment == 0],rep(time_cap,5),
           data1$time_index[data1$treatment == 1],rep(time_cap,5))
  
  n_t2 = c(data2$total[data2$treatment == 0],rep(100,5),
           data2$total[data2$treatment == 1],rep(100,5))
  
  duf_t2 = c(data2$unf_dead[data2$treatment == 0],rep(0,5),
             data2$unf_dead[data2$treatment == 1],rep(0,5))
  
  df_t2 = c(data2$bf_dead[data2$treatment == 0],rep(0,5),
            data2$bf_dead[data2$treatment == 1],rep(0,5))
  
  fed_t2 = c(data2$bf_live[data2$treatment == 0],rep(100,5),
             data2$bf_live[data2$treatment == 1],rep(100,5))
  
  unfed_t2 = n_t2 - duf_t2 - df_t2 - fed_t2 
  
  deterrence2 = ifelse(c(data2$total[data2$treatment == 0] - data2$total[data2$treatment == 1])<0,0,
                       c(data2$total[data2$treatment == 0] - data2$total[data2$treatment == 1]))
  
  deterrence_total2 = c(data2$total[data2$treatment == 0] + data2$total[data2$treatment == 1])
  
  
  treatments2 = c(data2$treatment[data2$treatment == 0],rep(0,5),
                  data2$treatment[data2$treatment == 1],rep(1,5) )
  
  time2 = c(data2$time_index[data2$treatment == 0],rep(time_cap,5),
            data2$time_index[data2$treatment == 1],rep(time_cap,5))
  
  return(list(  
    N1 = length(n_t1),
    # N1d = length(n_t1)/2,
    n_t1 = n_t1,
    duf_t1 = duf_t1,
    df_t1 = df_t1,
    fed_t1 = fed_t1,
    unfed_t1 = unfed_t1,
    deterrence1 = c(deterrence1,rep(0,5),deterrence1,rep(0,5)),
    deterrence_total1 = c(deterrence_total1,rep(10,5),deterrence_total1,rep(10,5)),
    treatments1 = treatments1,
    time1 = time1,
    
    N2 = length(n_t2),
    # N2d = length(n_t2)/2,
    n_t2 = n_t2,
    duf_t2 = duf_t2,
    df_t2 = df_t2,
    fed_t2 = fed_t2,
    unfed_t2 = unfed_t2,
    deterrence2 = c(deterrence2,rep(0,5),deterrence2,rep(0,5)),
    deterrence_total2 = c(deterrence_total2,rep(10,5),deterrence_total2,rep(10,5)),
    treatments2 = treatments2,
    time2 = time2)) 
}


## Sarah Moore
dtaMS = readxl::read_excel("C:/Users/esherrar/OneDrive - LSTM/Spatial Repellent Review/SRs Systematic Review/data colated.xlsx", sheet = "Mosquito Shield-EHTs" )

dtaMS$date2 = lubridate::ymd(dtaMS$date)
dtaMS <- dtaMS[order(dtaMS$date2), ]   # ensure sorted by date
dtaMS$time_index <- as.integer(dtaMS$date2 - min(dtaMS$date2)) + 1

for(i in 1:nrow(dtaMS)){
  dtaMS$total[i] = dtaMS$`unfed alive 24`[i] + dtaMS$`fed alive 24`[i] +
    dtaMS$`unfed dead 24`[i] + dtaMS$`fed dead 24`[i]
}

dtaMS <- dtaMS[order(dtaMS$treatment_ID, 
                     dtaMS$time_index, 
                     dtaMS$species), ]

dtaMS_summary <- dtaMS %>%
  group_by(treatment_ID, time_index, species) %>%
  summarise(
    across(c("unfed alive 24",
             "fed alive 24",
             "unfed dead 24",
             "fed dead 24",
             total), sum, na.rm = TRUE)
  )


##Ngufor (gambiae)
dtbMS = readxl::read_excel("C:/Users/esherrar/Documents/RProjects/spatial-emanators/data/Ngufor - MosquitoShield EHT/Benin_Mosquito Shield_Hut data_CN.xlsx", sheet = "An. gambiae SMC" )

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

dtbMS1_summary <- dtbMS1 %>%
  group_by(treatment, time_index) %>%
  summarise(
    across(c("unf_live",
             "bf_live",
             "unf_dead",
             "bf_dead",
             total), sum, na.rm = TRUE)
  )

dtbMS2_summary <- dtbMS2 %>%
  group_by(treatment, time_index) %>%
  summarise(
    across(c("unf_live",
             "bf_live",
             "unf_dead",
             "bf_dead",
             total), sum, na.rm = TRUE)
  )
## MODEL 
## data list for stan
ms_data = cleaner_baseA(dtaMS_summary, time_cap = c(40))
ms_datb = cleaner_baseB(dtbMS1_summary, dtbMS2_summary, time_cap = c(40))

## plot to ensure correct
probufd = ms_data$duf_t[ms_data$treatments == 1] / ms_data$n_t[ms_data$treatments == 1]
probfd = ms_data$df_t[ms_data$treatments == 1] / ms_data$n_t[ms_data$treatments == 1]
probf = ms_data$fed_t[ms_data$treatments == 1] / ms_data$n_t[ms_data$treatments == 1]
probdet = ms_data$deterrence[ms_data$treatments == 1] / ms_data$deterrence_total[ms_data$treatments == 1]
times_plot = ms_data$time[ms_data$treatments == 1 ]

probufdB1 = ms_datb$duf_t1[ms_datb$treatments1 == 1] / ms_datb$n_t1[ms_datb$treatments1 == 1]
probfdB1 = ms_datb$df_t1[ms_datb$treatments1 == 1] / ms_datb$n_t1[ms_datb$treatments1 == 1]
probfB1 = ms_datb$fed_t1[ms_datb$treatments1 == 1] / ms_datb$n_t1[ms_datb$treatments1 == 1]
probdetB1 = ms_datb$deterrence1[ms_datb$treatments1 == 1] / ms_datb$deterrence_total1[ms_datb$treatments1 == 1]
times_plotB1 = ms_datb$time1[ms_datb$treatments1 == 1 ]

probufdB2 =    ms_datb$duf_t2[ms_datb$treatments2 == 1] / ms_datb$n_t2[ms_datb$treatments2 == 1]
probfdB2 =     ms_datb$df_t2[ms_datb$treatments2 == 1] / ms_datb$n_t2[ms_datb$treatments2 == 1]
probfB2 =      ms_datb$fed_t2[ms_datb$treatments2 == 1] / ms_datb$n_t2[ms_datb$treatments2 == 1]
probdetB2 =    ms_datb$deterrence2[ms_datb$treatments2 == 1] / ms_datb$deterrence_total2[ms_datb$treatments2 == 1]
times_plotB2 = ms_datb$time2[ms_datb$treatments2 == 1 ]

probufdB1C = ms_datb$duf_t1[ms_datb$treatments1 == 0] / ms_datb$n_t1[ms_datb$treatments1 == 0]
probfdB1C = ms_datb$df_t1[ms_datb$treatments1 == 0] / ms_datb$n_t1[ms_datb$treatments1 == 0]
probfB1C = ms_datb$fed_t1[ms_datb$treatments1 == 0] / ms_datb$n_t1[ms_datb$treatments1 == 0]
probdetB1C = ms_datb$deterrence1[ms_datb$treatments1 == 0] / ms_datb$deterrence_total1[ms_datb$treatments1 == 0]
times_plotB1C = ms_datb$time1[ms_datb$treatments1 ==0  ]

probufdB2C =    ms_datb$duf_t2[ms_datb$treatments2 == 0] / ms_datb$n_t2[ms_datb$treatments2 == 0]
probfdB2C =     ms_datb$df_t2[ms_datb$treatments2 == 0] / ms_datb$n_t2[ms_datb$treatments2 == 0]
probfB2C =      ms_datb$fed_t2[ms_datb$treatments2 == 0] / ms_datb$n_t2[ms_datb$treatments2 == 0]
probdetB2C =    ms_datb$deterrence2[ms_datb$treatments2 == 0] / ms_datb$deterrence_total2[ms_datb$treatments2 == 0]
times_plotB2C = ms_datb$time2[ms_datb$treatments2 == 0]

plot(probdetB1 ~ times_plotB1)

stan_MS <- stan(
  # file = "C:/Users/esherrar/Documents/RProjects/spatial-emanators/stan code/probability_estimates_srs_interaction_v3.stan",
  file = "C:/Users/esherrar/Documents/RProjects/spatial-emanators/stan code/probability_estimates_srs_interaction_v3_multistudy.stan",
  # file = "C:/Users/esherrar/Documents/RProjects/spatial-emanators/stan code/probability_estimates_srs_interaction_v4_multistudy.stan",
  data = ms_datb, 
  warmup = 2000,
  control = list(adapt_delta = 0.98,
                 max_treedepth = 20),
  iter = 3000, chains = 4)

stan_MStanz <- stan(
  file = "C:/Users/esherrar/Documents/RProjects/spatial-emanators/stan code/probability_estimates_srs_interaction_v3.stan",
  # file = "C:/Users/esherrar/Documents/RProjects/spatial-emanators/stan code/probability_estimates_srs_interaction_v3_multistudy.stan",
  # file = "C:/Users/esherrar/Documents/RProjects/spatial-emanators/stan code/probability_estimates_srs_interaction_v4_multistudy.stan",
  data = ms_data,
  warmup = 2000,
  control = list(adapt_delta = 0.98,
                 max_treedepth = 20),
  iter = 3000, chains = 4)

rstan::traceplot(stan_MS,pars = c("alpha1","alpha2","alpha3",
                                  "alphaf1","alphaf2","alphaf3",
                                  "beta1","beta2","beta3",
                                  "omega1","omega2", "phi_det"))

MSb <- rstan::extract(stan_MS)
MSbtanz <- rstan::extract(stan_MStanz)
summary(stan_MS, parameters = c("alpha1","alpha2","alpha3",
                                "alphaf1","alphaf1","alphaf1",
                                "beta1","beta2","beta3",
                                "omega1","omega2"))
# saveRDS(MSb, "C:/Users/esherrar/Documents/RProjects/spatial-emanators/manuscript/simulations/MSb_stan_model_outputNgufor")
# saveRDS(MSbtanz, "C:/Users/esherrar/Documents/RProjects/spatial-emanators/manuscript/simulations/MSbtanz_stan_model_outputTanzania")
time = 1:365 ## days
# rstan::traceplot(stan_MS)


eta1uf = eta1f = eta2 = array(dim = c(4000,length(time),2))
trt = c(0,1)
eta3 = array(dim = c(4000,length(time)))

for(tr in 1:length(trt)){
  for(t in 1:length(time)){
    for(i in 1:4000){
      
      eta1uf[i,t,tr] <-MSb$alpha1[i] +  MSb$alpha3[i] *  trt[tr] + MSb$alpha2[i] * time[t]  # vector of draws (length = posterior samples)
      eta1f[i,t,tr] <- MSb$alphaf1[i] + MSb$alphaf3[i] * trt[tr] + MSb$alphaf2[i] * time[t]  # vector of draws (length = posterior samples)
      eta2[i,t,tr] <- MSb$beta1[i] + MSb$beta3[i] * trt[tr] + MSb$beta2[i] * time[t] #+ MSb$beta4[i] * trt[tr] * time[t]  # vector of draws (length = posterior samples)
      eta3[i,t] <- MSb$omega1[i] + MSb$alpha2[i] * time[t]  # vector of draws (length = posterior samples)
      
    }  
  }
  
}

# inverse logit
mortuf_pred_draws <- plogis(eta1uf[,,2])
mortf_pred_draws <- plogis(eta1f[,,2])                         # inverse logit
fed_pred_draws <- plogis(eta2[,,2])                         # inverse logit
det_pred_draws <- plogis(eta3)                         # inverse logit

# mortuf_pred_drawsC <- plogis(eta1uf[,,1])                         # inverse logit
mortuf_pred_drawsC <- plogis(eta1uf[,,1])
mortf_pred_drawsC <- plogis(eta1f[,,1])                         # inverse logit
fed_pred_drawsC <- plogis(eta2[,,1])                         # inverse logit

## model parameters... 
extract_params_f = function(p1,p2,p3){
  combo = numeric(length = 4000)
  for(i in 1:length(p1)){
    combo[i] = p1[i] + p3[i]
    
  }
  combi = as.numeric(quantile(combo,c(0.1,0.5,0.9)))
  a1 = as.numeric(quantile(p1,c(0.1,0.5,0.9)))
  a2 = as.numeric(quantile(p2,c(0.1,0.5,0.9)))
  ## we want just the coef that indicates difference in outcome due to SE:
  a3 = as.numeric(quantile(p3,c(0.1,0.5,0.9)))
  
  return(list(combi,a1,a2,a3))
}
mort_unfed = extract_params_f(MSb$alpha1,MSb$alpha2,MSb$alpha3)
mort_fed = extract_params_f(MSb$alphaf1,MSb$alphaf2,MSb$alphaf3)
successful_fed = extract_params_f(MSb$beta1,MSb$beta2,MSb$beta3)

d1 = as.numeric(quantile(MSb$omega1,c(0.1,0.5,0.9)))
d2 = as.numeric(quantile(MSb$omega2,c(0.1,0.5,0.9)))
dettered = list(d1,d2)


#########################################
##
## MALARIASIMULATION BELOW
##
##
###########################################
colMedians = function(x){
  X2 = numeric(length = ncol(x))
  for(i in 1:ncol(x)){
    X2[i] = median(x[,i])
  }
  return(X2)
}

mortuf_pred_mean <- colMedians(mortuf_pred_draws)                # posterior mean
mortf_pred_mean <- colMedians(mortf_pred_draws)                # posterior mean
fed_pred_mean <- colMedians(fed_pred_draws)                # posterior mean
det_pred_mean <- colMedians(det_pred_draws)                # posterior mean

# plot(mortuf_pred_mean ~ time, 
#      type = "l",
#      ylim = c(0,1), 
#      xlim = c(1,100), 
#      ylab = "Probability")
# ## model parameters... 
# extract_params_f = function(p1,p2,p3){
#   combo = numeric(length = 4000)
#   for(i in 1:length(p1)){
#     combo[i] = p1[i] + p3[i]
#     
#   }
#   a1 = as.numeric(quantile(combo,c(0.1,0.5,0.9)))
#   a2 = as.numeric(quantile(p2,c(0.1,0.5,0.9)))
#   ## we want just the coef that indicates difference in outcome due to SE:
#   a3 = as.numeric(quantile(p3,c(0.1,0.5,0.9)))
#   
#   return(list(a1,a2))
# }
# mort_unfed = extract_params_f(MSb$alpha1,MSb$alpha2,MSb$alpha3)
# mort_fed = extract_params_f(MSb$alphaf1,MSb$alphaf2,MSb$alphaf3)
# successful_fed = extract_params_f(MSb$beta1,MSb$beta2,MSb$beta3)
# 
# d1 = as.numeric(quantile(MSb$omega1,c(0.1,0.5,0.9)))
# d2 = as.numeric(quantile(MSb$omega2,c(0.1,0.5,0.9)))
# dettered = list(d1,d2)
# 
# 
# test_lineUFD = 1/ (1 + exp(-(mort_unfed[[1]][3] + mort_unfed[[2]][2]*time)))
# test_lineFD = 1/(1 + exp(-(mort_fed[[1]][3] + mort_fed[[2]][2]* time)))
# test_linefed = 0.699 * 1/(1 + exp(-(successful_fed[[1]][3] + successful_fed[[2]][2]* time)))
# js_test = 1 - test_lineFD - test_lineUFD - test_linefed
# test_deterred = 1/(1 + exp(-(dettered[[1]][2] + mort_unfed[[2]][2]* time)))
# 
# plot(test_lineUFD ~ time, xlim = c(1,90), 
#      ylim = c(0,1),
#      type = "l", col = "darkblue")
# lines(c(test_lineFD+test_lineUFD) ~ time, col = "blue")
# lines(c(test_linefed+test_lineFD+test_lineUFD) ~ time, col = "red")
# lines(c(js_test+test_linefed+test_lineFD+test_lineUFD) ~ time, col = "orange")
# 
# re_doing_monthly_UFD = rep(test_lineUFD[1:30],12)
# re_doing_monthly_FD = rep(test_lineFD[1:30],12)
# re_doing_monthly_FED = rep(test_linefed[1:30],12)
# re_doing_monthly_JS = rep(js_test[1:30],12)
# lines(re_doing_monthly_UFD ~ time[1:360], lwd = 2, col = "darkblue")
# lines(c(re_doing_monthly_UFD+re_doing_monthly_FD) ~ time[1:360], lwd = 2, col = "blue")
# lines(c(re_doing_monthly_FED+re_doing_monthly_UFD+re_doing_monthly_FD) ~ time[1:360], lwd = 2, col = "red")
# lines(c(re_doing_monthly_JS+re_doing_monthly_FED+re_doing_monthly_UFD+re_doing_monthly_FD) ~ time[1:360], lwd = 2, col = "orange")
# 
# 
# ## Adjusting for deterrence
# mse_comp <- (1 - test_deterred)
# lUFse_prime <- test_lineUFD * mse_comp
# lFse_prime <- test_lineFD * mse_comp
# kse_prime <- test_linefed * mse_comp
# jse_prime <- js_test * mse_comp + test_deterred
# lse_prime <- (test_lineUFD + test_linefed) * mse_comp
# 
# rse_in =
#   (1 - kse_prime / 0.699) * (jse_prime / (lse_prime + jse_prime))
# 
# sse_in =
#   kse_prime / 0.699
# 
# dfse_in = lFse_prime
# 
# dse_in <- 1 - rse_in - sse_in - dfse_in ## any dead 
# 
# rse_in_comp = 1 - rse_in
# 
# plot(dse_in ~ time, xlim = c(1,90), 
#      ylim = c(0,1),
#      type = "l", col = "darkblue")
# lines(c(dse_in+dfse_in) ~ time, col = "blue")
# lines(c(dse_in+dfse_in+sse_in) ~ time, col = "red")
# lines(c(rse_in+dse_in+dfse_in+sse_in) ~ time, col = "orange")
# 
# 
# phi_indoors = 0.75
# phi_bednets = 0.7
# sse_out = 1
# rse_out_comp = 1 - 0
# rse_out = 0
# rs <- 0
# rs_comp <- 1
# ss <- 1
# sn <- 1
# rn <- 0
# 
# 
# prob_bitten_survives = (
#   (1 - phi_indoors) * sse_out + # yes
#     
#     # The indoor and bed biting should also be multiplied by rse_in_comp (they aren't repelled by se)
#     # rse_in_comp must take into account those that die before biting (- dse)
#     # And if rse + dse + dfse + sse = 1, we need the conditional probability that
#     # a mosquito survives, given that it hasn't been repelled or has already died:
#     # sse_in/(rse_in_comp - dse)
#     (phi_indoors - phi_bednets) * rs_comp * (rse_in_comp - dse_in) * ss * sse_in +
#     phi_bednets * rs_comp * (rse_in_comp - dse_in) * sn * ss * sse_in/(rse_in_comp - dse_in)
#   
# )
# 
# prob_bitten = (
#   (1 - phi_indoors) * rse_out_comp + # yes, where rse_out + dse_out + sse_out = 1 and all that are not repelled will bite
#     
#     # IRS only has the later killing stage, but I think we need to account for the pre-bite deaths here
#     (phi_indoors - phi_bednets) * rs_comp * (rse_in_comp - dse_in) +
#     phi_bednets * rs_comp * (rse_in_comp - dse_in) * sn
#   
#   ## IRS parameters: 
#   ## rs + rs_comp = 1, ds + ss = 1, rs + ds + ss != 1
#   
#   ## Similarly, se parameters: 
#   ## rse + rse_comp = 1, rse + dse + x (survive first stage to bite) = 1, dfse + sse = 1, rse + ds + ss != 1
#   
#   ## If we have rse + dse + dfse + sse = 1, we need to group them carefully
#   ## rse + rse_comp = 1, where rse_comp = dse + dfse + sse
#   ## If a mosquito makes it through to biting, the probability is: 1 - rse - dse = dfse + sse = rse_comp - dse
#   ## Any of these should be fine to use in the conditional probability that you survive, given that you have bitten.
#   ## e.g. sse_in/(1 - rse - dse) = sse_in/(dfse + sse) = sse_in/(rse_comp - dse)
# )
# 
# prob_repelled = (
#   (1 - phi_indoors) * rse_out + # (I've moved this line to the top row to match the other probs)
#     
#     # We need to make sure that we take into account the additive marginal fractional impact
#     # There are three options for indoor repellancy: spray and se, spray-no se, no spray-se
#     ## To write it out in full, we could have:
#     phi_indoors * (rs + rse_in - rs * rse_in) +
#     # which I think will work because rs = 0 when spray is off and rse_in = 0 when se is off
#     
#     # Alternatively, I think this can be replaced with:
#     # phi_indoors * (1 - (rs_comp * rse_in_comp)) +
#     # Which gives the same thing: when spray is off rs_comp = 1, when se is off rse_in_comp = 1
#     
#     phi_bednets * rs_comp * rse_in_comp * rn # (I've swapped rse_in_comp with rn to reflect order)
# )
# 
# plot(prob_bitten_survives ~ time, xlim = c(1,90), 
#      ylim = c(0,1),
#      type = "l", col = "darkblue")
# lines(prob_bitten ~ time, col = "red")
# lines(prob_repelled ~ time, col = "orange")

#########################################
##
## MALARIASIMULATION ABOVE
##
##
###########################################

# 
# mortuf_pred_meanC <- colMedians(mortuf_pred_drawsC)                # posterior mean
# mortf_pred_meanC <- colMedians(mortf_pred_drawsC)                # posterior mean
# fed_pred_meanC <- colMedians(fed_pred_drawsC)                # posterior mean
# colMedians(mortuf_pred_draws)               
mortuf_pred_ci = mortf_pred_ci = fed_pred_ci = det_pred_ci = array(dim = c(length(time),2))
# mortufC_pred_ci = mortfC_pred_ci = fedC_pred_ci = array(dim = c(length(time),2))
for(t in 1:length(time)){
  mortuf_pred_ci[t,] <- as.numeric(quantile(mortuf_pred_draws[,t], c(0.025, 0.975)) )  # 95% CI
  mortf_pred_ci[t,] <- as.numeric(quantile(mortf_pred_draws[,t], c(0.025, 0.975)) )  # 95% CI
  fed_pred_ci[t,] <-  as.numeric(quantile(fed_pred_draws[,t], c(0.025, 0.975)) )  # 95% CI
  det_pred_ci[t,] <-  as.numeric(quantile(det_pred_draws[,t], c(0.025, 0.975)) )  # 95% CI
  
  # mortufC_pred_ci[t,] <- as.numeric(quantile(mortuf_pred_drawsC[,t], c(0.025, 0.975)) )  # 95% CI
  # mortfC_pred_ci[t,] <- as.numeric(quantile(mortf_pred_drawsC[,t], c(0.025, 0.975)) )  # 95% CI
  # fedC_pred_ci[t,] <- as.numeric(quantile(fed_pred_drawsC[,t], c(0.025, 0.975)) )  # 95% CI
  
}

probufdC = ms_data$duf_t[ms_data$treatments == 0] / ms_data$n_t[ms_data$treatments == 0]
probfdC = ms_data$df_t[ms_data$treatments == 0] / ms_data$n_t[ms_data$treatments == 0]
probfC = ms_data$fed_t[ms_data$treatments == 0] / ms_data$n_t[ms_data$treatments == 0]

times_plot = ms_data$time[ms_data$treatments == 1]

## Figure XX2
par(mfrow = c(2,4),
    oma = c(0, 0, 2, 0))
par(mar = c(4,4,1,1))
## Add "A-F"
labels <- LETTERS[1:4]
panel <- 0

## mortality data with model fit... 
panel <- panel + 1
plot(mortuf_pred_mean ~ time, # Unfed dead is darkblue
     ylim = c(0,1),cex.axis=1.3,cex.main=1.3,cex.lab=1.3,
     xlim = c(1,38),type = "l",lwd = 2,
     xlab = "Time in days",yaxt="n",
     ylab = "Probability outcome",col = "darkblue")
axis(2,las=2,at=seq(0,1,0.2),labels = seq(0,1,0.2),cex.axis = 1.3)
text(30,1,"Dead")
for(i in 1:2){
  lines(mortuf_pred_ci[,i] ~ time, lty = 2, lwd = 2,col = "darkblue")
}

mortuf_pred_mean[1] ## initial mortality
mortuf_pred_ci[1,]  ## initial mortality uncertainty

points(probufdB1 ~ 
         ms_datb$time1[ms_datb$treatments1 == 1], pch = 20,
       cex = 0.5+sqrt(ms_datb$n_t1[ms_datb$treatments1 == 1]/
                        mean(ms_datb$n_t1[ms_datb$treatments1 == 1])),
       col = "darkblue")
points(probufdB2 ~ 
         ms_datb$time2[ms_datb$treatments2 == 1], pch = 1,
       cex = 0.5+sqrt(ms_datb$n_t2[ms_datb$treatments2 == 1]/
                        mean(ms_datb$n_t2[ms_datb$treatments2 == 1])),
       col = "darkblue")
points(probufdB1C ~ 
         ms_datb$time1[ms_datb$treatments1 == 0], pch = 20,
       cex = 0.5+sqrt(ms_datb$n_t1[ms_datb$treatments1 == 0]/
                        mean(ms_datb$n_t1[ms_datb$treatments1 == 0])),
       col = "purple")

points(probufdB2C ~ 
         ms_datb$time2[ms_datb$treatments2 == 0], pch = 1,
       cex = 0.5+sqrt(ms_datb$n_t2[ms_datb$treatments2 == 0]/
                        mean(ms_datb$n_t2[ms_datb$treatments2 == 0])),
       col = "purple")

mtext(labels[panel],
      side = 3, line = 0.5, adj = 0,
      font = 2, cex = 1)

panel <- panel + 1
plot(fed_pred_mean ~ time,ylim = c(0,1),
     xlim = c(1,36),type = "l", lwd = 2,
     xlab = "Time in days",yaxt = "n",
     ylab = "Probability outcome",cex.axis=1.3,cex.main=1.3,cex.lab=1.3,
     col = "darkred")
axis(2,las=2,at=seq(0,1,0.2),labels = seq(0,1,0.2),cex.axis = 1.3)
text(30,1,"Success feeding")
for(i in 1:2){
  lines(fed_pred_ci[,i] ~ time, lty = 2, lwd = 2, col = "darkred")
  # lines(fedC_pred_ci[,i] ~ time, lty = 2, lwd = 2, col = "grey")
}

probf_new1 <- c(scales::rescale(probfB1[1:32], to = c(0, 0.699)),rep(0.699,5))
probf_new2 <- c(scales::rescale(probfB2[1:32], to = c(0, 0.699)),rep(0.699,5))
points(probf_new1 ~ 
         ms_datb$time1[ms_datb$treatments1 == 1], pch = 20, col = "darkred",
       cex = 0.5+sqrt(ms_datb$n_t1[ms_datb$treatments1 == 1]/
                        mean(ms_datb$n_t1[ms_datb$treatments1 == 1])))
points(probf_new2 ~ 
         ms_datb$time2[ms_datb$treatments2 == 1], pch = 1, col = "darkred",
       cex = 0.5+sqrt(ms_datb$n_t2[ms_datb$treatments2 == 1]/
                        mean(ms_datb$n_t2[ms_datb$treatments2 == 1])))

# max(probfC[1:32])
# probfC_new <- c(scales::rescale(probfC[1:32], to = c(0, 0.699)),rep(0.699,5))
# points(probfC_new ~
#          ms_data$time[ms_data$treatments == 0],
#        pch = 1,
#        col = adegenet::transp("grey24",0.4),
#        cex = 0.5+sqrt(ms_data$n_t[ms_data$treatments == 0]/
#                  mean(ms_data$n_t[ms_data$treatments == 0])))
mtext(labels[panel],
      side = 3, line = 0.5, adj = 0,
      font = 2, cex = 1)

panel <- panel + 1
plot(det_pred_mean ~ time,ylim = c(0,1),
     xlim = c(1,36),type = "l",cex.axis=1.3,cex.main=1.3,cex.lab=1.3,
     xlab = "Time in days",yaxt = "n",
     ylab = "Probability outcome",
     col = "orange")
axis(2,las=2,at=seq(0,1,0.2),labels = seq(0,1,0.2),cex.axis = 1.3)
text(30,1,"Deterrence")
for(i in 1:2){
  lines(det_pred_ci[,i] ~ time, lty = 2, lwd = 2, col = "orange")
}

points(probdetB1 ~ 
         ms_datb$time1[ms_datb$treatments1 == 1], pch = 20, 
       col = "orange",
       cex = 0.5*sqrt(ms_datb$deterrence1[ms_datb$treatments1 == 1]/
                      mean(ms_datb$deterrence_total1[ms_datb$treatments1 == 1])))

points(probdetB2 ~ 
         ms_datb$time1[ms_datb$treatments1 == 1], pch = 1, 
       col = "orange",
       cex = 4*sqrt(ms_datb$deterrence1[ms_datb$treatments1 == 1]/
                      mean(ms_datb$deterrence_total1[ms_datb$treatments1 == 1])))

mtext(labels[panel],
      side = 3, line = 0.5, adj = 0,
      font = 2, cex = 1)



### Combining with deterrence to estimates model parameters
## indoor use
mortf <- mortf_pred_mean
mortany <- mortuf_pred_mean + mortf
fed1 <- (1 - mortany) * (fed_pred_mean) * (1 - det_pred_mean)
dead1 <- mortany * (1 - det_pred_mean)
dead1fed <- mortf * (1 - det_pred_mean)
rep1 <- (1 - (dead1 + fed1)) * (1 - det_pred_mean)
det1 = det_pred_mean
# probability_plot
TOTS = dead1 + rep1 + fed1 + det1

first_line = fed1 / TOTS
second_line = (fed1 + rep1) / TOTS
third_line = (fed1 + rep1 + det1 ) / TOTS
fourth_line = (fed1 + rep1 + det1 + dead1fed) / TOTS

Time = 1:365
Time2 = rev(Time)
minimal = rep(0,length(Time))
maximal = rep(1,length(Time))

# panel <- panel + 1
# plot(rev(first_line) ~ Time2,ylim=c(0,1),yaxt="n",
#      ylab= "Probability outcomes",xaxt="n",
#      main="",xlab = "Time in days",
#      cex.axis=1.3,cex.main=1.3,cex.lab=1.3,bty="n",pch="")
# axis(2,las=2,at=seq(0,1,0.2),labels = seq(0,1,0.2),cex.axis = 1.3)
# axis(1,at=seq(0,365,120),labels = seq(0,365,120),cex.axis = 1.3)
# polygon(c(Time2,rev(Time2)),c(rev(first_line),rev(minimal)),col=transp("red",0.4),border = NA)
# polygon(c(Time2,rev(Time2)),c(rev(second_line),first_line),col=transp("orange",0.6),border = NA)
# polygon(c(Time2,rev(Time2)),c(rev(third_line),second_line),col=transp("darkgreen",0.4),border = NA)
# polygon(c(Time2,rev(Time2)),c(maximal,third_line),col=transp("blue",0.3),border = NA)
# polygon(c(Time2,rev(Time2)),c(maximal,fourth_line),col=transp("royalblue",0.6),border = NA)
# abline(v = 28,lty = 2)
# 
# mtext(labels[panel],
#       side = 3, line = 0.5, adj = 0,
#       font = 2, cex = 1)

##adjust for resetting MS every 28 days as per Kenya trial

panel <- panel + 1
first_line = rep(first_line[1:28],13)
second_line = rep(second_line[1:28],13)
third_line = rep(third_line[1:28],13)
fourth_line = rep(fourth_line[1:28],13)

first_line = first_line[1:365] 
second_line = second_line[1:365] 
third_line = third_line[1:365] 
fourth_line = fourth_line[1:365] 

Time = 1:365
Time2 = rev(Time)
minimal = rep(0,length(Time))
maximal = rep(1,length(Time))

plot(rev(first_line) ~ Time2,ylim=c(0,1),yaxt="n",
     ylab= "Probability outcomes",xaxt="n",
     main="",xlab = "Time in days",
     cex.axis=1.3,cex.main=1.3,cex.lab=1.3,bty="n",pch="")
axis(2,las=2,at=seq(0,1,0.2),labels = seq(0,1,0.2),cex.axis = 1.3)
axis(1,at=seq(0,365,120),labels = seq(0,365,120),cex.axis = 1.3)
polygon(c(Time2,rev(Time2)),c(rev(first_line),rev(minimal)),col=transp("red",0.4),border = NA)
polygon(c(Time2,rev(Time2)),c(rev(second_line),first_line),col=transp("orange",0.6),border = NA)
polygon(c(Time2,rev(Time2)),c(rev(third_line),second_line),col=transp("darkgreen",0.4),border = NA)
polygon(c(Time2,rev(Time2)),c(maximal,third_line),col=transp("blue",0.3),border = NA)
polygon(c(Time2,rev(Time2)),c(maximal,fourth_line),col=transp("royalblue",0.6),border = NA)
abline(v = seq(1,365, by = 28),lty = 2)
mtext(labels[panel],
      side = 3, line = 0.5, adj = 0,
      font = 2, cex = 1)

## JUMP TO 0_manuscript_analysis1_MosquitoShield_tanzania outcomes.R to complete figure 2


##
## Simulate the Kenta_simulations_v1.R code
## lines 239 to 365 (plus devtools bits!)
## Then loop through a bunch of random draws from the statistical model
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

pyr_only_res_20 = read.csv("C:/Users/esherrar/Documents/RProjects/spatial-emanators/data/ochomo 2024/pyr_only_res_20.csv")
pyr_only_res_50 = read.csv("C:/Users/esherrar/Documents/RProjects/spatial-emanators/data/ochomo 2024/pyr_only_res_50.csv")
pyr_pbo_res_80 = read.csv("C:/Users/esherrar/Documents/RProjects/spatial-emanators/data/ochomo 2024/pyr_PBO_res_80.csv")
# MSb = readRDS("C:/Users/esherrar/Documents/RProjects/spatial-emanators/manuscript/simulations/MSb_stan_model_output")

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
det_a2 = sts$alpha2[rr]

Tfed_a1 = stsT$beta1[rr] + stsT$beta3[rr]
Tfed_a2 = stsT$beta2[rr]
Tfeddead_a1 = stsT$alphaf1[rr] + stsT$alphaf3[rr]
Tfeddead_a2 = stsT$alphaf2[rr]
Tunfeddead_a1 = stsT$alpha1[rr] + stsT$alpha3[rr]
Tunfeddead_a2 =  stsT$alpha2[rr]
Tdet_a1 = stsT$omega1[rr]
Tdet_a2 = stsT$alpha2[rr]

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
dtsens = data.frame(reps = 1:50,
                    fed_a1 = rep(NA,50),
                    fed_a2 = rep(NA,50),
                    feddead_a1 = rep(NA,50),
                    feddead_a2 = rep(NA,50),
                    unfeddead_a1 = rep(NA,50),
                    unfeddead_a2 = rep(NA,50),
                    det_a1 = rep(NA,50),
                    det_a2 = rep(NA,50),
                    
                    Tfed_a1 = rep(NA,50),
                    Tfed_a2 = rep(NA,50),
                    Tfeddead_a1 = rep(NA,50),
                    Tfeddead_a2 = rep(NA,50),
                    Tunfeddead_a1 = rep(NA,50),
                    Tunfeddead_a2 = rep(NA,50),
                    Tdet_a1 = rep(NA,50),
                    Tdet_a2 = rep(NA,50))
for(i in 1:50){
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
for(i in 1:50){
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
  
  i_draw <- sample(x = 1:4000, size = 1, replace = FALSE)  # or sample a random draw if you like
  
  spatial_emanator_timesteps = round(seq(310+1095+1095, 310+1095+1095+730, by = 30),0) # Each round of Mosquito Shield spatial emanator is implemented in the 1st and second year.
  
  spatial_emanator_only_params <- set_spatial_emanator(
    bednetparams,
    timesteps = spatial_emanator_timesteps,
    coverages = rep(0.8, length(spatial_emanator_timesteps)), # # The first round covers 30% of the population and the second covers 80%. 
    kse_in_a1 = matrix(c(rep(dtsens$Tbeta1[i_draw],length(spatial_emanator_timesteps)),
                         rep(dtsens$beta1[i_draw],length(spatial_emanator_timesteps)),
                         rep(dtsens$beta1[i_draw],length(spatial_emanator_timesteps))), ncol=3), # Matrix of impact parameters indoors
    kse_in_a2 = matrix(c(rep(dtsens$Tbeta2[i_draw],length(spatial_emanator_timesteps)),
                         rep(dtsens$beta2[i_draw],length(spatial_emanator_timesteps)),
                         rep(dtsens$beta2[i_draw],length(spatial_emanator_timesteps))),ncol=3), # Matrix of impact parameters indoors
    dfse_in_a1 = matrix(c(rep(dtsens$Tfeddead_a1[i_draw],length(spatial_emanator_timesteps)),
                          rep(dtsens$feddead_a1[i_draw],length(spatial_emanator_timesteps)),
                          rep(dtsens$feddead_a1[i_draw],length(spatial_emanator_timesteps))),ncol=3), # Matrix of impact parameters indoors
    dfse_in_a2 = matrix(c(rep(dtsens$Tfeddead_a2[i_draw],length(spatial_emanator_timesteps)),
                          rep(dtsens$feddead_a2[i_draw],length(spatial_emanator_timesteps)),
                          rep(dtsens$feddead_a2[i_draw],length(spatial_emanator_timesteps))),ncol=3), # Matrix of impact parameters indoors
    dufse_in_a1 = matrix(c(rep(dtsens$Tunfeddead_a1[i_draw],length(spatial_emanator_timesteps)),
                           rep(dtsens$unfeddead_a1[i_draw],length(spatial_emanator_timesteps)),
                           rep(dtsens$unfeddead_a1[i_draw],length(spatial_emanator_timesteps))), ncol=3), # Matrix of impact parameters indoors
    dufse_in_a2 = matrix(c(rep(dtsens$Tunfeddead_a2[i_draw],length(spatial_emanator_timesteps)),
                           rep(dtsens$unfeddead_a2[i_draw],length(spatial_emanator_timesteps)),
                           rep(dtsens$unfeddead_a2[i_draw],length(spatial_emanator_timesteps))),ncol=3), # Matrix of impact parameters indoors
    detse_in_a1 = matrix(c(rep(dtsens$Tdet_a1[i_draw],length(spatial_emanator_timesteps)),
                           rep(dtsens$det_a1[i_draw],length(spatial_emanator_timesteps)),
                           rep(dtsens$det_a1[i_draw],length(spatial_emanator_timesteps))), ncol=3), # Matrix of impact parameters indoors
    detse_in_a2 = matrix(c(rep(dtsens$Tdet_a2[i_draw],length(spatial_emanator_timesteps)),
                           rep(dtsens$det_a2[i_draw],length(spatial_emanator_timesteps)),
                           rep(dtsens$det_a2[i_draw],length(spatial_emanator_timesteps))),ncol=3) # Matrix of impact parameters indoors
    
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
  
  print(i_draw)
}
# saveRDS(store1,"C:/Users/esherrar/Documents/RProjects/spatial-emanators/manuscript/simulations/store1")
# saveRDS(store_inc1,"C:/Users/esherrar/Documents/RProjects/spatial-emanators/manuscript/simulations/store_inc1")
# saveRDS(store_ara1,"C:/Users/esherrar/Documents/RProjects/spatial-emanators/manuscript/simulations/store_ara1")
# saveRDS(store_gam1,"C:/Users/esherrar/Documents/RProjects/spatial-emanators/manuscript/simulations/store_gam1")
# saveRDS(store_fun1,"C:/Users/esherrar/Documents/RProjects/spatial-emanators/manuscript/simulations/store_fun1")
# saveRDS(store_eir1gam,"C:/Users/esherrar/Documents/RProjects/spatial-emanators/manuscript/simulations/store_eir1gam")
# saveRDS(store_eir1ara,"C:/Users/esherrar/Documents/RProjects/spatial-emanators/manuscript/simulations/store_eir1ara")
# saveRDS(store_eir1fun,"C:/Users/esherrar/Documents/RProjects/spatial-emanators/manuscript/simulations/store_eir1fun")
# saveRDS(store2,"C:/Users/esherrar/Documents/RProjects/spatial-emanators/manuscript/simulations/store2")
# saveRDS(store_inc2,"C:/Users/esherrar/Documents/RProjects/spatial-emanators/manuscript/simulations/store_inc2")
# saveRDS(store_ara2,"C:/Users/esherrar/Documents/RProjects/spatial-emanators/manuscript/simulations/store_ara2")
# saveRDS(store_gam2,"C:/Users/esherrar/Documents/RProjects/spatial-emanators/manuscript/simulations/store_gam2")
# saveRDS(store_fun2,"C:/Users/esherrar/Documents/RProjects/spatial-emanators/manuscript/simulations/store_fun2")
# saveRDS(store_eir2gam,"C:/Users/esherrar/Documents/RProjects/spatial-emanators/manuscript/simulations/store_eir2gam")
# saveRDS(store_eir2ara,"C:/Users/esherrar/Documents/RProjects/spatial-emanators/manuscript/simulations/store_eir2ara")
# saveRDS(store_eir2fun,"C:/Users/esherrar/Documents/RProjects/spatial-emanators/manuscript/simulations/store_eir2fun")

store1 = readRDS("C:/Users/esherrar/Documents/RProjects/spatial-emanators/manuscript/simulations/store1")
store_inc1 = readRDS("C:/Users/esherrar/Documents/RProjects/spatial-emanators/manuscript/simulations/store_inc1")
store_ara1 = readRDS("C:/Users/esherrar/Documents/RProjects/spatial-emanators/manuscript/simulations/store_ara1")
store_gam1 = readRDS("C:/Users/esherrar/Documents/RProjects/spatial-emanators/manuscript/simulations/store_gam1")
store_fun1 = readRDS("C:/Users/esherrar/Documents/RProjects/spatial-emanators/manuscript/simulations/store_fun1")
store_eir1gam = readRDS("C:/Users/esherrar/Documents/RProjects/spatial-emanators/manuscript/simulations/store_eir1gam")
store_eir1ara = readRDS("C:/Users/esherrar/Documents/RProjects/spatial-emanators/manuscript/simulations/store_eir1ara")
store_eir1fun = readRDS("C:/Users/esherrar/Documents/RProjects/spatial-emanators/manuscript/simulations/store_eir1fun")

store2 = readRDS("C:/Users/esherrar/Documents/RProjects/spatial-emanators/manuscript/simulations/store2")
store_inc2 = readRDS("C:/Users/esherrar/Documents/RProjects/spatial-emanators/manuscript/simulations/store_inc2")
store_ara2 = readRDS("C:/Users/esherrar/Documents/RProjects/spatial-emanators/manuscript/simulations/store_ara2")
store_gam2 = readRDS("C:/Users/esherrar/Documents/RProjects/spatial-emanators/manuscript/simulations/store_gam2")
store_fun2 = readRDS("C:/Users/esherrar/Documents/RProjects/spatial-emanators/manuscript/simulations/store_fun2")
store_eir2gam = readRDS("C:/Users/esherrar/Documents/RProjects/spatial-emanators/manuscript/simulations/store_eir2gam")
store_eir2ara = readRDS("C:/Users/esherrar/Documents/RProjects/spatial-emanators/manuscript/simulations/store_eir2ara")
store_eir2fun = readRDS("C:/Users/esherrar/Documents/RProjects/spatial-emanators/manuscript/simulations/store_eir2fun")


plot(store_ara1[,2] ~ output1$timestep, type = "l", 
     col = "white",lwd = 2,
     xlim = c(2000,3400), ylim = c(0,400000),
     ylab = "Mosquito densities",
     xlab = "Date",xaxt = "n", yaxt = "n")
axis(1, at = seq(1,10*365,by = 365),labels = seq(2015,2024,by = 1))
# axis(2, las = 2, at = seq(0,400000,by = 100000),labels = seq(0,400,by = 100))
# lines(output1$prpf ~ output1$timestep, col = "blue", lwd = 2)

for(i in 1:50){
  
  lines(store_ara1[,i+1] ~ output1$timestep, col = adegenet::transp("darkblue",0.2), lwd = 2)
  lines(store_ara2[,i+1] ~ output1$timestep, col = adegenet::transp("darkred",0.2), lwd = 2)
  lines(store_gam1[,i+1] ~ output0$timestep, col = adegenet::transp("blue",0.2), lwd = 1)
  lines(store_gam2[,i+1] ~ output0$timestep, col = adegenet::transp("red",0.2), lwd = 1)
  lines(store_fun1[,i+1] ~ output0$timestep, col = adegenet::transp("lightblue",0.2), lwd = 1)
  lines(store_fun2[,i+1] ~ output0$timestep, col = adegenet::transp("pink",0.2), lwd = 1)
  
}
abline(v = spatial_emanator_timesteps[1], lty = 2, col = "grey80")
abline(v = spatial_emanator_timesteps[1]+365, lty = 2, col = "grey80")
# abline(v = spatial_emanator_timesteps[1]+365+365, lty = 2, col = "grey80")
abline(v = bednetstimesteps, lty = 2, col = "grey30")


## Prepare figures 
## Final analysis for MS paper... 
## Figure XX2
par(mfrow = c(1,3),
    oma = c(0, 0, 2, 0))
par(mar = c(8,5,1,1))
## Add "A-F"
labels <- LETTERS[1:3]
panel <- 0

panel <- panel + 1


## for A plot the relative reduction in mosquitoes
summary_fun = function(store_objectA,
                       store_objectB){
  fc1 = ft1 = array(dim = c(365,3))
  
  for(t in 2501:2865){
    fc1[t-2500,] = quantile(as.numeric(store_objectA[t,2:51]),c(0.05,0.5,0.95)) 
    ft1[t-2500,] = quantile(as.numeric(store_objectB[t,2:51]),c(0.05,0.5,0.95))
  }
  rel_low = (fc1[1:365,2] - ft1[1:365,1])/fc1[1:365,2]
  rel_mean = (fc1[1:365,2] - ft1[1:365,2])/fc1[1:365,2]
  rel_upp = (fc1[1:365,2] - ft1[1:365,3])/fc1[1:365,2]
  
  return(c(
    mean(rel_mean),
    mean(rel_upp),
    mean(rel_low)
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
    fc1[t-2865,] = quantile(as.numeric(store_objectA[t,2:51]),c(0.05,0.5,0.95)) 
    ft1[t-2865,] = quantile(as.numeric(store_objectB[t,2:51]),c(0.05,0.5,0.95))
  }
  rel_low = (fc1[1:365,2] - ft1[1:365,1])/fc1[1:365,2]
  rel_mean = (fc1[1:365,2] - ft1[1:365,2])/fc1[1:365,2]
  rel_upp = (fc1[1:365,2] - ft1[1:365,3])/fc1[1:365,2]
  
  return(c(
    mean(rel_mean),
    mean(rel_upp),
    mean(rel_low)
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
    fc1[t-2500,] = quantile(as.numeric(store_objectA[t,2:51]),c(0.05,0.5,0.95)) 
    ft1[t-2500,] = quantile(as.numeric(store_objectB[t,2:51]),c(0.05,0.5,0.95))
  }
  rel_low = (fc1[1:730,2] - ft1[1:730,1])/fc1[1:730,2]
  rel_mean = (fc1[1:730,2] - ft1[1:730,2])/fc1[1:730,2]
  rel_upp = (fc1[1:730,2] - ft1[1:730,3])/fc1[1:730,2]
  
  return(c(
    mean(rel_mean),
    mean(rel_upp),
    mean(rel_low)
  ))
}
summary_fun_both(store1,store2)
summary_fun_both(store_inc1,store_inc2)
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

# manuscript = c(11.5 (-41.6, 112.9))

plot(c(0.115,NA,MEAN_EMPIRICL_results[1:2],reductions[1:3],NA,
       MEAN_EMPIRICL_results[3],  reductions[4:6]) ~ c(1:12),
     ylim = c(-0.2,0.8),
     ylab = "% reduction in mosquito density",
     yaxt = "n",
     xaxt = "n", xlab = "",
     pch = c(19,NA,19,19,16,17,18,NA,19,16,17,18),
     cex = 1.5, 
     col = c("grey30",NA,"blue", "blue", rep("darkred",3),"white", "blue", rep("darkred",3)))


abline(v = 8, lty = 2, col = "grey50", lwd = 2)
abline(v = 2, lty = 2, col = "grey50", lwd = 2)
# text(1,0.78,"Trial")
text(5,0.78,"Follow up year 1")
text(10,0.78,"Follow up year 2")
axis(2, las = 2, at = seq(-0.2,0.8,by = 0.2),
     labels = seq(-20,80,by = 20))

lower_list = c(-0.416,NA,low_EMPIRICL_results[1:2],lower[1:3],NA,
               low_EMPIRICL_results[3],  lower[4:6])
upper_list = c(1.129,NA,upp_EMPIRICL_results[1:2],upper[1:3],NA,
               upp_EMPIRICL_results[3],  upper[4:6])
cols_list =  col = c("grey30",NA,"blue", "blue", rep("darkred",3),"white", "blue", rep("darkred",3))
for(i in 1:12){
  segments(x0 = i, x1 = i,
           y0 = lower_list[i],
           y1 = upper_list[i],
           lty = 2, col = cols_list[i])}

labs <- c("Core region overall","","HLC","light traps","simulated gambiae","arabiensis","funestus","",
          "light traps","simulated gambiae","arabiensis","funestus")

x <- seq(1, 12, by = 1)

# draw axis without labels
axis(1, at = x, labels = FALSE)

# add rotated labels
text(x = x,
     y = par("usr")[3] - 0.05 * diff(par("usr")[3:4]),
     labels = labs,
     srt = 45,
     adj = 1,
     xpd = TRUE)
axis(2, las = 2, at = seq(0,12,by = 2))

legend("bottomright",
       legend = c("Trial reported",
                  "DiD-estimate",
                  "Model simulation"),
       col = c("grey30","blue","darkred"),
       bty = "n",
       pch = 19)


mtext(labels[panel],
      side = 3, line = 0.5, adj = 0,
      font = 2, cex = 1)

#
panel <- panel + 1
plot(output1$prpf ~ output1$timestep, type = "l", 
     col = "white",lwd = 2,
     xlim = c(2000,3400), ylim = c(0,0.6),
     ylab = "Model simulated estimate for prevalence in 0.5 - 10 years (%)",
     xlab = "",xaxt = "n", yaxt = "n")
axis(1, at = seq(1,10*365,by = 365),labels = seq(2015,2024,by = 1))
axis(2, las = 2, at = seq(0,1,by = 0.2),labels = seq(0,100,by = 20))
# lines(output1$prpf ~ output1$timestep, col = "blue", lwd = 2)
abline(v = spatial_emanator_timesteps, lty = 2, col = "grey80")
# abline(v = spatial_emanator_timesteps[1]+365, lty = 2, col = "grey80")
# abline(v = spatial_emanator_timesteps[1]+365+365, lty = 2, col = "grey80")
abline(v = bednetstimesteps, lty = 2, col = "grey30")

for(i in 1:50){
  
  lines(store1[,i+1] ~ output1$timestep, col = adegenet::transp("blue",0.2), lwd = 2)
  lines(store2[,i+1] ~ output0$timestep, col = adegenet::transp("darkred",0.2), lwd = 1)
  
}

mtext(labels[panel],
      side = 3, line = 0.5, adj = 0,
      font = 2, cex = 1)

## prevalence reductions
quantile(unlist(store1[2501:2865,]),c(0.025,0.5,0.975))
quantile(unlist(store2[2501:2865,]),c(0.025,0.5,0.975))

(quantile(unlist(store1[2501:2865,]),0.5)-
    quantile(unlist(store2[2501:2865,]),c(0.1,0.5,0.9)))/
  quantile(unlist(store1[2501:2865,]),0.5)

(quantile(unlist(store1[2866:3230,]),0.5)-
    quantile(unlist(store2[2866:3230,]),c(0.1,0.5,0.9)))/
  quantile(unlist(store1[2866:3230,]),0.5)

(quantile(unlist(store1[2501:3230,]),0.5)-
    quantile(unlist(store2[2501:3230,]),c(0.1,0.5,0.9)))/
  quantile(unlist(store1[2501:3230,]),0.5)

panel <- panel + 1
# par(mar = c(6,4,1,1))
## yearly incidence estimates
## baseline year 2136:2500
## year 1 follow up 2501:2865
## year 2 follow up 2866:323
baseC = c(store_inc1[2136:2500,i+1])
baseT = c(store_inc2[2136:2500,i+1])
FOL_C1 = c(store_inc1[2501:2865,i+1])
FOL_T1 = c(store_inc2[2501:2865,i+1])
FOL_C2 = c(store_inc1[2866:3230,i+1])
FOL_T2 = c(store_inc2[2866:3230,i+1])

(quantile(unlist(store_inc1[2501:2865,]),0.5) -
    quantile(unlist(store_inc2[2501:2865,]),c(0.1,0.5,0.9)))/
  quantile(unlist(store_inc1[2501:2865,]),0.5)

(quantile(unlist(store_inc1[2866:3230,]),0.5) -
    quantile(unlist(store_inc2[2866:3230,]),c(0.1,0.5,0.9)))/
  quantile(unlist(store_inc1[2866:3230,]),0.5)

boxplot(baseC,baseT,NA,FOL_C1,FOL_T1,NA,FOL_C2,FOL_T2, 
        col = c("grey80","grey30","white",
                "darkblue","darkred","white",
                "darkblue","darkred"),
        # xlim = c(2000,3400),#ylim = c(0,1),
        ylab = "Model simulated clinical incidence per 1,000 children (0.5 - 10 years)",
        xlab = "",xaxt = "n", yaxt = "n")
labs <- c("Baseline control","Baseline SE","",
          "Year 1 control","Year 1 SE","","Year 2 control","Year 2 SE")

x <- seq(1, 8, by = 1)

# draw axis without labels
axis(1, at = x, labels = FALSE)

# add rotated labels
text(x = x,
     y = par("usr")[3] - 0.05 * diff(par("usr")[3:4]),
     labels = labs,
     srt = 45,
     adj = 1,
     xpd = TRUE)
axis(2, las = 2, at = seq(0,10,by = 2))


mtext(labels[panel],
      side = 3, line = 0.5, adj = 0,
      font = 2, cex = 1)


