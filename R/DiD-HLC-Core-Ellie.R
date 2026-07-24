# --- Environment Setup ---
rm(list = ls(all.names = TRUE))  # Clear all variables from the workspace
gc()  # Run garbage collection to free up memory
if (!is.null(dev.list())) dev.off()  # Close any open graphics devices

library(rstanarm)
library(boot)
library(dplyr)
set.seed(123)
# install.packages("readxl") # CRAN version
df_ento = readxl::read_excel("C:/Users/esherrar/Documents/RProjects/spatial-emanators/data/DiD-data.xlsx", sheet = "A_D4_HumanLandingRate")
df_ento$date = as.Date(df_ento$date, format = "%d-%m-%y")
df_ento$year = lubridate::year(df_ento$date)

# Clusters that were given SEs interventions
target_numbers <- c(3, 12, 20, 42, 55, 57)
replacement_letters <- c("a", "b", "c", "d", "e", "f") 
mapping <- setNames(replacement_letters, target_numbers)
df_ento$site_score <- ifelse(df_ento$cluster_number %in% target_numbers,
                             mapping[as.character(df_ento$cluster_number)],
                             "g")
df_ento$arm <- ifelse(df_ento$cluster_number %in% target_numbers, 1, 0)
df_ento$site_score = relevel(as.factor(df_ento$site_score), ref = "g")
df_cohort1 <- df_ento %>% filter(!study_phase %in% "Cohort 2")
df_cohort1$binary = ifelse(df_cohort1$study_phase == "Baseline", 0, 1)
location = ifelse(df_cohort1$location == "in", 1, 0)
modL2 = data.frame(count_bites = df_cohort1$anoph_count,
                   arm = df_cohort1$arm,
                   binary = df_cohort1$binary,
                   sitea = df_cohort1$cluster_number,
                   household_ID = df_cohort1$household_ID,
                   site_score = df_cohort1$site_score,
                   location = df_cohort1$location)
# modL2 = na.omit(modL2)
###########################################################################
########## 1- did with clusters as a random effect #################
stan2b <- rstanarm::stan_glmer.nb(
   count_bites ~ binary * arm + (1|sitea),#
   data = modL2, link = "log", chains = 4, iter = 2000)

print(stan2b)
plot(stan2b)
plot(stan2b, plotfun = "trace")
summary(stan2b)
prior_summary(stan2b)
pp_check(stan2b, plotfun = "stat", stat = "sd")
pp_check(stan2b, plotfun = "dens_overlay")
yy_pp = posterior_predict(stan2b)
n_sims = nrow(yy_pp)
subset = sample(n_sims, 100)

translate_f = function(pre_or_post, arm){
   cnt = (stan2b$coefficients[1] + 
             stan2b$coefficients[2] * pre_or_post +
             stan2b$coefficients[3] * arm +
            stan2b$coefficients[4] * pre_or_post * arm)
   return(as.numeric(cnt))
   }

B_pre_con = translate_f(pre_or_post = 0, arm = 0)
A1_pre_trt = translate_f(pre_or_post = 0, arm = 1)
D_post_con = translate_f(pre_or_post = 1, arm = 0)
C1_post_trt = translate_f(pre_or_post = 1, arm = 1)

ccv1 = exp((C1_post_trt - A1_pre_trt) - (D_post_con - B_pre_con))
100*(1 - ccv1)

### Generate uncertainty
translate_uncert_f = function(num_reps){
   B_pre_con = A_pre_trt = D_post_con = C_post_trt = numeric(num_reps)
   translate_f = function(pre_or_post, arm, rr){
      cnt = (sims2[rr,1] + 
                sims2[rr,2] * pre_or_post +
                sims2[rr,3] * arm +
                sims2[rr,4] * pre_or_post * arm)
      return(as.numeric(cnt))
      }
   sims2 = as.matrix(stan2b)
   n_sims = nrow(sims2)
   random_draw = sample(n_sims, num_reps)
   for(j in 1:num_reps){
      B_pre_con[j] = translate_f(pre_or_post = 0, arm = 0, rr = random_draw[j])
      A1_pre_trt[j] = translate_f(pre_or_post = 0, arm = 1, rr = random_draw[j])
     
      D_post_con[j] = translate_f(pre_or_post = 1, arm = 0, rr = random_draw[j])
      C1_post_trt[j] = translate_f(pre_or_post = 1, arm = 1, rr = random_draw[j])
   }
   ccc1 = exp((C1_post_trt - A1_pre_trt) - (D_post_con - B_pre_con))
   did1 = 100*(1 - ccc1)
  
   return(data.frame(did1 = did1,
                     B_pre_con = exp(B_pre_con),
                     A1_pre_trt = exp(A1_pre_trt),
                     
                     D_post_con = exp(D_post_con),
                     C1_post_trt = exp(C1_post_trt)
                     ))}
rrr = translate_uncert_f(num_reps = 1000)
quantile(rrr[,1],c(0.05, 0.5, 0.95))

boxplot(sqrt(modL2$count_bites[modL2$binary == 0 & modL2$arm == 0]),
        sqrt(modL2$count_bites[modL2$binary == 0 & modL2$arm == 1]), NA,
        sqrt(modL2$count_bites[modL2$binary == 1 & modL2$arm == 0]),
        sqrt(modL2$count_bites[modL2$binary == 1 & modL2$arm == 1]), 
        yaxt = "n", border = "grey", ylim = sqrt(c(0, 160)),
        col = adegenet::transp(c("red", "blue", NA, "red", "blue"), 0.4),
        ylab = "Mosquito counts", xlab = "Treatment", lwd = 3)
text(x = 1.5, y = 11.5, "(A): HLC-Core", cex = 1)
axis(2, las = 2, at = sqrt(seq(0, 160, 20)), labels = seq(0, 160, 20))
abline(v = 3, lty = 1, col = "grey34", lwd = 3)
mtext("Baseline                        Follow up Cohort 1", 
      cex = 0.8, side = 1, line = 1)
segments(x0 = 1.5, x1 = 4.5,
         y0 = sqrt(as.numeric(quantile(rrr[,2], 0.5))),
         y1 = sqrt(as.numeric(quantile(rrr[,4], 0.5))),
         col = "darkred", lwd = 3)
diff_before = sqrt(as.numeric(quantile(rrr[,2],0.5))) - sqrt(as.numeric(quantile(rrr[,3],0.5)))
segments(x0 = 1.5, x1 = 4.5,
         y0 = sqrt(as.numeric(quantile(rrr[,3],0.5))),
         y1 = sqrt(as.numeric(quantile(rrr[,5],0.5))),
         lty = 1, col = "darkblue", lwd = 3)
y2 = (3*sqrt(as.numeric(quantile(rrr[,5],0.5))) + diff_before) 
segments(x0 = 3, x1 = 4.5, y0 = y2,
         y1 = sqrt(as.numeric(quantile(rrr[,5],0.5))),
         lty = 2, col = "darkblue", lwd = 3)
# points(sqrt(c(as.numeric(quantile(rrr[,7], 0.5)),
#               as.numeric(quantile(rrr[,8], 0.5)),
#               as.numeric(quantile(rrr[,14], 0.5)),
#               as.numeric(quantile(rrr[,15], 0.5)))) ~ c(1,2,4,5), col = c("red", "blue", "red", "blue"), pch = 19) 
# mins_stan4 = sqrt(c(as.numeric(quantile(rrr[,7], 0.05)),
#                     as.numeric(quantile(rrr[,8], 0.05)), NA,
#                     as.numeric(quantile(rrr[,14], 0.05)),
#                     as.numeric(quantile(rrr[,15], 0.05))) )#
# maxs_stan4 = sqrt(c(as.numeric(quantile(rrr[,7], 0.95)),
#                     as.numeric(quantile(rrr[,8], 0.95)), NA,
#                     as.numeric(quantile(rrr[,14], 0.95)),
#                     as.numeric(quantile(rrr[,15], 0.95))) )
# col_segs = c("red", "blue", NA, "red", "blue")
# for(i in c(1,2,4,5)){
#    segments(x0 = i, x1 = i, y0 = mins_stan4[i],
#             y1 = maxs_stan4[i], col = col_segs[i], lwd = 2)}
legend("topright",legend = c("Control","SR"),
       col = c("red", "blue"), pch = 15, bty = "n")

##############################
##
## CDC light trap... 
df_ento = readxl::read_excel("C:/Users/esherrar/Documents/RProjects/spatial-emanators/data/ochomo 2024/DiD-data.xlsx", sheet = "A_D5_LT_anoph_density_Core")
df_ento$date = as.Date(df_ento$date, format = "%d-%m-%y")
df_ento$year = lubridate::year(df_ento$date)

# Clusters that were given SEs interventions
target_numbers <- c(3, 10, 12, 20, 25, 28, 39, 42, 55, 57)
df_ento$arm <- ifelse(df_ento$cluster_number %in% target_numbers, 1, 0)
df_cohort1 <- subset(df_ento, df_ento$study_phase != "Cohort 2")
df_cohort1$binary = ifelse(df_cohort1$study_phase == "Baseline", 0, 1)
modL2 = data.frame(count_bites = df_cohort1$anoph_count,
                   arm = df_cohort1$arm,
                   binary = df_cohort1$binary,
                   sitea = df_cohort1$cluster_number,
                   household_ID = df_cohort1$household_ID)
# modL2 = na.omit(modL2)
###########################################################################
########## 1- did with clusters as a random effect #################
stan2b <- rstanarm::stan_glmer.nb(
  count_bites ~ binary * arm + (1|sitea),#
  data = modL2, link = "log", chains = 4, iter = 2000)

print(stan2b)
plot(stan2b)
plot(stan2b, plotfun = "trace")
summary(stan2b)
prior_summary(stan2b)
pp_check(stan2b, plotfun = "stat", stat = "sd")
pp_check(stan2b, plotfun = "dens_overlay")
yy_pp = posterior_predict(stan2b)
n_sims = nrow(yy_pp)
subset = sample(n_sims, 100)

translate_f = function(pre_or_post, arm){
  cnt = (stan2b$coefficients[1] + 
           stan2b$coefficients[2] * pre_or_post +
           stan2b$coefficients[3] * arm +
           stan2b$coefficients[4] * pre_or_post * arm)
  return(as.numeric(cnt))
}

B_pre_con = translate_f(pre_or_post = 0, arm = 0)
A1_pre_trt = translate_f(pre_or_post = 0, arm = 1)
D_post_con = translate_f(pre_or_post = 1, arm = 0)
C1_post_trt = translate_f(pre_or_post = 1, arm = 1)

ccv1 = exp((C1_post_trt - A1_pre_trt) - (D_post_con - B_pre_con))
100*(1 - ccv1)

### Generate uncertainty
translate_uncert_f = function(num_reps){
  B_pre_con = A_pre_trt = D_post_con = C_post_trt = numeric(num_reps)
  translate_f = function(pre_or_post, arm, rr){
    cnt = (sims2[rr,1] + 
             sims2[rr,2] * pre_or_post +
             sims2[rr,3] * arm +
             sims2[rr,4] * pre_or_post * arm)
    return(as.numeric(cnt))
  }
  sims2 = as.matrix(stan2b)
  n_sims = nrow(sims2)
  random_draw = sample(n_sims, num_reps)
  for(j in 1:num_reps){
    B_pre_con[j] = translate_f(pre_or_post = 0, arm = 0, rr = random_draw[j])
    A1_pre_trt[j] = translate_f(pre_or_post = 0, arm = 1, rr = random_draw[j])
    
    D_post_con[j] = translate_f(pre_or_post = 1, arm = 0, rr = random_draw[j])
    C1_post_trt[j] = translate_f(pre_or_post = 1, arm = 1, rr = random_draw[j])
  }
  ccc1 = exp((C1_post_trt - A1_pre_trt) - (D_post_con - B_pre_con))
  did1 = 100*(1 - ccc1)
  
  return(data.frame(did1 = did1,
                    B_pre_con = exp(B_pre_con),
                    A1_pre_trt = exp(A1_pre_trt),
                    
                    D_post_con = exp(D_post_con),
                    C1_post_trt = exp(C1_post_trt)
  ))}
rrr = translate_uncert_f(num_reps = 1000)
quantile(rrr[,1],c(0.05, 0.5, 0.95))

## repeat for cohrot 2
df_cohort2 <- subset(df_ento, df_ento$study_phase != "Cohort 1")
df_cohort2$binary = ifelse(df_cohort2$study_phase == "Baseline", 0, 1)
modL2 = data.frame(count_bites = df_cohort2$anoph_count,
                   arm = df_cohort2$arm,
                   binary = df_cohort2$binary,
                   sitea = df_cohort2$cluster_number,
                   household_ID = df_cohort2$household_ID)
# modL2 = na.omit(modL2)
###########################################################################
########## 1- did with clusters as a random effect #################
stan2b <- rstanarm::stan_glmer.nb(
  count_bites ~ binary * arm + (1|sitea),#
  data = modL2, link = "log", chains = 4, iter = 2000)

print(stan2b)
plot(stan2b)
plot(stan2b, plotfun = "trace")
summary(stan2b)
prior_summary(stan2b)
pp_check(stan2b, plotfun = "stat", stat = "sd")
pp_check(stan2b, plotfun = "dens_overlay")
yy_pp = posterior_predict(stan2b)
n_sims = nrow(yy_pp)
subset = sample(n_sims, 100)

translate_f = function(pre_or_post, arm){
  cnt = (stan2b$coefficients[1] + 
           stan2b$coefficients[2] * pre_or_post +
           stan2b$coefficients[3] * arm +
           stan2b$coefficients[4] * pre_or_post * arm)
  return(as.numeric(cnt))
}

B_pre_con = translate_f(pre_or_post = 0, arm = 0)
A1_pre_trt = translate_f(pre_or_post = 0, arm = 1)
D_post_con = translate_f(pre_or_post = 1, arm = 0)
C1_post_trt = translate_f(pre_or_post = 1, arm = 1)

ccv1 = exp((C1_post_trt - A1_pre_trt) - (D_post_con - B_pre_con))
100*(1 - ccv1)

### Generate uncertainty
translate_uncert_f = function(num_reps){
  B_pre_con = A_pre_trt = D_post_con = C_post_trt = numeric(num_reps)
  translate_f = function(pre_or_post, arm, rr){
    cnt = (sims2[rr,1] + 
             sims2[rr,2] * pre_or_post +
             sims2[rr,3] * arm +
             sims2[rr,4] * pre_or_post * arm)
    return(as.numeric(cnt))
  }
  sims2 = as.matrix(stan2b)
  n_sims = nrow(sims2)
  random_draw = sample(n_sims, num_reps)
  for(j in 1:num_reps){
    B_pre_con[j] = translate_f(pre_or_post = 0, arm = 0, rr = random_draw[j])
    A1_pre_trt[j] = translate_f(pre_or_post = 0, arm = 1, rr = random_draw[j])
    
    D_post_con[j] = translate_f(pre_or_post = 1, arm = 0, rr = random_draw[j])
    C1_post_trt[j] = translate_f(pre_or_post = 1, arm = 1, rr = random_draw[j])
  }
  ccc1 = exp((C1_post_trt - A1_pre_trt) - (D_post_con - B_pre_con))
  did1 = 100*(1 - ccc1)
  
  return(data.frame(did1 = did1,
                    B_pre_con = exp(B_pre_con),
                    A1_pre_trt = exp(A1_pre_trt),
                    
                    D_post_con = exp(D_post_con),
                    C1_post_trt = exp(C1_post_trt)
  ))}
rrr = translate_uncert_f(num_reps = 1000)
quantile(rrr[,1],c(0.05, 0.5, 0.95))
