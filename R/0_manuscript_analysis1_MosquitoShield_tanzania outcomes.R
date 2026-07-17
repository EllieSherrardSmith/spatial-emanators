time = 1:365 ## days
# rstan::traceplot(stan_MS)


eta1uf = eta1f = eta2 = array(dim = c(4000,length(time),2))
trt = c(0,1)
eta3 = array(dim = c(4000,length(time)))

for(tr in 1:length(trt)){
  for(t in 1:length(time)){
    for(i in 1:4000){
      
      eta1uf[i,t,tr] <-MSbtanz$alpha1[i] +  MSbtanz$alpha3[i] *  trt[tr] + MSbtanz$alpha2[i] * time[t]  # vector of draws (length = posterior samples)
      eta1f[i,t,tr] <- MSbtanz$alphaf1[i] + MSbtanz$alphaf3[i] * trt[tr] + MSbtanz$alphaf2[i] * time[t]  # vector of draws (length = posterior samples)
      eta2[i,t,tr] <- MSbtanz$beta1[i] + MSbtanz$beta3[i] * trt[tr] + MSbtanz$beta2[i] * time[t] #+ MSbtanz$beta4[i] * trt[tr] * time[t]  # vector of draws (length = posterior samples)
      eta3[i,t] <- MSbtanz$omega1[i] + MSbtanz$alpha2[i] * time[t]  # vector of draws (length = posterior samples)
      
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
mort_unfed = extract_params_f(MSbtanz$alpha1,MSbtanz$alpha2,MSbtanz$alpha3)
mort_fed = extract_params_f(MSbtanz$alphaf1,MSbtanz$alphaf2,MSbtanz$alphaf3)
successful_fed = extract_params_f(MSbtanz$beta1,MSbtanz$beta2,MSbtanz$beta3)

d1 = as.numeric(quantile(MSbtanz$omega1,c(0.1,0.5,0.9)))
d2 = as.numeric(quantile(MSbtanz$omega2,c(0.1,0.5,0.9)))
dettered = list(d1,d2)

mortuf_pred_mean <- colMedians(mortuf_pred_draws)                # posterior mean
mortf_pred_mean <- colMedians(mortf_pred_draws)                # posterior mean
fed_pred_mean <- colMedians(fed_pred_draws)                # posterior mean
det_pred_mean <- colMedians(det_pred_draws)                # posterior mean

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
## Figure XX2
# par(mfrow = c(2,2),
    # oma = c(0, 0, 2, 0))
# par(mar = c(4,4,1,1))
## Add "A-F"
labels <- LETTERS[5:8]
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

mortuf_pred_mean[1]
mortuf_pred_ci[1,]
# lines(mortf_pred_mean ~ time,ylim = c(0,1),
#       xlim = c(1,40),lwd = 2, lty = 1,col ="lightblue")
probufd = ms_data$duf_t[ms_data$treatments == 1] / ms_data$n_t[ms_data$treatments == 1]
probfd = ms_data$df_t[ms_data$treatments == 1] / ms_data$n_t[ms_data$treatments == 1]
probf = ms_data$fed_t[ms_data$treatments == 1] / ms_data$n_t[ms_data$treatments == 1]
probdet = ms_data$deterrence[ms_data$treatments == 1] / ms_data$deterrence_total[ms_data$treatments == 1]
times_plot = ms_data$time[ms_data$treatments == 1 ]

probufdC = ms_data$duf_t[ms_data$treatments == 0] / ms_data$n_t[ms_data$treatments == 0]
probfdC = ms_data$df_t[ms_data$treatments == 0] / ms_data$n_t[ms_data$treatments == 0]
probfC = ms_data$fed_t[ms_data$treatments == 0] / ms_data$n_t[ms_data$treatments == 0]
probdetC = ms_data$deterrence[ms_data$treatments == 0] / ms_data$deterrence_total[ms_data$treatments == 0]

points(probufd ~ 
         ms_data$time[ms_data$treatments == 1], pch = 20,
       cex = 0.5+sqrt(ms_data$n_t[ms_data$treatments == 1]/
                        mean(ms_data$n_t[ms_data$treatments == 1])),
       col = "darkblue")

points(probufdC ~ 
         ms_data$time[ms_data$treatments == 1], pch = 20,
       cex = 0.5+sqrt(ms_data$n_t[ms_data$treatments == 1]/
                        mean(ms_data$n_t[ms_data$treatments == 1])),
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

probf_new1 <- c(scales::rescale(probf[1:32], to = c(0, 0.699)),rep(0.699,5))
points(probf_new1 ~ 
         ms_data$time[ms_data$treatments == 1], pch = 20, col = "darkred",
       cex = 0.5+sqrt(ms_data$n_t[ms_data$treatments == 1]/
                        mean(ms_data$n_t[ms_data$treatments == 1])))

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

points(probdet ~ 
         ms_data$time[ms_data$treatments == 1], pch = 20, 
       col = "orange",
       cex = 4*sqrt(ms_data$deterrence[ms_data$treatments == 1]/
                      mean(ms_data$deterrence_total[ms_data$treatments == 1])))


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

