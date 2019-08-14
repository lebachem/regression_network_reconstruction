#-----------------------------------------------------------#
# File that creates tables for comparing L_1 and L_2 errors #
#-----------------------------------------------------------#

# remove the old stuff
rm(list=ls())


library(stargazer)

# set working directory----
setwd("/Workflow/")
load("1) Data Read In/data.RData")

# Load the IPFP Model----

load("3) Model/IPFP Model/results/ipfp.RData")
r1<-rowbuilder(X_hat_ipfp,X_inform)$mean_row

# load the Regression model with random Effects----

X_hat_raneff<-c()
for (tau in 1:length(t)){

load(paste("3) Model/Regression with random effects/results/res",tau,".RData",sep=""))
 X_hat_raneff<-cbind(X_hat_raneff,esol) 
}  

r2<-rowbuilder(X_hat_raneff,X_inform)$mean_row

# load the Regression model with gdp_i, gdp_j, trade_ij----

X_hat_gdp_trade<-c()
for (tau in 1:length(t)){
  
  load(paste("3) Model/Regression with gdp and trade/results/res",tau,".RData",sep=""))
  X_hat_gdp_trade<-cbind(X_hat_gdp_trade,esol) 
}  

r3<-rowbuilder(X_hat_gdp_trade,X_inform)$mean_row

# load the Regression model with gdp_i, gdp_j, trade_ij and random effects----

X_hat_gdp_trade_raneff<-c()
for (tau in 1:length(t)){
  
  load(paste("3) Model/Regression with gdp and trade and random effects/results/res",tau,".RData",sep=""))
  X_hat_gdp_trade_raneff<-cbind(X_hat_gdp_trade_raneff,esol) 
}  

r4<-rowbuilder(X_hat_gdp_trade_raneff,X_inform)$mean_row

# load the Regression model with gdp_i, gdp_j, dist_ij----

X_hat_gdp_dist<-c()
for (tau in 1:length(t)){
  
  load(paste("3) Model/Regression with gdp and distance/results_without_bootstrap/res",tau,".RData",sep=""))
  X_hat_gdp_dist<-cbind(X_hat_gdp_dist,esol) 
}  

r4<-rowbuilder(X_hat_gdp_dist,X_inform)$mean_row

# load the Regression model with gdp_i, gdp_j, dist_ij and random effects----

X_hat_gdp_dist_raneff<-c()
for (tau in 1:length(t)){
  
  load(paste("3) Model/Regression with gdp and distance and random effects/results/res",tau,".RData",sep=""))
  X_hat_gdp_dist_raneff<-cbind(X_hat_gdp_dist_raneff,esol) 
}  

r5<-rowbuilder(X_hat_gdp_dist_raneff,X_inform)$mean_row

# load the Regression model with gdp_i, gdp_j and lagged values----

X_hat_gdp_lag<-c()
for (tau in 1:length(t)){
  
  load(paste("3) Model/Regression with gdp and lagged values/results_without_bootstrap/res",tau,".RData",sep=""))
  X_hat_gdp_lag<-cbind(X_hat_gdp_lag,esol) 
}  

r6<-rowbuilder(X_hat_gdp_lag,X_inform)$mean_row

# load the Regression model with gdp_i, gdp_j and lagged values and random effects----

X_hat_gdp_lag_raneff<-c()
for (tau in 1:length(t)){
  
  load(paste("3) Model/Regression with gdp and lagged values and random effects/results/res",tau,".RData",sep=""))
  X_hat_gdp_lag_raneff<-cbind(X_hat_gdp_lag_raneff,esol) 
}  

r7<-rowbuilder(X_hat_gdp_lag_raneff,X_inform)$mean_row

# Results - Tabele ----

rbind(r1,r2,r3,r4,r5,r6,r7)

