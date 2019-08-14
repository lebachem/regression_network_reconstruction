#-----------------------------------------------------------#
# File that creates tables for comparing L_1 and L_2 errors #
#-----------------------------------------------------------#

# remove the old stuff
rm(list=ls())

library(stargazer)

# set working directory----
setwd("/Workflow/")
load("1) Data Read In/data.RData")

# load the Regression model with gdp_i, gdp_j, trade_ij----
X_hat_gdp_trade<-c()
for (tau in 1:length(t)){
  
  load(paste("3) Model/Regression with gdp and trade/results/res",tau,".RData",sep=""))
  X_hat_gdp_trade<-cbind(X_hat_gdp_trade,esol) 
}  

r1<-rowbuilder(X_hat_gdp_trade,X_inform)$mean_row

# load the Gravity model----
load("3) Model/Gravity Model/results/gravity.RData")

r2<-rowbuilder(X_hat_gravity,X_inform)$mean_row

# load the Tomogravity model----
load("3) Model/Tomogravity Model/results/tomogravity.RData")

r3<-rowbuilder(X_hat_tomogravity,X_inform)$mean_row

# load the LASSO model----
load("3) Model/LASSO Model/results/lasso.RData")

r4<-rowbuilder(X_hat_lasso,X_inform)$mean_row

# load the Ecological Regression models----
load("3) Model/Ecological Regression/results/eireg.RData")

r5<-rowbuilder(X_hat_EIreg1,X_inform)$mean_row
r6<-rowbuilder(X_hat_EIreg2,X_inform)$mean_row

# load the hierarchical models models----
load("3) Model/Hierarchical Models/results/hierarch.RData")

r7<-rowbuilder(meanX_hat_ER,X_inform)$mean_row
r8<-rowbuilder(meanX_hat_fitness,X_inform)$mean_row

# Results - Tabele ----

rbind(r1,r2,r3,r4,r5,r6,r7,r8)
