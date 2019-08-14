#-----------------------------------------------------------#
# File that implements Ecological Regression                #
#-----------------------------------------------------------#

# remove the old stuff
rm(list=ls())

# Load the packges
library(glmnet)

# set working directory----
setwd("/Workflow/")
load("1) Data Read In/data.RData")

# Create the routing matrix
rM<-routing_mat(n)

# Extract the matrix and the index set
A<-rM$A

### Fit a EI Regression model in the direction Y_i+|Y_+i----
Y_i_plus=Y[,1:21]
Y_plus_i=Y[,22:42]
X_hat_EIreg<-c()

for (i in 1:n){
  coef_i<- lm(Y_i_plus[,i]~-1+Y_plus_i[,-i])$coef
  coef_i<-abs(min(coef_i))+coef_i
  coef_i<-coef_i/sum(coef_i)
  X_hat_EIreg<-rbind(X_hat_EIreg,t(Y_plus_i[,-i]*coef_i))
}
X_hat_EIreg1<-X_hat_EIreg

### Fit a EI Regression model in the direction Y_+i|Y_i+----
## Fit a model----
Y_i_plus=Y[,22:42]
Y_plus_i=Y[,1:21]
X_hat_EIreg<-c()

for (i in 1:n){
  coef_i<- lm(Y_i_plus[,i]~-1+Y_plus_i[,-i])$coef
  coef_i<-abs(min(coef_i))+coef_i
  coef_i<-coef_i/sum(coef_i)
  X_hat_EIreg<-rbind(X_hat_EIreg,t(Y_plus_i[,-i]*coef_i))
}
X_hat_EIreg2<-X_hat_EIreg



## Save the results----
rm(list=ls()[-c(which(ls()=="X_hat_EIreg1"),which(ls()=="X_hat_EIreg2"))])
save.image("3) Model/Ecological Regression/results/eireg.RData")
