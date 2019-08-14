#-----------------------------------------------------------#
# File that implements the non-negative LASSO               #
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



### Search for optimal penalty term  ----
grid<-seq(0,100000,by=0.1)
lasso<-c()
up<-1
for (l in grid){
  X_hat_lasso<-c()
  for (tau in 1:length(t)){
    
    X_hat_lasso<-cbind(X_hat_lasso,as.matrix(glmnet(x=A,y=Y[tau,],lower.limits=0,lambda=l,standardize = FALSE)$beta))
    # print(tau)
  }
  mse<-rowbuilder(X_hat_lasso,X_inform)$mean_L2
  lasso<-c(lasso,mse)
  print(mse)
  up<-up+1
  print(l)
}

### Plot results of hyperparameter optimization ----
plot(lasso,xaxt="n")
axis(1,1:length(lasso),grid)
abline(v=which(lasso==min(lasso)))


### Fit the model with optimal hyperparameter ----
X_hat_lasso<-c()
for (tau in 1:length(t)){
  
  X_hat_lasso<-cbind(X_hat_lasso,as.matrix(glmnet(x=A,y=Y[tau,],lower.limits=0,lambda=grid[which(lasso==min(lasso))],standardize = FALSE)$beta))
  # print(tau)
}


### Save the results ----
rm(list=ls()[-which(ls()=="X_hat_lasso")])
save.image("3) Model/LASSO Model/results/lasso.RData")