#-----------------------------------------------------------#
# File that implements the IPFP Model                       #
#-----------------------------------------------------------#

# remove the old stuff
rm(list=ls())

# Load the packges
library(ipfp)

# set working directory----
setwd("/Workflow/")
load("1) Data Read In/data.RData")

# Create the routing matrix
rM<-routing_mat(n)

# Extract the matrix and the index set
A<-rM$A

## Fit a model----
X_hat_entrop<-c()
for (tau in 1:length(t)){
  X_hat_entrop<-rbind(X_hat_entrop,ipfp(y=Y[tau,],A=A,x0=rep(1,n*(n-1))))
  print(tau)
}

X_hat_ipfp<-t(X_hat_entrop)

## Save the results----
rm(list=ls()[-which(ls()=="X_hat_ipfp")])
save.image("3) Model/IPFP Model/results/ipfp.RData")