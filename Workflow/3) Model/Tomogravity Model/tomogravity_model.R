#-----------------------------------------------------------#
# File that implements the Tomogravity model                #
#-----------------------------------------------------------#

# remove the old stuff
rm(list=ls())
# Load the packges
library(networkTomography)

# set working directory----
setwd("/Workflow/")
load("1) Data Read In/data.RData")

# Create the routing matrix
rM<-routing_mat(n)

# Extract the matrix and the index set
A<-rM$A
iS<-rM$index_set

# find the optimal threshold value for the tomogravity model----

tomo_list<-list()
mse_tomo<-c()
grid<-seq(0,1,by=0.001)
up<-1
for (l in grid){
  print(l)
  fit<-t(tomogravity(Y,A,lambda=l)$Xhat)
  tomo_list[[up]]<-fit
  mse<-rowbuilder(fit,X_inform)$mean_L2
  mse_tomo<-c(mse_tomo,mse)
  print(mse)
up<-up+1
}

fit<-t(tomogravity(Y,A,lambda=0.011)$Xhat)

# save results----
rm(list=ls()[-which(ls()=="X_hat_tomogravity")])
save.image("3) Model/Tomogravity Model/results/Tomogravity.RData")


