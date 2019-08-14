#-----------------------------------------------------------#
# File that implements the Gravity Model                    #
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

## Fit the model----
# Calculate the 
X_hat_gravity<-t(gravity(Y,getSrcDstIndices(A)))

## Save the results----
rm(list=ls()[-which(ls()=="X_hat_gravity")])
save.image("3) Model/Gravity Model/results/gravity.RData")