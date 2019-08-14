#--------------------------------------------------------------#
# File that implements the Hiearchical Models (EmpiricalBayes) #
#--------------------------------------------------------------#

# remove the old stuff
rm(list=ls())

# Load the packges
library(ipfp)
library(rlist)
library(nloptr)
library(pracma)
library(truncdist)
library(systemicrisk)


set.seed(123)

# set working directory----
setwd("/Workflow/")
load("1) Data Read In/data.RData")
N<-n*(n-1)

rep<-100

# Run 100 sequences of the ER and fitness model
X_hat_fitness<-list()
X_hat_ER<-list()
for (tau in 1:length(t)){
  print(tau)
  
  l<-Y[tau,1:n]
  l_index<-which(l==0)
  
  a<-Y[tau,(n+1):(2*n)]
  a_index<-which(a==0)
  
  dens<-sum(as.numeric(X_inform[,tau]>0))/N
  
  
  
  L_fixed<-matrix(NA,ncol=n,nrow=n)
  
  if (length(a_index)>0&length(l_index)>0){
    L_fixed[l_index,a_index]<-0
  }
  
  if (length(a_index)>0&length(l_index)==0){
    L_fixed[,a_index]<-0
  }
  
  if (length(a_index)==0&length(l_index)>0){
    L_fixed[l_index,]<-0
  }
  diag(L_fixed)<-0
  
  ## Fitness Model----
  model1<-calibrate_FitnessEmp(l,a,L_fixed = L_fixed,targetdensity = dens)
  Lsamp1 <- sample_HierarchicalModel(l=l,a=a,L_fixed = L_fixed,
                                     model=model1,nsamples=rep,thin=1e2)
  
  ## ER Model----
  model2<-calibrate_ER(l,a,L_fixed = L_fixed,targetdensity = dens)
  Lsamp2 <- sample_HierarchicalModel(l=l,a=a,L_fixed = L_fixed,
                                     model=model2,nsamples=rep,thin=1e2)
  
  
  vec1<-c()
  vec2<-c()
  for (r in 1:length(Lsamp1$L)){
    
    inf1<-c()
    inf2<-c()
    for (i in 1:n){
      for (j in 1:n){
        if (i !=j){
          inf1<-c(inf1,Lsamp1$L[[r]][i,j])
          inf2<-c(inf2,Lsamp2$L[[r]][i,j])
        }
      }
    }
    

    vec1<-cbind(vec1,inf1)
    vec2<-cbind(vec2,inf2)
    
    
  }
  
  
  X_hat_fitness[[tau]]<-vec1
  X_hat_ER[[tau]]<-vec2
  
}

# Evaluate E[X_hat]----

meanX_hat_fitness<-c()
meanX_hat_ER<-c()

for (tau in 1:length(t)){
  mt1<-c()
  mt2<-c()
  for (r in 1:rep){
    mt1<-cbind(mt1,X_hat_fitness[[tau]][,r])
    mt2<-cbind(mt2,X_hat_ER[[tau]][,r])
    
  }
  meanX_hat_fitness<-cbind(meanX_hat_fitness,rowMeans(mt1))  
  meanX_hat_ER<-cbind(meanX_hat_ER,rowMeans(mt2))  
  
}



## save the results----
rm(list=ls()[-c(which(ls()=="meanX_hat_ER"),which(ls()=="meanX_hat_fitness"))])
save.image("3) Model/Hierarchical Models/results/hierarch.RData")


