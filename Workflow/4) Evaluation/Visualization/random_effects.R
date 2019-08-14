#---------------------------------------------------------------#
# File that extracts the random effects from the fitted models  #
#---------------------------------------------------------------#

# remove the old stuff
rm(list=ls())

# set working directory----
setwd("/Workflow/")

incl<-c(1:52)
delta<-c()
gamma<-c()
sigma2<-c()


for (tau in incl){
  load(paste("3) Model/Regression with random effects/results/res",tau,".RData",sep=""))

  delta<-cbind(delta,parms[1:n,dim(parms)[2]])
  gamma<-cbind(gamma,parms[(n+1):(2*n),dim(parms)[2]])
  sigma2<-cbind(sigma2,c(exp(parms[(2*n+2),dim(parms)[2]]),exp(parms[(2*n+3),dim(parms)[2]]),parms[(2*n+4),dim(parms)[2]]))
  
  }


pdf("4) Evaluation/Visualization/sigmas_raneff.pdf",height=3,width=8)
# The Variance components----
par (mfrow=c(1,3))
plot(sigma2[1,],xaxt="n",type="l",lwd=2,xlab="t",ylab="var(delta)",ylim =c(1,13),cex.lab=1.2)
axis(1,1:length(t),t)

plot(sigma2[2,],xaxt="n",type="l",lwd=2,xlab="t",ylab="var(gamma)",ylim =c(1,13),cex.lab=1.2)
axis(1,1:length(t),t)

plot(sigma2[3,]/(sqrt(sigma2[1,])*sqrt(sigma2[2,])),xaxt="n",type="l",lwd=2,xlab="t",ylab="correlation(gamma,delta)",cex.lab=1.2)
axis(1,1:length(t),t)

dev.off()



