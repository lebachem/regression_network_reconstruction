rm(list=ls())

# set working directory----
setwd("/Workflow/")

sequence_new<-round(seq(-4,4,0.1),1)
length(sequence_new)
relative_median<-c()
relative_mean<-c()
coef<-c()

# Extract the data----

for (beta in sequence_new){
  load(paste("5) Simulations/results/beta=",beta,".RData"))

  relative_mean<-c(relative_mean,mean( (mse_ipfp2 )/mse_em2))
  relative_median<-c(relative_median,median( ( mse_ipfp2)/mse_em2))
  coef<-cbind(coef,beta_s)
  
}

# Plot the data----
pdf("5) Simulations/results/compare.pdf",height=5,width=10)
par (mfrow=c(1,1))
plot(relative_mean,type="l",xaxt="n",col="gray",lwd=2.5,ylim=c(0,12),ylab="RSS",xlab="beta",lty=1,cex.lab=1.2)
lines(relative_median,lwd=2.5,lty=1)
abline(h=1,lty=2)
abline(v=median(1:length(sequence_new)),lty=2)
axis(1,1:length(sequence_new),sequence_new)

dev.off()
