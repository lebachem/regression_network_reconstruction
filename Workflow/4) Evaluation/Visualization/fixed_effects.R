#-----------------------------------------------------------#
# File that extracts the regression coefficients            #
#-----------------------------------------------------------#

# remove the old stuff
rm(list=ls())

# set working directory----
setwd("/Workflow/")

incl<-c(1:52)
delta<-c()
gamma<-c()
beta<-c()
predictions<-c()
for (tau in incl){
  load(paste("3) Model/Regression with gdp and trade/results/res",tau,".RData",sep=""))
  
  predictions<-cbind(predictions,esol)
  delta<-cbind(delta,parms[1:n,dim(parms)[2]])
  gamma<-cbind(gamma,parms[(n+1):(2*n),dim(parms)[2]])
  beta<-cbind(beta,parms[(2*n+1):(2*n+3),dim(parms)[2]])

  }


# Plot the fitted values----

pdf("4) Evaluation/Visualization/fitted.pdf",width = 10,height = 5)


plot(predictions[1,],type="l",ylim=c(min(predictions),max(predictions)),col="gray",xlim=c(-5,60),xaxt="n",xlab="t",main="Estimated bilateral Transfers (million USD)",ylab="",cex.lab=1.2)
for (i in 1:dim(predictions)[1]){
  lines(predictions[i,],type="l",col="lightgray")
}

axis(1,1:length(t),t,las=1)

rM<-rowMeans(predictions)
numb<-4
upper<-head(sort(rM,decreasing = T),numb)

for (i in 1:numb){
  lines(predictions[which(rM==upper[i]),],col=i,lwd=3,lty=i)
  text(y=predictions[which(rM==upper[i]),length(t)],x=58,combinations[which(rM==upper[i])] ,col=i)
  text(y=predictions[which(rM==upper[i]),1],x=-3,combinations[which(rM==upper[i])] ,col=i)
  
}
dev.off()

# Plot the coefficients----

pdf("4) Evaluation/Visualization/coefficients.pdf",height=4,width=8)
par(mfrow=c(1,3))
plot(beta[1,],type="l",lwd=2,main="GDP(i), sender",ylab="Beta",xlab="t",xaxt="n",ylim=c(0,0.75),cex.lab=1.2)
axis(1,1:length(t),t)
plot(beta[2,],type="l",lwd=2,main="GDP(j), receiver",ylab="Beta",xlab="t",xaxt="n",ylim=c(0,0.75),cex.lab=1.2)
axis(1,1:length(t),t)
plot(beta[3,],type="l",lwd=2,main="Trade(i,j), sender-receiver",ylab="Beta",xlab="t",xaxt="n",ylim=c(0,0.75),cex.lab=1.2)
axis(1,1:length(t),t)
dev.off()

pdf("4) Evaluation/Visualization/intercepts.pdf",height=5,width=9)

par(mfrow=c(1,2))
plot(delta[1,],type="l",lwd=1,main="Row-specific Intercept",ylab="delta",xlab="t",xaxt="n",ylim=c(-2,1.5),xlim=c(-3,60),cex.lab=1.2)
for (i in 1:n){
  lines(delta[i,],type="l",lwd=2,lty=i)
}
for (i in 1:n){
  lines(delta[i,],type="l",lwd=2,lty=i)
  text(y=delta[i,length(t)],x=58,countries[i],cex=1)
  text(y=delta[i,1],x=-3,countries[i],cex=1)
}
axis(1,1:length(t),t)



plot(gamma[1,],type="l",lwd=1,main="Column-specific Intercept",ylab="gamma",xlab="t",xaxt="n",ylim=c(-5,2),xlim=c(-3,60),cex.lab=1.2)
for (i in 1:n){
  lines(gamma[i,],type="l",lwd=2,lty=i)
  text(y=gamma[i,length(t)],x=58,countries[i],cex=1)
  text(y=gamma[i,1],x=-3,countries[i],cex=1)
}
axis(1,1:length(t),t)

dev.off()








