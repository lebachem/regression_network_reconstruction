#---------------------------------------------------------------#
# File that extracts predicton intervals from the fitted models #
#---------------------------------------------------------------#

# remove the old stuff
rm(list=ls())
# set working directory----
setwd("/Workflow/")

cover<-c()
for (tau in c(1:52)){
  load(paste("3) Model/Regression with gdp and trade/results/res",tau,".RData",sep=""))

comb<-X_inform[,tau]

Y_long<-routing_mat(n)$A%*%comb
if (length(row)>0){
  comb<-comb[-c(which(iS_old[,1]%in%row))]
}
if (length(column)>0){
  comb<-comb[-c(which(iS_old[,2]%in%column))]
}


CI_coverage<-c()
for (i in 1:dim(iS)[1]){
  q_low<- quantile(PI[i,],0.005)
  if (q_low<0){q_low=0}
  q_high<- quantile(PI[i,],0.955)
  if (q_high>min(Y_stack[i,])){q_high<-min(Y_stack[i,])}                             
  CI_coverage<-rbind(CI_coverage,c(comb[i],q_low,q_high))
}
false<-which(data.table::between(CI_coverage[,1],CI_coverage[,2],CI_coverage[,3])==F)
col=rep(1,dim(iS)[1])
col[false]<-2


cover<-c(cover,sum(data.table::between(CI_coverage[,1],CI_coverage[,2],CI_coverage[,3])==T)/dim(iS)[1])
}

# Plot the prediction intervals----

pdf("4) Evaluation/Visualization/coverage.pdf",height=6,width=11)
par(mfrow=c(1,2))
plot(sol_init,CI_coverage[,1],ylim=c(0,max(CI_coverage[,1])),xlim=c(0,max(CI_coverage[,3])),xlab="Estimates + 95%-Prediction Interval in Gray",ylab="Real Values",lwd=1,col=col,main=paste("Share covered in",t[tau],"=",round(print(sum(data.table::between(CI_coverage[,1],CI_coverage[,2],CI_coverage[,3])==T)/dim(iS)[1] ),2),sep=" "),cex.lab=1.3)
abline(a=0,b=1)
segments(x0=CI_coverage[,2], y0=CI_coverage[,1], x1 = CI_coverage[,3], y1 = CI_coverage[,1],col="gray")
points(sol_init[false],CI_coverage[false,1],col=2,lwd=2)
high<-sort(CI_coverage[,1],decreasing = T)[1:6]
text(sol_init[which(CI_coverage[,1]%in%high)],CI_coverage[which(CI_coverage[,1]%in%high),1],paste(countries[iS[which(CI_coverage[,1]%in%high),1]],"-",countries[iS[which(CI_coverage[,1]%in%high),2]],sep=""))

plot(cover,type="l",lwd=2,main=paste("Real values within prediction intervals. Mean=",round(mean(cover),2),sep = " "),ylab="Share",xlab="t",xaxt="n",ylim=c(0.90,1),cex.lab=1.3)
abline(h=0.95,lty=2)
axis(1,1:length(t),t)

dev.off()
