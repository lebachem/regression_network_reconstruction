#-----------------------------------------------------------#
# File that creates some descriptives                       #
#-----------------------------------------------------------#

# remove the old stuff
rm(list=ls())

# set working directory----
setwd("/Workflow/")
load("1) Data Read In/data.RData")


# Density of the network:
d<-c()
total<-c()
for (tau in 1:length(t)){
  d_t<-sum(as.numeric(X_inform[,tau]>0))
  total<-c(total,sum(X_inform[,tau]))
  d<-c(d,d_t)
}

# Number of zero-valued marginals

m<-c()

for (tau in 1:length(t)){
m<-c(m,sum(Y[tau,]==0))
}



# plot the network density and the aggregated development
## Plot the density and aggregated quantities----
pdf("2) Descriptives/density_and_total.pdf",width = 9,height = 4.5)
par(mfrow=c(1,3))
plot(d/(n*(n-1)),type="l",xaxt="n",lwd=2,main="Density of the network",ylab="",xlab="t",cex.lab=1.5,cex.axis=1.2)
axis(1,1:length(t),t,las=1,cex.axis=1.2)

plot(m/42,type="l",xaxt="n",lwd=2,main="Share of zero-valued marginals",ylab="",xlab="t",cex.lab=1.5,cex.axis=1.2)
axis(1,1:length(t),t,las=1,cex.axis=1.2)

total<-rowSums(Y[,1:n])/2
plot(total,type = "l",xaxt="n",main="Development of total exposures",ylab = "",ylim=c(0,max(total)*1.2),lwd=2,las=2,cex.axis=1.2,cex.lab=1.5,xlab="t")
axis(1,1:length(t),t,las=1,cex.axis=1.2)
dev.off()

# calculate the correlations between the variables

## Correlations between the variables----

X_long<-c()
Dist_long<-c()
GDP_longi<-c()
GDP_longj<-c()
Trade_long<-c()
X_t_1_long<-c()

for (tau in 1:52){
  up<-1
  for (i in 1:n){
    for (j in 1:n){
      if (i!=j){
        X_long<-c(X_long,X_inform[up,tau])
        Dist_long<-c(Dist_long,exp(log_Dist[i,j]))
        GDP_longi<-c(GDP_longi,exp(log_GDP[i,tau]))
        GDP_longj<-c(GDP_longj,exp(log_GDP[j,tau]))
        Trade_long<-c(Trade_long,exp(log_Trade[[tau]][i,j]))
        if (tau>1){
        X_t_1_long<-c(X_t_1_long,X_inform[up,tau-1])
          
        }
        up<-up+1
      }
    }
  }
  

}

# Show the correlation between X_ij and the covariates
cor(X_long,GDP_longi)
cor(X_long,GDP_longj)
cor(X_long,Trade_long)
cor(X_long,Dist_long)
cor(X_long[-c(1:(n*(n-1)))],X_t_1_long)

# remove everything
rm(list=ls())
