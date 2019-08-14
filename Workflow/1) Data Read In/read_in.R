#-----------------------------------------------------------#
# File that reads the data into R                           #
#-----------------------------------------------------------#

# remove the old stuff
rm(list=ls())

# load packages
library(stargazer)
library(rlist)

# set working directory----
setwd("/Workflow/")

#-----------------------------------------------------------#
#  load the BIS data (available form www.bis.org)           #
#-----------------------------------------------------------#
## BIS data----
Data<-read.csv("Data/BIS_bilateral_exposures.csv",sep=",")

# identify the different time-points
t<-unique(Data[,1])

# identify the countries included
countries<-unique(Data[,3])

# define n as the number of countries
n<-length(countries)

# different country-combinations
combinations<-c()
for (i in 1:n){
  for (j in 1:n){
    if (i!=j){
      combinations<-c(combinations,paste(countries[i],"-",countries[j],sep=""))
    }
  } 
}

#-----------------------------------------------------------#
# split the data inro informative data X_inform that        #
# contains all dyads and a traning set Y that only contains # 
# the row and column sums                                   #
#-----------------------------------------------------------#

# open a container for the informative data
X_inform<-matrix(0,nrow=n*(n-1),ncol=length(t))

# give names to rows and colums
colnames(X_inform)<-t
rownames(X_inform)<-combinations

# Include the dyadic values in the informative set
for (tau in 1:length(t)){
  for (i in 1:length(combinations)){
    if (length(Data[Data$net_id==t[tau]&Data$arc_id==combinations[i],5])>0){
    X_inform[i,tau]<- Data[Data$net_id==t[tau]&Data$arc_id==combinations[i],5]
    }
    else{X_inform[i,tau]<-0}
  }
}



# create a (Tx2n) matrix containing n=21 claims and n=21 liabilities in each row (standing for one year)
Y<-matrix(0,nrow=length(t),ncol=2*n)
row.names(Y)<-t
colnames(Y)<-c(paste("from_",countries,sep=""),paste("to_",countries,sep=""))


# Now we fill the matrix with "from" and "to"
for (tau in 1:length(t)){
  for (i in 1:length(countries)){
    if (length(Data[Data$net_id==t[tau]&Data$from_id==countries[i],5])>0){
    Y[tau,i]<-sum(Data[Data$net_id==t[tau]&Data$from_id==countries[i],5])
    }
    else{Y[tau,i]=0}
    
    if (length(Data[Data$net_id==t[tau]&Data$to_id==countries[i],5])>0){
    Y[tau,i+n]<-sum(Data[Data$net_id==t[tau]&Data$to_id==countries[i],5])
    }
    else{Y[tau,i+n]<-0}
  }
}

#-----------------------------------------------------------#
# define a function that constructs the rounting matrix A   #
# (AX=Y) and a corresponding index set, based on the number #
# of nodes n.                                               #
#-----------------------------------------------------------#
## Routing matrix----
routing_mat<-function(n_nodes){
  
  n_edges<-n_nodes*(n_nodes-1)
  n_flows<-2*n_nodes
  A<-matrix(0,ncol=n_edges,nrow=n_flows)
  
  index_set<-cbind(rep(1:n_nodes,each=n_nodes),rep(1:n_nodes,n_nodes))
  index_set<-index_set[-which(index_set[,1]==index_set[,2]),]
  
  for (i in 1:n_nodes){
    for (j in 1:n_edges){
      if(index_set[j,1]==i){
        A[i,j]<-1
      }
    }
  }
  
  for (i in (n_nodes+1):n_flows ){
    for (j in 1:n_edges){
      if(index_set[j,2]==(i-n_nodes)){
        A[i,j]<-1
      }
    }
  }
  return(list(A=A,index_set=index_set))
}

#-----------------------------------------------------------#
# Remove all the stuff that is not needed anymore           #
#-----------------------------------------------------------#

rm(list=ls()[-which(ls()%in%c("Y","countries","t","combinations","routing_mat","n","X_inform","X_fna_entrop","X_fna_min_dens1","X_fna_min_dens2","X_fna_min_dens3","X_fna_min_dens4","X_fna_min_dens5","X_","d","logY"))])

#-----------------------------------------------------------#
# Now we read in the exogenuous information                 #
#-----------------------------------------------------------#

# Define different conversions for country identifier
country_convert_BIS<-c("AT","AU","BE","CA","CH","CL","DE","ES","FI","FR","GB","GR","IE","IT","JP","KR","NL","SE","TR","TW","US")
country_convert_IMF<-c("AUT","AUS","BEL","CAN","CHE","CHL","DEU","ESP","FIN","FRA","GBR","GRC","IRL","ITA","JPN","KOR","NLD","SWE","TUR","TWN","USA")
country_convert_COW<-c("Austria","Australia","Belgium","Canada","Switzerland","Chile","Germany","Spain","Finland","France","United Kingdom","Greece","Ireland","Italy","Japan","South Korea","Netherlands","Sweden","Turkey","Taiwan","United States of America")

count<-cbind(cbind(country_convert_BIS[1:7],country_convert_COW[1:7]),
             cbind(country_convert_BIS[8:14],country_convert_COW[8:14]),
             cbind(country_convert_BIS[15:21],country_convert_COW[15:21])
)

# Table with country codes and full country names from the publication
stargazer(count,type="text")

# Define the years to be included in the analysis
time_frame<-2005:2017

# We import the IMF GDP data
# Source: https://www.imf.org/en/Data
## GDP data----
GDP_data<-read.csv("Data/gdp.csv")

# Save the logarithmic GDP Data quarterly as node-Specific Variable
log_GDP<-c()
for (i in 1:13){
  log_GDP<-cbind(log_GDP,cbind(log(GDP_data[,i]),log(GDP_data[,i]),log(GDP_data[,i]),log(GDP_data[,i])))
}


# We import the COW dyadic trade data
# Source: http://www.correlatesofwar.org/data-sets/bilateral-trade
## Trade data----
trade<-read.csv("Data/Dyadic_COW_4.0.csv")

# Put Everything in a list
trade_list<-list()

for (tau in 2005:2014){
  trade_t<-trade[trade[,3]==tau,]
  Trade_mat<-matrix(0,ncol=length(country_convert_COW),nrow=length(country_convert_COW))
  for (i in 1:length(country_convert_COW)){
    for (j in 1:length(country_convert_COW)){
      if (i!=j){
        
        receiver<-country_convert_COW[i]
        sender<-country_convert_COW[j]
        if (length(trade_t[trade_t[,4]==sender&trade_t[,5]==receiver,6])>0){
          Trade_mat[i,j]<-trade_t[trade_t[,4]==sender&trade_t[,5]==receiver,6]
        }
        if (length(trade_t[trade_t[,4]==sender&trade_t[,5]==receiver,6])==0){
          Trade_mat[i,j]<-trade_t[trade_t[,4]==receiver&trade_t[,5]==sender,7]
        }
      }
    } 
  }
  trade_list[[tau-2004]]<-Trade_mat
}


# Extend the data to the periods 2015,2016,2017

trade_list[[11]]<-matrix(0,ncol=length(country_convert_COW),nrow=length(country_convert_COW))
trade_list[[12]]<-matrix(0,ncol=length(country_convert_COW),nrow=length(country_convert_COW))
trade_list[[13]]<-matrix(0,ncol=length(country_convert_COW),nrow=length(country_convert_COW))

# Use an autoregressive Model to impute future values of the Trade
for (i in 1:21){
  for (j in 1:21){
    if (i!=j){
      
      trade_series<-c()
      for (tau in 1:10){
        trade_series<-c(trade_series,trade_list[[tau]][i,j])
      }
      X<-cbind(1,trade_series[-length(trade_series)])
      b<-solve(t(X)%*%X)%*%t(X)%*%trade_series[-1]
      hat2015<-b[1]+b[2]*trade_series[length(trade_series)]
      trade_list[[11]][i,j]<-hat2015
      hat2016<-b[1]+b[2]*hat2015
      trade_list[[12]][i,j]<-hat2016
      hat2017<-b[1]+b[2]*hat2016
      trade_list[[13]][i,j]<-hat2017
    }
    
  }
}



log_Trade<-list()
up<-1
for (i in 1:13){

    log_Trade<- list.append(log_Trade ,log(1+trade_list[[i]]))
    log_Trade<- list.append(log_Trade ,log(1+trade_list[[i]]))
    log_Trade<- list.append(log_Trade ,log(1+trade_list[[i]]))
    log_Trade<- list.append(log_Trade ,log(1+trade_list[[i]]))
}


country_convert_COW<-c("AUS","AUL","BEL","CAN","SWZ","CHL","GMY","SPN","FIN","FRN","UK","GRC","IRE","ITA","JPN","ROK","NTH","SWD","TUR","TAW","USA")

# We import the distance data
# Source: http://ksgleditsch.com/data-3.html
## Distance data----
dist<-read.csv("Data/capdist.csv")
dist<-data.frame(dist)
dist[,2]<-as.character(dist[,2])
dist[,4]<-as.character(dist[,4])

D<-matrix(0,ncol=21,nrow=21)

for (i in country_convert_COW){
  for (j in country_convert_COW){
    if (i !=j){
      if (length(which(dist[,2]==i&dist[,4]==j)>0)){
        a<-dist[which(dist[,2]==i&dist[,4]==j),5]
      } else {
        a<-dist[which(dist[,4]==i&dist[,2]==j),5]}
      D[which(country_convert_COW==i),which(country_convert_COW==j)]<-a
    }
  } 
}

log_Dist<-log(D)

## Function to evaluate the L_1 and L_2 Error----
rowbuilder<-function(X_hat,X_inform,divide=100000){
  len=dim(X_hat)[2]
  L1<-c()
  L2<-c()
  
  for (tau in 1:len){
    L2<-c(L2,sqrt( sum(  (c(X_hat[,tau]/divide-X_inform[,tau]/divide))^2 )  ))
    L1<-c(L1, sum(  (abs(c(X_hat[,tau]/divide-X_inform[,tau]/divide)) )  ))
  }
  
  overall_L2=sqrt( sum(  (c( c(X_hat/divide)-c(X_inform/divide)))^2 )  )
  overall_L1=sum(  (abs(c(c(X_hat/divide)-c(X_inform/divide))) )  )
  mean_row=c(overall_L1,overall_L2,mean(L1),sd(L1),mean(L2),sd(L2))
  return(list(mean_row=mean_row, mean_L1=mean(L1),mean_L2=mean(L2), sd_L1=sd(L1),sd_L2=sd(L2),overall_L1=overall_L1,overall_L2=overall_L2))
  
}

rm(list=ls()[-which(ls()%in%c("Y","countries","t","combinations","routing_mat","n","X_inform","log_Dist","log_GDP","log_Trade","rowbuilder"))])

## Save the data----
save.image("1) Data Read In/data.RData")

# remove everything
rm(list=ls())
