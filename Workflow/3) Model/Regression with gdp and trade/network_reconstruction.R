#-----------------------------------------------------------#
# File that implements the Regression-based approach        #
#-----------------------------------------------------------#

# remove the old stuff
rm(list=ls())

# Load the packges
library(ipfp)
library(nloptr)
library(pracma)


# set working directory----
setwd("/Workflow/")
load("1) Data Read In/data.RData")

set.seed(1)

# Define the functions that are needed----

ll<-function(theta){
  logl<-sum(-Z%*%theta-vec_E/exp(Z%*%theta))
  return(-logl)
}

llgradient<-function(theta){
  emp_mat<-matrix(0,nrow=length(vec_E),ncol=length(vec_E))  
  diag(emp_mat)<-vec_E/exp(Z%*%theta)-1
  ob<-emp_mat%*%Z
  return(colSums(ob))
}

constraint<-function(theta){
  
  return(A%*%exp(Z%*%theta)-Y)
  
}

cgradient<-function(theta){
  A%*%diag(c(exp(Z%*%theta)))%*%Z
}

# Define the paramters for the non-linear optimization----
local_opts <- list( "algorithm" = "NLOPT_LD_AUGLAG",
                    "xtol_rel" = 1.0e-20 )
opts <- list( "algorithm" = "NLOPT_LD_AUGLAG",
              "xtol_rel" = 1.0e-20,
              "maxeval" = 1000,"local_opts" = local_opts )

EM_solution<-c()
parms<-c()
coverage<-c()
mse_old<-c()
coverage_sd_parameters<-c()

for (tau in 1:length(t)){

# Version with the Covariates gdp_i, gdp_j and trade_ij
Xi<-log_GDP[,tau]
Xj<-log_GDP[,tau]
Xij<-log_Trade[[tau]]
# Save the real values
comb<-X_inform[,tau]

# Create the routing matrix
A<-routing_mat(n)$A
# Calculate the marginals
Y=A%*%comb

# Remove the zero valued margins
killer<-which(Y==0)
row<-killer[which(killer<=n)]
column<-killer[which(killer>n)]-n

if (length(killer)>0){
  A<-A[-killer,]
}

iS<-routing_mat(n)$index_set
iS_old<-routing_mat(n)$index_set


if (length(row)>0&length(column)==0){
  A<-A[,-c(which(iS[,1]%in%row))]
}
if (length(column)>0&length(row)==0){
  A<-A[,-c(which(iS[,2]%in%column))]
}
if (length(column)>0&length(row)>0){
  A<-A[,-c(which(iS[,1]%in%row,which(iS[,2]%in%column)))]
}


if (length(row)>0&length(column)==0){
  iS<-iS[-c(which(iS_old[,1]%in%row)),]
}
if (length(column)>0&length(row)==0){
  iS<-iS[-c(which(iS_old[,2]%in%column)),]
}
if (length(column)>0&length(row)>0){
  iS<-iS[-c(which(iS[,1]%in%row,which(iS[,2]%in%column))),]
}




n_sender<-n-length(row)
n_receiver<-n-length(column)
  u<-matrix(0,nrow=dim(iS)[1],ncol=n)
  v<-matrix(0,nrow=dim(iS)[1],ncol=n)
  X<-matrix(0,nrow=dim(iS)[1],ncol=3)

  Y_stack<-c()
  for (i in 1:dim(iS)[1]){
  sender<-iS[i,1]
  receiver<-iS[i,2]
  u[i,sender]<-1
  v[i,receiver]<-1
  X[i,1]<-Xi[sender]
  X[i,2]<-Xj[receiver]
  X[i,3]<-Xij[sender,receiver]
  Y_stack<-rbind(Y_stack,c(Y[sender],Y[21+receiver]))

  } 
 
if (length(row)>0){
  u<-u[,-row]
}
if (length(column)>0){
  v<-v[,-column]
}

if (length(killer)>0){
Y<-Y[-killer]
}
Z<-cbind(u,v,X)

# IPFP solution as starting values
X_hat<-ipfp(y=Y,A=A,x0=rep(1,dim(A)[2]))
vec_E<-X_hat

# First non-linear optimizatin step----
nl<-nloptr(x0=rep(0,n_sender+n_receiver+3), eval_jac_g_eq =  cgradient, eval_f = ll,eval_grad_f = llgradient,eval_g_eq =constraint,opts = opts)
old<-nl$solution

# Iterate until convergence----
crit<-100
while(crit>0){
  # E-step
  vec_E<-exp(Z%*%nl$solution)
  # M-Step
  nl<-nloptr(x0=nl$solution, eval_jac_g_eq =  cgradient, eval_f = ll,eval_grad_f = llgradient,eval_g_eq =constraint,opts = opts)

  crit<-t(nl$solution-old)%*%(nl$solution-old)
  print(crit)
  old<-nl$solution

}

# Save the initial predictions
sol_init<-exp(Z%*%nl$solution)
delta_est<-nl$solution[1:n_sender]
gamma_est<-nl$solution[(n_sender+1):(n_sender+n_receiver)]
beta_est<-nl$solution[(n_sender+n_receiver+1):length(nl$solution)]
Y_save<-Y
N<-dim(iS)[1]

# Draw 100 Bootstrap samples----
B<-100
theta_sim<-c()
sol_sim<-c()
PI<-c()


for (b in 1:B){
  
  comb<-rexp(N,rate=1/exp(Z%*%c(delta_est,gamma_est,beta_est)))
  Y=A%*%comb

  
  # IPFP solution as starting values
  X_hat<-ipfp(y=Y,A=A,x0=rep(1,dim(A)[2]))
  vec_E<-X_hat
  
  # First E-Step
  nl<-nloptr(x0=rnorm(n_sender+n_receiver+3,0,1), eval_jac_g_eq =  cgradient, eval_f = ll,eval_grad_f = llgradient,eval_g_eq =constraint,opts = opts)
  old<-nl$solution
  
  
  crit<-100
  while(crit>0){
    # E-Step
    vec_E<-exp(Z%*%nl$solution)
    # M-Step
    nl<-nloptr(x0=nl$solution, eval_jac_g_eq =  cgradient, eval_f = ll,eval_grad_f = llgradient,eval_g_eq =constraint,opts = opts)
    
    crit<-t(nl$solution-old)%*%(nl$solution-old)
    old<-nl$solution
    
  }
  
  # Save the paramters and predictions from the bootstrap
  theta_sim<-cbind(theta_sim,nl$solution)
  # Save the fitted values
  soli<-exp(Z%*%nl$solution)
  sol_sim<-cbind(sol_sim,soli)
  # We construct the predictin error as x*-\mu*
  PI<-cbind(PI, sol_init +comb-soli)
  
  Y_sim<-c()
  for (i in 1:dim(sol_sim)[2]){
    Y_sim<-cbind(Y_sim,A%*%sol_sim[,i])
  }
print(b)
}


plot(rowMeans(Y_sim)~Y_save)
abline(lm(rowMeans(Y_sim)~Y_save))
abline(a=0,b=1)
###
comb<-X_inform[,tau]

if (length(row)>0){
  comb<-comb[-c(which(iS_old[,1]%in%row))]
}
if (length(column)>0){
  comb<-comb[-c(which(iS_old[,2]%in%column))]
}

CI_coverage<-c()
for (i in 1:dim(iS)[1]){
  q_low<- quantile(PI[i,],0.025)

  q_high<- quantile(PI[i,],0.975)
  if (q_high>min(Y_stack[i,])){q_high<-min(Y_stack[i,])}      
  CI_coverage<-rbind(CI_coverage,c(comb[i],q_low,q_high))
}


esol<-sol_init

X_hat_imp<-c()
esol_imp<-c()
for (i in 1:dim(iS_old)[1]){
  sender<-iS_old[i,1]
  receiver<-iS_old[i,2]
  
  if (length(which(iS[,1]==sender&iS[,2]==receiver))>0){
    X_hat_imp<-c(X_hat_imp,X_hat[which(iS[,1]==sender&iS[,2]==receiver)])
    esol_imp<-c(esol_imp,esol[which(iS[,1]==sender&iS[,2]==receiver)])
  }
  if (length(which(iS[,1]==sender&iS[,2]==receiver))==0){
    X_hat_imp<-c(X_hat_imp,0)
    esol_imp<-c(esol_imp,0)
  }
  
  
}
esol<-esol_imp
X_hat<-X_hat_imp

init_par<-c(delta_est,gamma_est,beta_est)

sd_parameters<-c()
for (i in 1:(n_sender+n_receiver+3)){
  
  sd_parameters<-c(sd_parameters,sd(theta_sim[i,]))

  }


comb<-X_inform[,tau]

EM_solution<-cbind(EM_solution,esol)

parameter_na<-rep(0,2*n+3)
sd_parameter_na<-rep(0,2*n+3)

parameter_na[row]<-NA
parameter_na[column+n]<-NA
sd_parameter_na[row]<-NA
sd_parameter_na[column+n]<-NA

up<-1
for (i in 1:length(parameter_na)){
  if (is.na(parameter_na[i])==F){
    parameter_na[i]<-init_par[up]
    sd_parameter_na[i]<-sd_parameters[up]
    up<-up+1
  }
}


parms<-cbind(parms,parameter_na)
coverage_sd_parameters<-cbind(coverage_sd_parameters,sd_parameter_na)

coverage<-c(coverage,sum(data.table::between(CI_coverage[,1],CI_coverage[,2],CI_coverage[,3])==T)/(dim(iS)[1]) )
mse_old<-c(mse_old,mean((comb-esol)^2))
print(tau)
# Save the results----

save.image(paste("3) Model/Regression with gdp and trade/results/res",tau,".RData",sep=""))

}


