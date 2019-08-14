#-----------------------------------------------------------#
# File that implements the Regression-based approach        #
#-----------------------------------------------------------#

# remove the old stuff
rm(list=ls())

# Load the packges
library(ipfp)
library(nloptr)
library(pracma)

set.seed(1)

# set working directory----
setwd("/Workflow/")
load("1) Data Read In/data.RData")



# Define the paramters for the non-linear optimization----
local_opts <- list( "algorithm" = "NLOPT_LD_AUGLAG",
                    "xtol_rel" = 1.0e-20 )
opts <- list( "algorithm" = "NLOPT_LD_AUGLAG",
              "xtol_rel" = 1.0e-20,
              "maxeval" = 1000,"local_opts" = local_opts )

EM_solution<-c()
parms<-c()
mse_old<-c()

for (tau in 1:length(t)){

  
  comb<-X_inform[,tau]
  
  A<-routing_mat(n)$A
  Y=A%*%comb
  
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
  
  
  
  
  N<-dim(iS)[1]
  n_sender<-n-length(row)
  n_receiver<-n-length(column)
  u<-matrix(0,nrow=dim(iS)[1],ncol=n)
  v<-matrix(0,nrow=dim(iS)[1],ncol=n)
  X<-matrix(0,nrow=dim(iS)[1],ncol=1)
  Y_stack<-c()
  for (i in 1:dim(iS)[1]){
    sender<-iS[i,1]
    receiver<-iS[i,2]
    u[i,sender]<-1
    v[i,receiver]<-1
    X[i,1]<-1

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
  U<-cbind(u,v)
  
  
  G_full<-matrix(0,ncol=2*n,nrow=2*n)
  
  
  for (p in 1:dim(iS_old)[1]){
    G_full[iS_old[p,1],n+iS_old[p,2]]<-1
  }
  
  G_reduced<-G_full+t(G_full)  
  G_reduced[G_reduced!=0]<-1
  
  if (length(row)>0){
    G_reduced<-G_reduced[-row,-row]
  }
  if (length(column)>0){
    G_reduced<-G_reduced[-(column+n),-(column+n)]
  }
  
  
  
  
  X_hat<-ipfp(y=Y,A=A,x0=rep(1,dim(iS)[1]))
  vec_E<-X_hat
  
  ## Now we prepare the EM-Algorithm----
  
  ll<-function(theta){
    s<-theta[1:n_sender]
    r<-theta[(n_sender+1):(n_sender+n_receiver)]
    theta_tilde<-theta[(n_sender+n_receiver+1)]
    
    sigma_s<-exp(sol_2[1])
    sigma_r<-exp(sol_2[2])
    sigma_rs<-sol_2[3]
    G<-diag(c(rep(sigma_s,n_sender),rep(sigma_r,n_receiver)))+G_reduced*sigma_rs
    
    logl<-sum(-Z%*%c(s,r,theta_tilde)-vec_E/exp(Z%*%c(s,r,theta_tilde)))-t(c(s,r))%*%solve(G)%*%c(s,r)/2
    
    return(-logl)
  }
  
  llgradient<-function(theta){
    grad(f=ll,x0=theta)
    
  }
  
 

  constraint<-function(theta){
    
    return(A%*%exp(Z%*%theta)-Y)
    
  }
  
  cgradient<-function(theta){
    A%*%diag(c(exp(Z%*%theta)))%*%Z
  }
  
  
  # Solve for beta and the random effects given G----
  sol_2<-c(0,0,0)
  nl<-nloptr(x0=rep(1,n_sender+n_receiver+1), eval_jac_g_eq =  cgradient, eval_f = ll,eval_grad_f = llgradient,eval_g_eq =constraint,opts = opts)
  sol_1<-nl$solution
  
  # Define the restricted likelihood----
  l_R<-function(sigma){
    sigma_s<-exp(sigma[1])
    sigma_r<-exp(sigma[2])
    sigma_rs<-sigma[3]
    
    Sigma=diag( c(  rep(sigma_s,n_sender),rep(sigma_r,n_receiver) ) )+G_reduced*sigma_rs
    V=diag(rep(1,N))+U%*%Sigma%*%t(U)
    Y_tilde=Z%*%sol_1+diag(c(exp(-Z%*%sol_1)))%*%(vec_E-exp(Z%*%sol_1)  )
    re=-log(det(V))-log( det( t(X)%*%solve(V)%*%X )  ) - t(Y_tilde-X%*%sol_1[(n_sender+n_receiver+1):(length(sol_1))])%*%solve(V)%*%(Y_tilde-X%*%sol_1[(n_sender+n_receiver+1):(length(sol_1))])
    
        return(-re)  
  }

  
  # Solve for G(vartheta) given beta and the random effects----
  sol_2<-optim(par=sol_2,fn=l_R,method = "Nelder-Mead")$par
  vec_E<-exp(Z%*%sol_1)

  
  old<-nl$solution
  
  crit<-100
  while(crit>0){
    nl<-nloptr(x0=sol_1, eval_jac_g_eq =  cgradient, eval_f = ll,eval_grad_f = llgradient,eval_g_eq =constraint,opts = opts)
    sol_1<-nl$solution
    sol_2<-optim(par=sol_2,fn=l_R,method="Nelder-Mead")$par
    vec_E<-exp(Z%*%sol_1)

    crit<-t(nl$solution-old)%*%(nl$solution-old)
    old<-nl$solution
    
  }
  
  
  esol<-exp(Z%*%nl$solution)
  
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
  
  
  Y=routing_mat(n)$A%*%comb
  Y_stack<-c()
  for (i in 1:dim(iS_old)[1]){
    sender<-iS_old[i,1]
    receiver<-iS_old[i,2]
    
    Y_stack<-rbind(Y_stack,c(Y[sender],Y[21+receiver]))
   
  } 
  

  EM_solution<-cbind(EM_solution,esol)
  
  parameter_na<-rep(0,2*n+1)
  parameter_na[row]<-NA
  parameter_na[column+n]<-NA
  
  up<-1
  for (i in 1:length(parameter_na)){
    if (is.na(parameter_na[i])==F){
      parameter_na[i]<-nl$solution[up]
      up<-up+1
    }
  }
  
  
  parms<-cbind(parms,c(parameter_na,sol_2))
  
  # save results----
  
  save.image(paste("3) Model/Regression with random effects/results/res",tau,".RData",sep=""))
  
  mse_old<-c(mse_old,mean((comb-esol)^2))
  print(tau)
}

