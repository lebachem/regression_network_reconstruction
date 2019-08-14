#-----------------------------------------------------------#
# Simulation Study: RSS ratio between Regression and IPFP   #
#-----------------------------------------------------------#

rm(list=ls())

# Load the packges
library(ipfp)
library(nloptr)
library(pracma)


set.seed(1)

# set working directory----
setwd("/Workflow/")
load("1) Data Read In/data.RData")

# Define Parameters----
local_opts <- list( "algorithm" = "NLOPT_LD_AUGLAG",
                    "xtol_rel" = 1.0e-20 )
opts <- list( "algorithm" = "NLOPT_LD_AUGLAG",
              "xtol_rel" = 1.0e-20,
              "maxeval" = 1000,"local_opts" = local_opts )

n=10
replications<-1000
sequence<-round( seq(-4,4,by=0.1),1)
N<-n*(n-1)


for (beta in sequence){
  
  
  mse_em2<-c()
  mse_ipfp2<-c()
  beta_s<-c()
  for (rep in 1:replications){    
    
    
    X<-rnorm(N,0,1)
    
    A<-routing_mat(n)$A
    
    iS<-routing_mat(n)$index_set
    
    
    
    sender<-as.factor(iS[,1])
    receiver<-as.factor(iS[,2])
    
    
    
    u<-matrix(0,nrow=dim(iS)[1],ncol=n)
    v<-matrix(0,nrow=dim(iS)[1],ncol=n)
    Y_stack<-c()
    for (i in 1:dim(iS)[1]){
      sender<-iS[i,1]
      receiver<-iS[i,2]
      u[i,sender]<-1
      v[i,receiver]<-1
      Y_stack<-rbind(Y_stack,c(Y[sender],Y[21+receiver]))

    } 

    
    Z<-cbind(u,v,X)
    
    
    
    print(beta)
    
    # Draw random network----
    
    delta<-rnorm(n,0,1)
    gamma<-rnorm(n,0,1)
    comb<-rexp(N,rate=1/exp(Z%*%c(delta,gamma,beta)))
    Y=A%*%comb
    
    
    
    Mat<-matrix(0,ncol=n,nrow=n)
    
    
    for (i in 1:dim(iS)[1]){
      Mat[iS[i,1],iS[i,2]]<-comb[i]
    }
    
    # Define functions----
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
    
    

    
    X_hat<-ipfp(y=Y,A=A,x0=rep(1,n*(n-1)))
    
    vec_E<-X_hat
    
    # First iteration----
    nl<-nloptr(x0=rnorm(2*n+1,0,1), eval_jac_g_eq =  cgradient, eval_f = ll,eval_grad_f = llgradient,eval_g_eq =constraint,opts = opts)

    old<-nl$solution
    
    # Iterate until convergence----
    crit<-100
    while(crit>0.001){
      # E-step
      vec_E<-exp(Z%*%nl$solution)
      # M-Step
      nl<-nloptr(x0=nl$solution, eval_jac_g_eq =  cgradient, eval_f = ll,eval_grad_f = llgradient,eval_g_eq =constraint,opts = opts)
      
      crit<-t(nl$solution-old)%*%(nl$solution-old)
      #print(crit)
      old<-nl$solution
      
    }
    
    sol_em<-exp(Z%*%old)
    mse_em2<-c(mse_em2,mean( (sol_em-comb)^2  ))
    mse_ipfp2<-c(mse_ipfp2,mean( (X_hat-comb)^2  ))
    beta_s<-c(beta_s,old[length(old)])
  }
  
  # save results----
  save.image(paste("5) Simulations/results/beta=",beta,".RData"))
  
}


