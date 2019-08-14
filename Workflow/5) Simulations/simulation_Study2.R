#-----------------------------------------------------------#
# Simulation Study: Bias in Regression and IPFP             #
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

# Define Parameters----

local_opts <- list( "algorithm" = "NLOPT_LD_AUGLAG",
                    "xtol_rel" = 1.0e-20 )
opts <- list( "algorithm" = "NLOPT_LD_AUGLAG",
              "xtol_rel" = 1.0e-20,
              "maxeval" = 1000,"local_opts" = local_opts )

n=10
replications<-1000
sequence<- c(-1,0,1)
N<-n*(n-1)


for (beta in sequence){
  
  
  mse_em2<-c()
  mse_ipfp2<-c()
  beta_s<-c()
  deviation<-c()
  deviation2<-c()
  ldeviation<-c()
  solution<-c()
  stat<-c()
  mean_dev<-c()
  
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
  
  
  # Simulate delta and gamma----
  print(beta)
  delta<-rnorm(n,0,1)
  gamma<-rnorm(n,0,1)  
  
  rM<-A%*%exp(Z%*%c(delta,gamma,beta)) 
  
  for (rep in 1:replications){    
    
    
    # draw a network----
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
    while(crit>0){
      # E-step
      vec_E<-exp(Z%*%nl$solution)
      # M-Step
      nl<-nloptr(x0=nl$solution, eval_jac_g_eq =  cgradient, eval_f = ll,eval_grad_f = llgradient,eval_g_eq =constraint,opts = opts)
      
      crit<-t(nl$solution-old)%*%(nl$solution-old)
      old<-nl$solution
      
    }
    
    # Save results----
    sol_em<-exp(Z%*%old)
    mse_em2<-c(mse_em2,mean( (sol_em-exp(Z%*%c(delta,gamma,beta)))^2  ))
    mse_ipfp2<-c(mse_ipfp2,mean( (X_hat-exp(Z%*%c(delta,gamma,beta)))^2  ))
    beta_s<-c(beta_s,old[length(old)])
    deviation<-cbind(deviation,sol_em-exp(Z%*%c(delta,gamma,beta)))
    deviation2<-cbind(deviation2,X_hat-exp(Z%*%c(delta,gamma,beta)))
    ldeviation<-cbind(ldeviation,log(sol_em)-Z%*%c(delta,gamma,beta))
    solution<-cbind(solution,sol_em)

    
    
  }
  
  # Plot the results----
  if (beta==0){
    low<- -4
    high<- 4
  }
  
  if (beta==1){
    low<- -5
    high<- 5
  }
  
  if (beta==-1){
    low<- -17
    high<- 5
  }
    pdf(paste("5) Simulations/results2/beta=",beta,".pdf"),width=10,height = 5)
    par (mfrow=c(1,2))
    t1<-t(deviation/apply(deviation,1,sd))
    t2<-t(deviation2/apply(deviation2,1,sd))
    boxplot(t1,ylim=c(low,high),main="Regression-based",ylab=paste("normalized difference, beta=",beta,sep=""),outline=F)
    
    abline(h=0)
    abline(h=1.96,lty=2)
    abline(h=-1.96,lty=2)
    boxplot(t2,ylim=c(low,high),main="IPFP",ylab=paste("normalized difference, beta=",beta,sep=""),outline=F)
    
    abline(h=0)
    abline(h=1.96,lty=2)
    abline(h=-1.96,lty=2)
    dev.off()
  
  
  
  # Save the results----
  save.image(paste("5) Simulations/results2/beta=",beta,".RData"))
  
}






