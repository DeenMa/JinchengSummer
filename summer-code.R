

#####################################################################
rm(list=ls())
library(MCMCpack)
library(bvsflex)
#####################################################################

#####################################################################
estimate <- function(E,T,X,D,Z,mcmc.samples=8000,burn.in = 3000,prior.sd,probability){

### Y=(E,T), efficacy, toxi response
### X=covariate matrix, n by p. no intercept
### Z[i]denote the treatment assignment,0->control; 1->treatment
### D[i]be the treatment assignment that patient i actually received
   X.orig<-matrix(1,nrow=nrow(X),ncol=((ncol(X)-2)*4+2))
   X.orig[,1]<-X[,1]
   X.orig[,2]<-X[,2]
   for(i in 3:ncol(X)){
      X.orig[,((i-2)*4-1):((i-2)*4+2)]<-cbind(X[,i]^(1/2),X[,i]^(1),X[,i]^(3/2),X[,i]^(2))
   }
   X.true<-matrix(1,nrow=nrow(X),ncol=sum(probability))
   tempp=1
   for(i in 1:ncol(X.orig)){
   if(probability[i]==1){
     X.true[,tempp]<-X.orig[,i]
     tempp=tempp+1
   }
   }
   p.beta <- 3*ncol(X.true)
   p.beta.all <- 3*ncol(X.orig) 
   n <- nrow(X)      #sample size
   y <- 1 + 2*E + T  # y=1,2,3,4 (0,0),(0,1),(1,0),(1,1),note that y is vector

## figure out the compliance status
## s= -1,missing; 0,n; 1,c; 2,a

s <- rep(-1,n)  #Set Initial Value

for(i in 1:n){
   if(Z[i]==1 & D[i]==0) s[i] <- 0 
   if(Z[i]==0 & D[i]==1) s[i] <- 2
}

miss.loc <- which(s==-1) #the index of missing value
   
   prior.mn <- rep(0,p.beta)

   n.miss <- length(miss.loc)
   p <- p.beta.all + 4*6 + n.miss 

   store.par<-matrix(0,mcmc.samples,p)       #store parameters
   acc<-rep(0,p.beta)                        #acceptance ratio monitor
   a.prior <- sum(D[Z==0])/length(D)    
   n.prior <- sum(D[Z==1]==0)/length(D)                                                                     
   c.prior <- 1 - a.prior - n.prior     
    # Initial values: 
       beta<-rnorm(p.beta,0,1)
       s.miss <- NULL
       for(i in 1:n){
       if(s[i]==-1 & Z[i]==0 & D[i]==0){
         s[i] <- rbinom(1,1,c.prior/(c.prior+n.prior))
         s.miss <- c(s.miss,s[i])
       }
       if(s[i]==-1 & Z[i]==1 & D[i]==1){
         s[i] <- 1 + rbinom(1,1,a.prior/(a.prior+c.prior))
         s.miss <- c(s.miss,s[i])
       }
       }
       theta.c.trt <- rdirichlet(1,rep(1,6)) #set independent Dirichlet prior
       theta.c.con <- rdirichlet(1,rep(1,6)) 
       theta.n.trt <- rdirichlet(1,rep(1,6))
       theta.n.con <- rdirichlet(1,rep(1,6))
       theta.a.trt <- rdirichlet(1,rep(1,6))
       theta.a.con <- rdirichlet(1,rep(1,6))

   for(i in 1:mcmc.samples){ 
   
    # update theta
    # y.1 <- y[s==1]
    # y.0 <- y[s==0]
    # y.2 <- y[s==2]
     y.1.trt <- y[Z==1 & s==1]
     y.1.con <- y[Z==0 & s==1]
     y.0.trt <- y[Z==1 & s==0]
     y.0.con <- y[Z==0 & s==0]
     y.2.trt <- y[Z==1 & s==2]
     y.2.con <- y[Z==0 & s==2]

     theta.u1.trt <- c(sum(y.1.trt==1),sum(y.1.trt==2),sum(y.1.trt==3),sum(y.1.trt==4)) + 1 
     theta.u1.con <- c(sum(y.1.con==1),sum(y.1.con==2),sum(y.1.con==3),sum(y.1.con==4)) + 1
     theta.u0.trt <- c(sum(y.0.trt==1),sum(y.0.trt==2),sum(y.0.trt==3),sum(y.0.trt==4)) + 1
     theta.u0.con <- c(sum(y.0.con==1),sum(y.0.con==2),sum(y.0.con==3),sum(y.0.con==4)) + 1
     theta.u2.trt <- c(sum(y.2.trt==1),sum(y.2.trt==2),sum(y.2.trt==3),sum(y.2.trt==4)) + 1
     theta.u2.con <- c(sum(y.2.con==1),sum(y.2.con==2),sum(y.2.con==3),sum(y.2.con==4)) + 1
     
     theta.c.trt <- rdirichlet(1,theta.u1.trt) #sampling theta(s,z)
     theta.c.con <- rdirichlet(1,theta.u1.con)
     theta.n.trt <- rdirichlet(1,theta.u0.trt)
     theta.n.con <- rdirichlet(1,theta.u0.con)
     theta.a.trt <- rdirichlet(1,theta.u2.trt)
     theta.a.con <- rdirichlet(1,theta.u2.con)
     
     theta <- c(theta.c.trt,theta.c.con,theta.n.trt,theta.n.con,theta.a.trt,theta.a.con) 
     theta.c <- c(theta.c.trt,theta.c.con)
     theta.n <- c(theta.n.trt,theta.n.con)
     theta.a <- c(theta.a.trt,theta.a.con)

     store.par[i,(p.beta.all+1):(p.beta.all+24)]<-theta #save the theta 
    
#############################################################    

    #Update beta using MH sampling: Gibbs sampling
     for(j in 1:p.beta){           #sampling beta
       canbeta<-beta
       canbeta[j]<-rnorm(1,beta[j],1) 
       post1 <- post(X.true,s,canbeta,0,prior.sd)
       post2 <- post(X.true,s,beta,0,prior.sd)
       R<-exp(post1-post2)   #compute acceptance ratio
       U<-runif(1)                          
       if(U<R){                    #Accept the candidate with prob min(R,1)
         beta<-canbeta
         acc[j]<-acc[j]+1
       }
     }
    beta.all<-rep(0,p.beta.all)
    beta.all[which(probability==1)]<-beta[1:(p.beta/3)]
    beta.all[which(probability==1)+(p.beta.all/3)]<-beta[(p.beta/3+1):(2*p.beta/3)]
    beta.all[which(probability==1)+(2*p.beta.all/3)]<-beta[(2*p.beta/3+1):p.beta]
    acc.all<-rep(0,p.beta.all)
    acc.all[which(probability==1)]<-acc[1:(p.beta/3)]
    acc.all[which(probability==1)+(p.beta.all/3)]<-acc[(p.beta/3+1):(2*p.beta/3)]
    acc.all[which(probability==1)+(2*p.beta.all/3)]<-acc[(2*p.beta/3+1):p.beta]
    store.par[i,1:p.beta.all]<-beta.all    #save the beta

#########################
#########################                                                                                                
xbeta.c <- ifelse(abs(X.true%*%beta[1:(p.beta/3)]) <=10,X.true%*%beta[1:(p.beta/3)],10*sign(X.true%*%beta[1:(p.beta/3)]))  #logistic regression
xbeta.n <- ifelse(abs(X.true%*%beta[(p.beta/3+1):(2*p.beta/3)]) <=10,X.true%*%beta[(p.beta/3+1):(2*p.beta/3)],10*sign(X.true%*%beta[(p.beta/3+1):(2*p.beta/3)]))
xbeta.a <- ifelse(abs(X.true%*%beta[(2*p.beta/3+1):p.beta]) <=10,X.true%*%beta[(2*p.beta/3+1):p.beta],10*sign(X.true%*%beta[(2*p.beta/3+1):p.beta]))
temp.c <- exp(xbeta.c)/(1+exp(xbeta.c)+exp(xbeta.n))
temp.n <- exp(xbeta.n)/(1+exp(xbeta.c)+exp(xbeta.n))
temp.a <- 1 - temp.c - temp.n
#########################
#########################

#############################################################
     # update s.miss
     s.miss <- NULL
     for(w in miss.loc){
      tempc <- temp.c[w]
      tempn <- temp.n[w]
      tempa <- temp.a[w]
      temp.i <- (1-Z[w])*4 + y[w]
      #prob.latent <-  s.prior * temp1 * theta.c[temp.i]/(s.prior * temp1 * theta.c[temp.i]+(1-s.prior) * (1-temp1) * theta.n[temp.i])
      if(Z[w]==0 & D[w]==0){
      p.i.c <- (c.prior/(c.prior+n.prior)) * tempc * theta.c[temp.i]/((c.prior/(c.prior+n.prior)) * tempc * theta.c[temp.i]+(1-(c.prior/(c.prior+n.prior))) * tempn * theta.n[temp.i])
      p.i.c <- ifelse(p.i.c > 1,1,p.i.c) 
      p.i.c <- ifelse(p.i.c < 0,0,p.i.c)
      s[w] <- rbinom(1,1,p.i.c)              
      }
      if(Z[w]==1 & D[w]==1){
      p.i.t <- (a.prior/(a.prior+c.prior)) * tempc * theta.c[temp.i]/((a.prior/(a.prior+c.prior)) * tempc * theta.c[temp.i]+(1-(a.prior/(a.prior+c.prior))) * tempa * theta.a[temp.i])
      p.i.t <- ifelse(p.i.t > 1,1,p.i.t) 
      p.i.t <- ifelse(p.i.t < 0,0,p.i.t)
      s[w] <- 1 + rbinom(1,1,1 - p.i.t) 
      }
     s.miss <- c(s.miss,s[w])
     }
     store.par[i,(p.beta.all+25):p] <- s.miss    #save the missing value
 

   }  # end loop i
   miss.mean <- apply(store.par[(burn.in+1):mcmc.samples,(p.beta.all+25):p]==1,2,mean) 
     
   list(par.all=store.par,par=store.par[(burn.in+1):mcmc.samples,1:(p.beta.all+24)],miss=miss.mean,acc=acc.all/(mcmc.samples+burn.in) ) #output

}
##################################################################

##################################################################
post <- function(X,s,beta,prior.mn,prior.sd){
    #full joint posterior 
    #used for MH sampling of beta
    #return the logarithm
    p.beta<-3*ncol(X)
    xbeta.c <-  ifelse(abs(X%*%beta[1:(p.beta/3)]) <=10,X%*%beta[1:(p.beta/3)],10*sign(X%*%beta[1:(p.beta/3)]))  #logistic regression
    xbeta.n <-  ifelse(abs(X%*%beta[(p.beta/3+1):(2*p.beta/3)]) <=10,X%*%beta[(p.beta/3+1):(2*p.beta/3)],10*sign(X%*%beta[(p.beta/3+1):(2*p.beta/3)]))
    xbeta.a <-  ifelse(abs(X%*%beta[(2*p.beta/3+1):p.beta]) <=10,X%*%beta[(2*p.beta/3+1):p.beta],10*sign(X%*%beta[(2*p.beta/3+1):p.beta]))
    like<-sum(log(((exp(xbeta.c)/(1+exp(xbeta.c)+exp(xbeta.n)))^(I(s==1)))*((exp(xbeta.n)/(1+exp(xbeta.c)+exp(xbeta.n)))^(I(s==0)))*(1-(exp(xbeta.c)/(1+exp(xbeta.c)+exp(xbeta.n)))-(exp(xbeta.n)/(1+exp(xbeta.c)+exp(xbeta.n))))^(I(s==2))))
    prior<-sum(dnorm(beta,prior.mn,prior.sd,log=T))
    like+prior
}
###################################################################

###################################################################
estimate.itt <- function(E,T,Z){         #traditional method
## prior (E,T) Dirichlet (1,1,1,1) for each assignment Z=0 or 1
  y <- 1 + 2*E + T  # y=1,2,3,4 (0,0),(0,1),(1,0),(1,1)
  y.1 <- y[Z==1]
  y.0 <- y[Z==0]
  theta.z1 <- c(length(which(y.1==1)),length(which(y.1==2)),length(which(y.1==3)),length(which(y.1==4))) + 1
  theta.z0 <- c(length(which(y.0==1)),length(which(y.0==2)),length(which(y.0==3)),length(which(y.0==4))) + 1
  matrix1 <- rdirichlet(5000,theta.z1)
  matrix2 <- rdirichlet(5000,theta.z0)
  return(cbind(matrix1,matrix2))
}
####################################################################

####################################################################
estimate.gold <- function(E,T,Z,comply){ #optimal result
## prior (E,T) Dirichlet (1,1,1,1) for each assignment Z=0 or 1 and comply=1
  y <- 1 + 2*E + T  # y=1,2,3,4 (0,0),(0,1),(1,0),(1,1)
  y.1 <- y[Z==1 & comply==1]
  y.0 <- y[Z==0 & comply==1]
  theta.z1 <- c(length(which(y.1==1)),length(which(y.1==2)),length(which(y.1==3)),length(which(y.1==4))) + 1
  theta.z0 <- c(length(which(y.0==1)),length(which(y.0==2)),length(which(y.0==3)),length(which(y.0==4))) + 1
  matrix1 <- rdirichlet(5000,theta.z1)
  matrix2 <- rdirichlet(5000,theta.z0)
  return(cbind(matrix1,matrix2))
}
#####################################################################




#####################################################################
### main function for trial continous monitoring

monitor <- function(n0, n.max, eff.new, eff.old, tox.new, tox.old, eff.cut, tox.max, tox.cut, eff.p, tox.p, beta.true){
### n0: initial recurits
### n.max: max number of patients.

## generate data
## covariates x1,x2


x1 <- rbinom(n.max,size=1,prob=.4)
x2 <- runif(n.max,min=0,max=1)
X <- cbind(rep(1,n.max),x1,x2)

X.true<-matrix(1,nrow=nrow(X),ncol=((ncol(X)-2)*4+2))
X.true[,1]<-X[,1]
X.true[,2]<-X[,2]
for(i in 3:ncol(X)){
   X.true[,((i-2)*4-1):((i-2)*4+2)]<-cbind(X[,i]^(1/2),X[,i]^(1),X[,i]^(3/2),X[,i]^(2))
}

p.beta<-3*ncol(X.true)

##############################
beta <- beta.true  # true beta
xbeta.c <-X.true%*%beta[1:(p.beta/3)]
xbeta.n <-X.true%*%beta[(p.beta/3+1):(2*p.beta/3)]
xbeta.a <-X.true%*%beta[(2*p.beta/3+1):p.beta]
prob.c <- exp(xbeta.c)/(1+exp(xbeta.c)+exp(xbeta.n))
prob.n <- exp(xbeta.n)/(1+exp(xbeta.c)+exp(xbeta.n))
prob.a <- 1 - prob.c - prob.n
##############################

comply <- NULL
c.plus.a <- rbinom(n.max,size=1,prob=prob.c+prob.a)
for(i in 1:n.max){
    if(c.plus.a[i]==0)  comply[i] <- 0
    if(c.plus.a[i]==1)  comply[i] <- c.plus.a[i] + rbinom(1,1,prob.a/(prob.a+prob.c))
}

p<-ncol(X.true)
n<-nrow(X.true)
d <- 3
B0 <- 1000
res<- logisticVS(X.true, comply, b=rep(0,p), h0=rep(1,p), g=1, block=diag(1,p), 
       mu=rep(d/p,p), phi=rep(p/2,p),
       MCMC=B0, thinn=1, seed=1234, 
       outdir=paste(getwd(), "example", sep="/"), 
       Piupdate=FALSE, MuPhiUpdate=FALSE, gupdate="none")
require(Matrix)
tmp <- read.table(paste(res, "GAM.txt", sep="/"))
GAM <- as.matrix(sparseMatrix(i = tmp[,1], j = tmp[,2], x = tmp[,3]))
K<-apply(GAM,1,mean)
probability<-(K>0.5)*1
if(sum(probability)<2){
  probability[order(K,decreasing=TRUE)[1:2]]<-c(1,1)
}

## trt assignment, z = 1, new trt
z <- rbinom(n.max,size=1,prob=0.5)

## actural taken assignment
d <- NULL
for(i in 1:n.max){
    if(comply[i]==2) d[i] <- 1
    if(comply[i]==0) d[i] <- 0
    if(comply[i]==1) d[i] <- z[i]
}

prob.e <-  d * eff.new + (1-d) * eff.old
prob.t <- d * tox.new + (1-d) * tox.old

eff <- rbinom(n.max,size=1,prob=prob.e)
tox <- rbinom(n.max,size=1,prob=prob.t)

earlyt.our <- 0  # early terminiation indicator
earlyt.itt <- 0
earlyt.gold <- 0

eff.coef <- c(0,0,1,1,0,0,-1,-1)
tox.coef <- c(0,1,0,1,0,0,0,0)
num.our <- n.max
num.itt <- n.max
num.gold <- n.max

for(ii in seq(from=n0,to=n.max,by=5) ){
    eff.ii <- eff[1:ii] 
    tox.ii <- tox[1:ii]
    X.ii   <- X[1:ii,]
    d.ii   <- d[1:ii]
    z.ii   <- z[1:ii]
    comply.ii <- comply[1:ii]
    comply.miss <- comply.ii[(z.ii==0 & d.ii==0) | (z.ii==1 & d.ii==1)]

    fit <- estimate(E=eff.ii,T=tox.ii,X=X.ii,D=d.ii,Z=z.ii,prior.sd=5,probability=probability)

    fit.use <- fit$par[,(p.beta+1):(p.beta+8)]   # extracts theta in compliance strata
               
    fit.beta <- apply(fit$par,2,mean)[1:p.beta]   ## fitted beta
    fit.acc <- fit$acc
    fit.eff <- c(fit.use %*% eff.coef)
    fit.max <- c(fit.use %*% tox.coef)

    fit.miss <- 1 - mean(abs( (fit$miss >= 0.5)*1 == comply.miss) )

    fit.theta.our <- apply(fit.use,2,mean)
     print(c("our prob",mean(fit.eff <= eff.cut)))
     print(c("our mean",fit.theta.our%*%eff.coef ))
     print(c("our miss",fit.miss))
   if(mean(fit.eff <= eff.cut) > eff.p){
      num.our <- ii
      earlyt.our <- 1
      cat("early termination due to efficacy(our method)")
      cat("number of patients=")
      cat(ii)
      cat("\n")
      break
   }

   
   if(mean(fit.max >= tox.max) > tox.p){
      num.our <- ii
      earlyt.our <- 1
      cat("early termination due to high toxicity(our method)")
      cat("number of patients=")
      cat(ii)
      cat("\n")
      break
   }
   
}
#################################################################

#################################################################
for(iii in seq(from=n0,to=n.max,by=5)){
    eff.iii <- eff[1:iii]
    tox.iii <- tox[1:iii]
    z.iii   <- z[1:iii]
    fit.itt <- estimate.itt(E=eff.iii,T=tox.iii,Z=z.iii)
 
    fit.eff <- c(fit.itt %*% eff.coef)

    fit.max <- c(fit.itt %*% tox.coef)
    fit.theta.itt <- apply(fit.itt,2,mean)
    print(c("itt prob",mean(fit.eff <= eff.cut)))
    print(c("itt mean",fit.theta.itt%*%eff.coef ))
    if(mean(fit.eff <= eff.cut) > eff.p){
      num.itt <- iii
      earlyt.itt <- 1
      cat("early termination due to efficacy(ITT method)")
      cat("number of patients=")
      cat(iii)
      cat("\n")
      break
     }

    if(mean(fit.max >= tox.max) > tox.p){
      num.itt <- iii
      earlyt.itt <- 1
      cat("early termination due to high toxicity(ITT method)")
      cat("number of patients=")
      cat(iii)
      cat("\n")
      break
   }
   
}
##################################################################

##################################################################
for(i in seq(from=n0,to=n.max,by=5)){
    eff.i <- eff[1:i]
    tox.i <- tox[1:i]
    d.i   <- d[1:i]
    comply.i <- comply[1:i]
    fit.gold <- estimate.gold(E=eff.i,T=tox.i,Z=d.i,comply=comply.i)
 
    
    fit.eff <- c(fit.gold %*% eff.coef)

    fit.max <- c(fit.gold %*% tox.coef)
    fit.theta.gold <- apply(fit.gold,2,mean)
    print(c("gold prob",mean(fit.eff <= eff.cut)))
    print(c("gold mean",fit.theta.gold%*%eff.coef ))
    if(mean(fit.eff <= eff.cut) > eff.p){
      num.gold <- i
      earlyt.gold <- 1
      cat("early termination due to efficacy(gold method)")
      cat("number of patients=")
      cat(i)
      cat("\n")
      break
     }

     if(mean(fit.max >= tox.max) > tox.p){
      num.gold <- i
      earlyt.gold <- 1
      cat("early termination due to high toxicity(gold method)")
      cat("number of patients=")
      cat(i)
      cat("\n")
      break
   }
   
}
####################################################################


   return(list(num.our=num.our,num.itt=num.itt,num.gold=num.gold,earlyt.our=earlyt.our,earlyt.itt=earlyt.itt,earlyt.gold=earlyt.gold,beta.est=fit.beta,diag.acc=fit.acc,theta.our=fit.theta.our,theta.itt=fit.theta.itt,theta.gold=fit.theta.gold,fit.miss =fit.miss ))


}## end of main function 
##################################################################

##################################################################
trial.sim <- function(B,n0,n.max,eff.new,eff.old,tox.new,tox.old,eff.cut,tox.max,eff.p,tox.p,beta.true){
## simulation function, B: number of replications
## n0: number of patients to start with, n.max: max number of patients

p.beta<-length(beta.true)
store <- matrix(0,B,6)  # store results
beta <- matrix(0,B,2*p.beta)   # store diagnostics, accpetance ratio, estimated beta.
theta <- matrix(0,B,8*3)  # store theta means
miss <- rep(0,B)

for(i in 1:B){

  result <- monitor(n0,n.max,eff.new,eff.old,tox.new,tox.old,eff.cut,tox.max,tox.cut,eff.p,tox.p,beta.true)

  store[i,] <- c(result$num.our,result$num.itt,result$num.gold,result$earlyt.our,result$earlyt.itt,result$earlyt.gold)
  beta[i,] <- c(result$beta.est,result$diag.acc)
  theta[i,] <- c(result$theta.our,result$theta.itt,result$theta.gold)
  miss[i] <- result$fit.miss
  }
  return(list(store=store,beta=beta,theta=theta,miss=miss))
}
####################################################################


source("trial_case1_main_speed.R")

set.seed(1234)


#result <- trial.sim(B=1000,n0=30,n.max=150,eff.new=.2,eff.old=.2,tox.new=.1,tox.old=.1,eff.cut=.15,tox.max=.2,eff.p=.9,tox.p=.9,beta.true=c(1.2,-8,3) )
result <- trial.sim(B=100,n0=30,n.max=150,eff.new=.2,eff.old=.2,tox.new=.1,tox.old=.1,eff.cut=.15,tox.max=.2,eff.p=.9,tox.p=.9,beta.true=rnorm(18,0,2) )


save(list=ls(),file = 'trial_case01_1.Rdata')




