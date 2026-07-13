source('lgcp/lgcp_sim_functions.R')
source('helper_functions.R')
library(RANN)
library(tidyverse)
library(spatstat)
library(sf)
library(fields)
library(rstan)
library(cmdstanr)
library(doParallel)
cl<-makeCluster(2)
registerDoParallel(cl)
registerDoParallel(cores=10)

#### SIMULATION SETTINGS (UNCHANGING)
nsim<-100
sig2<-1
phi<-0.1

#### SIMULATION SETTING (1)
rho<-300
mu<-log(rho)-sig2/2
wsize<-1

#### SIMULATE DATA
set.seed(0) # set seed
simulate_lgcp_data(nsim,mu,sig2,phi,wsize)

#### FIT PALM MODELS
set.seed(0) # set seed
R<-0.2
palm_lgcp_sim_study(nsim,mu,sig2,phi,wsize,R,nboot=100,empirical=TRUE)

set.seed(0) # set seed
R<-0.4
palm_lgcp_sim_study(nsim,mu,sig2,phi,wsize,R,nboot=100,empirical=TRUE)

#### FIT FULL LIKELIHOOD MODEL
S<-rLGCP(nsim=1000,model="exponential",mu=5.31,param=list(var=1.42,scale=0.13))
mean(unlist(lapply(S, \(x) x$n)))
var(unlist(lapply(S, \(x) x$n)))
hist(unlist(lapply(S, \(x) x$n)))
