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

# ------------------------------------------------------------------------------#
#### SIMULATION SETTING (1) -- expected 300 points, domain: [0,1]^2
rho<-300
mu<-log(rho)-sig2/2
wsize<-1

#### SIMULATE DATA
set.seed(0) # set seed
simulate_lgcp_data(nsim,mu,sig2,phi,wsize)

#### FIT PALM MODELS
# R = 0.2
R<-0.2
set.seed(0) # set seed
palm_lgcp_sim_study(nsim=10,mu,sig2,phi,wsize,R,nboot=100,empirical=TRUE)

set.seed(0) # set seed
disc_palm_lgcp_sim_study(nsim=nsim,mu,sig2,phi,wsize,R,nboot=100,empirical=TRUE)

# R = 0.4
R<-0.4
set.seed(0) # set seed
palm_lgcp_sim_study(nsim,mu,sig2,phi,wsize,R,nboot=100,empirical=TRUE)

set.seed(0) # set seed
disc_palm_lgcp_sim_study(nsim=nsim,mu,sig2,phi,wsize,R,nboot=100,empirical=TRUE)

#### FIT FULL LIKELIHOOD MODEL
full_lgcp_sim_study(nsim=nsim,mu,sig2,phi,wsize)


# ------------------------------------------------------------------------------#
#### SIMULATION SETTING (2) -- expected 1200 points, domain: [0,2]^2
rho<-300
mu<-log(rho)-sig2/2
wsize<-2

#### SIMULATE DATA
set.seed(0) # set seed
simulate_lgcp_data(nsim,mu,sig2,phi,wsize)

#### FIT PALM MODELS
# R = 0.2
R<-0.2
set.seed(0) # set seed
palm_lgcp_sim_study(nsim=nsim,mu,sig2,phi,wsize,R,nboot=100,empirical=TRUE)

set.seed(0) # set seed
disc_palm_lgcp_sim_study(nsim=nsim,mu,sig2,phi,wsize,R,nboot=100,empirical=TRUE)

# R = 0.4
R<-0.4
set.seed(0) # set seed
palm_lgcp_sim_study(nsim,mu,sig2,phi,wsize,R,nboot=100,empirical=TRUE)

set.seed(0) # set seed
disc_palm_lgcp_sim_study(nsim=nsim,mu,sig2,phi,wsize,R,nboot=100,empirical=TRUE)


# ------------------------------------------------------------------------------#
#### SIMULATION SETTING (3) -- expected 1200 points, domain: [0,1]^2
rho<-1200
mu<-log(rho)-sig2/2
wsize<-1

#### FIT PALM MODELS
# R = 0.2
R<-0.2
set.seed(0) # set seed
palm_lgcp_sim_study(nsim=nsim,mu,sig2,phi,wsize,R,nboot=100,empirical=TRUE)

set.seed(0) # set seed
disc_palm_lgcp_sim_study(nsim=nsim,mu,sig2,phi,wsize,R,nboot=100,empirical=TRUE)

# R = 0.4
R<-0.4
set.seed(0) # set seed
palm_lgcp_sim_study(nsim,mu,sig2,phi,wsize,R,nboot=100,empirical=TRUE)

set.seed(0) # set seed
disc_palm_lgcp_sim_study(nsim=nsim,mu,sig2,phi,wsize,R,nboot=100,empirical=TRUE)

