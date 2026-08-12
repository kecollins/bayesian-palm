source('dpp/dpp_sim_functions.R')
source('helper_functions.R')
library(RANN)
library(tidyverse)
library(spatstat)
library(sf)
library(fields)
library(rstan)
library(doParallel)
registerDoParallel(cores=10)

# ------------------------------------------------------------------------------#
#### SIMULATION SETTINGS
nsim<-100
rho<-75
lalpha<-log(0.05)
wsize<-2


#### SIMULATE DATA
set.seed(0)
simulate_dpp_data(nsim,rho,lalpha,wsize)

#### FIT PALM MODELS
# R = 0.2
R<-0.2
set.seed(0) # set seed
disc_palm_dpp_sim_study(nsim=nsim,mu,sig2,phi,wsize,R,nboot=100,empirical=TRUE) 

set.seed(0) # set seed
disc_palm_dpp_sim_study(nsim=nsim,mu,sig2,phi,wsize,R,nboot=100,empirical=FALSE)
# ------------------------------------------------------------------------------#

