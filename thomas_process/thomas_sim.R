source('thomas_process/thomas_sim_functions.R')
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
mu<-30
nu<-10
sig2<-0.0025
rho<-mu*nu
wsize<-1

#### SIMULATE DATA
set.seed(0)
#simulate_thomas_data(nsim=nsim,mu,nu,sig2,wsize)

#### FIT PALM MODELS
# R = 0.2
R<-0.2
set.seed(0) # set seed
disc_palm_thomas_sim_study(nsim=nsim,mu,sig2,phi,wsize,R,nboot=100,empirical=TRUE) 

set.seed(0) # set seed
disc_palm_thomas_sim_study(nsim=nsim,mu,sig2,phi,wsize,R,nboot=100,empirical=FALSE)
# ------------------------------------------------------------------------------#

