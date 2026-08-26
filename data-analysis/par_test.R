source('helper_functions.R')
library(RANN)
library(spatstat)
library(fields)
library(tidyverse)
library(rstan)
library(doParallel)
registerDoParallel(cores=10)

temp<-readRDS('data-analysis/bei_files.rds')


Z<-stan_data$t_X[,1:5000]
ZC<-stan_data$mat_counts[1:5000,]
system.time({
  foreach(k=c(1:5))%dopar%{
    t<-system.time({
      for(i in 1:100){
        (exp(c(-5.7,0.164,0.255)%*%Z))%*%(ZC%*%exp(sig2*exp(-stan_data$dG/60)))
      }
    })
    t
  }
})


