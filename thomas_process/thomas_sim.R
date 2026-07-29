source('thomas_process/thomas_sim_functions.R')
nsim<-100
mu<-30
nu<-10
sig2<-0.0025
wsize<-1

set.seed(0)
simulate_thomas_data(nsim=nsim,mu,nu,sig2,wsize)

R<-0.2
dr<-0.1/500
rho_sd<-100

model<-stan_model('thomas_process/thomas_palm.stan')


S<-rThomas(nsim=1,kappa=mu,mu=nu,scale=sqrt(sig2),win = owin(c(0,wsize),c(0,wsize)))


