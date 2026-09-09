#-----------------------------------------------------------------------#
# EXAMPLE CODE
# This script generates point patterns under an
# LGCP, Thomas process, and determinantal point process,
# and then fits each point pattern using our Bayesian
# Palm likelihood approach in Stan
#-----------------------------------------------------------------------#


# Load necessary packages
library(RANN)
library(tidyverse)
library(spatstat)
library(sf)
library(fields)
library(rstan)
library(doParallel)
registerDoParallel(cores=10)


source('helper_functions.R')
source('example_wrappers.R')

#-----------------------------------------------------------------------#
# LGCP SIMULATION

# Set parameters
rho<-300 # expected number of points
sig2<-1 # variance
phi<-0.1 # range parameter
mu<-log(rho)-sig2/2
wsize<-1 # we assume a square region of observation with side length wsize

set.seed(25)
S<-rLGCP(model="exponential",mu=mu,param=list(var=sig2,scale=phi),win = owin(c(0,wsize),c(0,wsize)))
lgcp_out<-fit_palm_lgcp(S,R=0.2,wsize,empirical=TRUE)

# Plot resulting trace plots
pdf("lgcp_example_traceplots.pdf")
par(mfrow=c(3,3))
plot(lgcp_out$palm_post[,1],type="l",main="mu",ylab="Unadjusted")
abline(h=mu,col="red")
plot(lgcp_out$palm_post[,2],type="l",main="log(sig2)",ylab="")
abline(h=log(sig2),col="red")
plot(lgcp_out$palm_post[,3],type="l",main="log(phi)",ylab="")
abline(h=log(phi),col="red")

plot(lgcp_out$cal1_post[,1],type="l",ylab="Adj. 1")
abline(h=mu,col="red")
plot(lgcp_out$cal1_post[,2],type="l",ylab="")
abline(h=log(sig2),col="red")
plot(lgcp_out$cal1_post[,3],type="l",ylab="")
abline(h=log(phi),col="red")

plot(lgcp_out$cal2_post[,1],type="l",ylab="Adj. 2",xlab="iterations")
abline(h=mu,col="red")
plot(lgcp_out$cal2_post[,2],type="l",ylab="",xlab="iterations")
abline(h=log(sig2),col="red")
plot(lgcp_out$cal2_post[,3],type="l",ylab="",xlab="iterations")
abline(h=log(phi),col="red")
dev.off()

#-----------------------------------------------------------------------#
# THOMAS PROCESS SIMULATION

#### Set parameters
mu<-30 # parent intensity
nu<-10 # offspring intensity
sig2<-0.0025 # variance of bivariate Gaussian
rho<-mu*nu # expected number of points
wsize<-1 # we assume a square region of observation with side length wsize

set.seed(25)
S<-rThomas(nsim=1,kappa=mu,mu=nu,scale=sqrt(sig2),win = owin(c(0,wsize),c(0,wsize)))
thomas_out<-fit_palm_thomas(S,R=0.2,wsize,empirical=FALSE)

# Plot resulting trace plots
pdf("thomas_example_traceplots.pdf")
par(mfrow=c(3,3))
plot(thomas_out$palm_post[,1],type="l",main="log(mu)",ylab="Unadjusted")
abline(h=log(mu),col="red")
plot(thomas_out$palm_post[,2],type="l",main="log(nu)",ylab="")
abline(h=log(nu),col="red")
plot(thomas_out$palm_post[,3],type="l",main="log(sig2)",ylab="")
abline(h=log(sig2),col="red")

plot(thomas_out$cal1_post[,1],type="l",ylab="Adj. 1")
abline(h=log(mu),col="red")
plot(thomas_out$cal1_post[,2],type="l",ylab="")
abline(h=log(nu),col="red")
plot(thomas_out$cal1_post[,3],type="l",ylab="")
abline(h=log(sig2),col="red")

plot(thomas_out$cal2_post[,1],type="l",ylab="Adj. 2",xlab="iterations")
abline(h=log(mu),col="red")
plot(thomas_out$cal2_post[,2],type="l",ylab="",xlab="iterations")
abline(h=log(nu),col="red")
plot(thomas_out$cal2_post[,3],type="l",ylab="",xlab="iterations")
abline(h=log(sig2),col="red")
dev.off()

#-----------------------------------------------------------------------#
# DPP SIMULATION

#### Set parameters
rho<-75 # expected number of points per unit square
lalpha<-log(0.05) # range parameter in Gaussian kernel
wsize<-2 # we assume a square region of observation with side length wsize

set.seed(25)
mG<-dppGauss(lambda=rho,alpha=exp(lalpha),d=2)
S<-simulate.dppm(nsim=1,mG,W=owin(c(0,wsize),c(0,wsize)))
dpp_out<-fit_palm_dpp(S,R=0.2,wsize,empirical=FALSE)

# Plot resulting trace plots
pdf("dpp_example_traceplots.pdf")
par(mfrow=c(3,2))
plot(dpp_out$palm_post[,1],type="l",main="rho",ylab="Unadjusted")
abline(h=rho,col="red")
plot(dpp_out$palm_post[,2],type="l",main="log(alpha)",ylab="")
abline(h=lalpha,col="red")

plot(dpp_out$cal1_post[,1],type="l",ylab="Adj. 1")
abline(h=rho,col="red")
plot(dpp_out$cal1_post[,2],type="l",ylab="")
abline(h=lalpha,col="red")

plot(dpp_out$cal2_post[,1],type="l",ylab="Adj. 2",xlab="iterations")
abline(h=rho,col="red")
plot(dpp_out$cal2_post[,2],type="l",ylab="",xlab="iterations")
abline(h=lalpha,col="red")
dev.off()
