#-----------------------------------------------------------------------#
# EXAMPLE CODE
# This script generates point patterns under an
# LGCP, Thomas process, and determinantal point process,
# and then fits each point pattern using our Bayesian
# Palm likelihood approach in Stan
#-----------------------------------------------------------------------#

# Comment out the next line if packages are already installed
install.packages('RANN','tidyverse','spatstat','sf','fields','rstan','doParallel')

# Load necessary packages
library(RANN, tidyverse, spatstat, sf, fields, rstan, doParallel)

# Initialize parallelization
cl<-makeCluster(10) # 10 cores
registerDoParallel(cl)

#-----------------------------------------------------------------------#
# LGCP SIMULATION

# Set parameters
rho<-300 # expected number of points
sig2<-1 # variance
phi<-0.1 # range parameter
mu<-log(rho)-sig2/2
wsize<-2 # we assume a square region of observation with side length wsize

set.seed(0)
S<-rLGCP(model="exponential",mu=mu,param=list(var=sig2,scale=phi),win = owin(c(0,wsize),c(0,wsize)))
lgcp_out<-fit_palm_lgcp(S,R,wsize,empirical=TRUE)

# Plot resulting trace plots after GPC adjustment
par(mfrow=c(1,3))
plot(lgcp_out$cal2_post[,1],type="l")
abline(h=mu,col="red")
plot(lgcp_out$cal2_post[,2],type="l")
abline(h=log(sig2),col="red")
plot(lgcp_out$cal2_post[,3],type="l")
abline(h=log(phi),col="red")

#-----------------------------------------------------------------------#
# THOMAS PROCESS SIMULATION

# Set parameters
rho<-300 # expected number of points
sig2<-1 # variance
phi<-0.1 # range parameter
mu<-log(rho)-sig2/2
wsize<-1 # we assume a square region of observation with side length wsize

set.seed(0)
S<-rLGCP(model="exponential",mu=mu,param=list(var=sig2,scale=phi))
lgcp_out<-fit_palm_lgcp(S,R,wsize,empirical=FALSE)

# Plot resulting trace plots after GPC adjustment
par(mfrow=c(1,3))
plot(lgcp_out$cal2_post[,1],type="l")
abline(h=mu,col="red")
plot(lgcp_out$cal2_post[,2],type="l")
abline(h=log(sig2),col="red")
plot(lgcp_out$cal2_post[,3],type="l")
abline(h=log(phi),col="red")

plot(density(lgcp_out$palm_post[,1]),xlim=c(3,6))
lines(density(lgcp_out$cal2_post[,1]),col="blue")
