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

model<-stan_model('thomas-process/thomas_palm.stan')


S<-rThomas(nsim=1,kappa=mu,mu=nu,scale=sqrt(sig2),win = owin(c(0,wsize),c(0,wsize)))
# ORGANIZE DATA
data<-data_clean(S=S,R=R,wsize,dr=dr)

# SETUP STAN DATA
stan_data<-list(
  rho_mean=S$n/area(owin(c(0,wsize),c(0,wsize))),
  rho_sd=rho_sd,
  dN=length(data$dG),
  dG=data$dG,
  dG_counts=data$dG_counts,
  dS_counts=data$dS_counts,
  dx=data$dx,
  eta=1
)

post <- sampling(
  model,  # Stan program compiled earlier
  data = stan_data,    # named list of data
  chains = 1,             # number of Markov chains
  warmup = burn,          # number of warmup iterations per chain
  iter = iters,            # total number of iterations per chain
  cores = 1,              # number of cores (could use one per chain)
  refresh = 0,             # no progress shown
  init=list(list(rho=rho,lsig2=0,lnu=10))#,
  #control = list(adapt_delta = 0.99)
)

traceplot(post)
kppm(S,cluster="Thomas",method="palm",rmax=0.4)




nu/(4*pi^2*sig2)exp(-stan_data$dG^2/(4*sig2))

