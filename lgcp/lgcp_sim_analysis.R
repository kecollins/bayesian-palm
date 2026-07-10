#### SIMULATION SETTINGS
rho<-300
sig2<-1
phi<-0.1
mu<-log(rho)-sig2/2

wsize<-1
nsim<-2

sim_output<-readRDS(paste0('sim_output/lgcp/sim_output_',R,'_',1,'_',round(mu),'_',wsize,'.rds'))

sim_output[[2]]$init_time
sim_output[[2]]$simulate_time+sim_output[[2]]$cal_time
