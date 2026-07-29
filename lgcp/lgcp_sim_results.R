
# SIMULATION (1) RESULTS
wsize<-1
sig2<-1
phi<-0.1
mu<-log(300)-1/2
S_all<-readRDS(paste0('sim_output/lgcp/sim_data_',round(mu),'_',wsize,'.rds'))
disc_palm_out_0.2<-readRDS(paste0('sim_output/lgcp/disc_palm_output_',0.2,'_',500,'_',1,'_',round(mu),'_',wsize,'.rds'))
disc_palm_out_0.2_NE<-readRDS(paste0('sim_output/lgcp/disc_palm_output_',0.2,'_',500,'_',100,'_',round(mu),'_',wsize,'.rds'))
disc_palm_out_0.4<-readRDS(paste0('sim_output/lgcp/disc_palm_output_',0.4,'_',500,'_',1,'_',round(mu),'_',wsize,'.rds'))
disc_palm_out_0.4_NE<-readRDS(paste0('sim_output/lgcp/disc_palm_output_',0.4,'_',500,'_',100,'_',round(mu),'_',wsize,'.rds'))
full_out<-readRDS(paste0('sim_output/lgcp/full_output_',round(mu),'_',wsize,'.rds'))


# ------------------------------------------------------------------------------#
### POINT ESTIMATION
# BIASES FROM POSTERIOR MEANS
round(apply(t(matrix(unlist(lapply(full_out, \(x) x$full_post_output[,1]-c(mu,log(sig2),log(phi)))),nrow=3)),2,mean),2)
round(apply(t(matrix(unlist(lapply(disc_palm_out_0.2, \(x) x$cal2_post_output[,1]-c(mu,log(sig2),log(phi)))),nrow=3)),2,mean),2)
round(apply(t(matrix(unlist(lapply(disc_palm_out_0.2_NE, \(x) x$cal2_post_output[,1]-c(mu,log(sig2),log(phi)))),nrow=3)),2,mean),2)
round(apply(t(matrix(unlist(lapply(disc_palm_out_0.4, \(x) x$cal2_post_output[,1]-c(mu,log(sig2),log(phi)))),nrow=3)),2,mean),2)
round(apply(t(matrix(unlist(lapply(disc_palm_out_0.4_NE, \(x) x$cal2_post_output[,1]-c(mu,log(sig2),log(phi)))),nrow=3)),2,mean),2)


# RMSES FROM POSTERIOR MEANS
round(sqrt(apply(t(matrix(unlist(lapply(full_out, \(x) (x$full_post_output[,1]-c(mu,log(sig2),log(phi)))^2)),nrow=3)),2,mean)),2)
round(sqrt(apply(t(matrix(unlist(lapply(disc_palm_out_0.2, \(x) (x$cal2_post_output[,1]-c(mu,log(sig2),log(phi)))^2)),nrow=3)),2,mean)),2)
round(sqrt(apply(t(matrix(unlist(lapply(disc_palm_out_0.2_NE, \(x) (x$cal2_post_output[,1]-c(mu,log(sig2),log(phi)))^2)),nrow=3)),2,mean)),2)
round(sqrt(apply(t(matrix(unlist(lapply(disc_palm_out_0.4, \(x) (x$cal2_post_output[,1]-c(mu,log(sig2),log(phi)))^2)),nrow=3)),2,mean)),2)
round(sqrt(apply(t(matrix(unlist(lapply(disc_palm_out_0.4_NE, \(x) (x$cal2_post_output[,1]-c(mu,log(sig2),log(phi)))^2)),nrow=3)),2,mean)),2)
# ------------------------------------------------------------------------------#

# ------------------------------------------------------------------------------#
### EMPIRICAL COVERAGE RATES
apply(t(matrix(unlist(lapply(full_out, \(x) (x$full_post_output[,4]<c(mu,log(sig2),log(phi)))*(x$full_post_output[,5]>c(mu,log(sig2),log(phi))))),nrow=3)),2,mean)
apply(t(matrix(unlist(lapply(disc_palm_out_0.2, \(x) (x$palm_post_output[,4]<c(mu,log(sig2),log(phi)))*(x$palm_post_output[,5]>c(mu,log(sig2),log(phi))))),nrow=3)),2,mean)
apply(t(matrix(unlist(lapply(disc_palm_out_0.2_NE, \(x) (x$palm_post_output[,4]<c(mu,log(sig2),log(phi)))*(x$palm_post_output[,5]>c(mu,log(sig2),log(phi))))),nrow=3)),2,mean)
apply(t(matrix(unlist(lapply(disc_palm_out_0.4, \(x) (x$palm_post_output[,4]<c(mu,log(sig2),log(phi)))*(x$palm_post_output[,5]>c(mu,log(sig2),log(phi))))),nrow=3)),2,mean)
apply(t(matrix(unlist(lapply(disc_palm_out_0.4_NE, \(x) (x$palm_post_output[,4]<c(mu,log(sig2),log(phi)))*(x$palm_post_output[,5]>c(mu,log(sig2),log(phi))))),nrow=3)),2,mean)

apply(t(matrix(unlist(lapply(disc_palm_out_0.2, \(x) (x$cal1_post_output[,4]<c(mu,log(sig2),log(phi)))*(x$cal1_post_output[,5]>c(mu,log(sig2),log(phi))))),nrow=3)),2,mean)
apply(t(matrix(unlist(lapply(disc_palm_out_0.2_NE, \(x) (x$cal1_post_output[,4]<c(mu,log(sig2),log(phi)))*(x$cal1_post_output[,5]>c(mu,log(sig2),log(phi))))),nrow=3)),2,mean)
apply(t(matrix(unlist(lapply(disc_palm_out_0.4, \(x) (x$cal1_post_output[,4]<c(mu,log(sig2),log(phi)))*(x$cal1_post_output[,5]>c(mu,log(sig2),log(phi))))),nrow=3)),2,mean)
apply(t(matrix(unlist(lapply(disc_palm_out_0.4_NE, \(x) (x$cal1_post_output[,4]<c(mu,log(sig2),log(phi)))*(x$cal1_post_output[,5]>c(mu,log(sig2),log(phi))))),nrow=3)),2,mean)

apply(t(matrix(unlist(lapply(disc_palm_out_0.2, \(x) (x$cal2_post_output[,2]<c(mu,log(sig2),log(phi)))*(x$cal2_post_output[,3]>c(mu,log(sig2),log(phi))))),nrow=3)),2,mean)
apply(t(matrix(unlist(lapply(disc_palm_out_0.2_NE, \(x) (x$cal2_post_output[,2]<c(mu,log(sig2),log(phi)))*(x$cal2_post_output[,3]>c(mu,log(sig2),log(phi))))),nrow=3)),2,mean)
apply(t(matrix(unlist(lapply(disc_palm_out_0.4, \(x) (x$cal2_post_output[,2]<c(mu,log(sig2),log(phi)))*(x$cal2_post_output[,3]>c(mu,log(sig2),log(phi))))),nrow=3)),2,mean)
apply(t(matrix(unlist(lapply(disc_palm_out_0.4_NE, \(x) (x$cal2_post_output[,2]<c(mu,log(sig2),log(phi)))*(x$cal2_post_output[,3]>c(mu,log(sig2),log(phi))))),nrow=3)),2,mean)

# CI LENGTH
round(apply(t(matrix(unlist(lapply(full_out, \(x) (x$full_post_output[,5]-x$full_post_output[,4]))),nrow=3)),2,mean),2)
round(apply(t(matrix(unlist(lapply(disc_palm_out_0.2, \(x) x$palm_post_output[,5]-x$palm_post_output[,4])),nrow=3)),2,mean),2)
round(apply(t(matrix(unlist(lapply(disc_palm_out_0.2_NE, \(x) x$palm_post_output[,5]-x$palm_post_output[,4])),nrow=3)),2,mean),2)
round(apply(t(matrix(unlist(lapply(disc_palm_out_0.4, \(x) x$palm_post_output[,5]-x$palm_post_output[,4])),nrow=3)),2,mean),2)
round(apply(t(matrix(unlist(lapply(disc_palm_out_0.4_NE, \(x) x$palm_post_output[,5]-x$palm_post_output[,4])),nrow=3)),2,mean),2)

round(apply(t(matrix(unlist(lapply(disc_palm_out_0.2, \(x) x$cal1_post_output[,5]-x$cal1_post_output[,4])),nrow=3)),2,mean),2)
round(apply(t(matrix(unlist(lapply(disc_palm_out_0.2_NE, \(x) x$cal1_post_output[,5]-x$cal1_post_output[,4])),nrow=3)),2,mean),2)
round(apply(t(matrix(unlist(lapply(disc_palm_out_0.4, \(x) x$cal1_post_output[,5]-x$cal1_post_output[,4])),nrow=3)),2,mean),2)
round(apply(t(matrix(unlist(lapply(disc_palm_out_0.4_NE, \(x) x$cal1_post_output[,5]-x$cal1_post_output[,4])),nrow=3)),2,mean),2)

round(apply(t(matrix(unlist(lapply(disc_palm_out_0.2, \(x) x$cal2_post_output[,3]-x$cal2_post_output[,2])),nrow=3)),2,mean),2)
round(apply(t(matrix(unlist(lapply(disc_palm_out_0.2_NE, \(x) x$cal2_post_output[,3]-x$cal2_post_output[,2])),nrow=3)),2,mean),2)
round(apply(t(matrix(unlist(lapply(disc_palm_out_0.4, \(x) x$cal2_post_output[,3]-x$cal2_post_output[,2])),nrow=3)),2,mean),2)
round(apply(t(matrix(unlist(lapply(disc_palm_out_0.4_NE, \(x) x$cal2_post_output[,3]-x$cal2_post_output[,2])),nrow=3)),2,mean),2)
# ------------------------------------------------------------------------------#


# SIMULATION (2) RESULTS
wsize<-1
mu<-log(1200)-1/2
S_all_Sim2<-readRDS(paste0('sim_output/lgcp/sim_data_',round(mu),'_',wsize,'.rds'))
disc_palm_out_0.2_Sim2<-readRDS(paste0('sim_output/lgcp/disc_palm_output_',0.2,'_',500,'_',1,'_',round(mu),'_',wsize,'.rds'))
disc_palm_out_0.4_Sim2<-readRDS(paste0('sim_output/lgcp/disc_palm_output_',0.4,'_',500,'_',1,'_',round(mu),'_',wsize,'.rds'))
full_out_Sim2<-readRDS(paste0('sim_output/lgcp/full_output_',round(mu),'_',wsize,'.rds'))


# ------------------------------------------------------------------------------#
### COMPUTATION TIMES IN MINUTES
# RUNTIMES
full_time<-round(mean(unlist(lapply(full_out, \(x) x$time/60))),2)
disc_palm_time_0.2<-round(apply(t(matrix(unlist(lapply(disc_palm_out_0.2[1:10], \(x) c(x$init_time,x$cal2_time,x$init_time+x$cal2_time))),nrow=3))/60,2,mean),2)
disc_palm_time_0.4<-round(apply(t(matrix(unlist(lapply(disc_palm_out_0.4[1:10], \(x) c(x$init_time,x$cal2_time,x$init_time+x$cal2_time))),nrow=3))/60,2,mean),2)

# ESS
full_ESS<-round(apply(t(matrix(unlist(lapply(palm_out, \(x) x$palm_post_output[,6]/(x$init_time+x$cal2_time))),nrow=3)),2,mean),2)
disc_0.2_ESS<-round(apply(t(matrix(unlist(lapply(disc_palm_out_0.2, \(x) x$palm_post_output[,6]/(x$init_time+x$cal2_time))),nrow=3)),2,mean),2)
disc_0.4_ESS<-round(apply(t(matrix(unlist(lapply(disc_palm_out_0.4, \(x) x$palm_post_output[,6]/(x$init_time+x$cal2_time))),nrow=3)),2,mean),2)
# ------------------------------------------------------------------------------#






#### RHO=1200
rho_sd<-1
rho<-1200
mu<-log(rho)-sig2/2
wsize<-1

disc_palm_out_0.2<-readRDS(paste0('sim_output/lgcp/disc_palm_output_',0.2,'_',500,'_',round(rho_sd),'_',round(mu),'_',wsize,'.rds'))
disc_palm_out_0.4<-readRDS(paste0('sim_output/lgcp/disc_palm_output_',0.4,'_',500,'_',round(rho_sd),'_',round(mu),'_',wsize,'.rds'))

apply(t(matrix(unlist(lapply(disc_palm_out_0.2, \(x) x$cal2_post_output[,1]-c(mu,log(sig2),log(phi)))),nrow=3)),2,mean)
apply(t(matrix(unlist(lapply(disc_palm_out_0.4, \(x) x$cal2_post_output[,1]-c(mu,log(sig2),log(phi)))),nrow=3)),2,mean)


# ------------------------------------------------------------------------------#
#### APPENDIX COMPARISON
palm_out<-readRDS(paste0('sim_output/lgcp/palm_output_',0.2,'_',round(rho_sd),'_',round(mu),'_',wsize,'.rds'))
disc_palm_out_0.2_250<-readRDS(paste0('sim_output/lgcp/disc_palm_output_',0.2,'_',250,'_',round(rho_sd),'_',round(mu),'_',wsize,'.rds'))
disc_palm_out_0.2_750<-readRDS(paste0('sim_output/lgcp/disc_palm_output_',0.2,'_',750,'_',round(rho_sd),'_',round(mu),'_',wsize,'.rds'))

palm_time<-round(apply(t(matrix(unlist(lapply(palm_out, \(x) c(x$init_time,x$cal2_time,x$init_time+x$cal2_time))),nrow=3))/60,2,mean),2)

palm_means<-t(matrix(unlist(lapply(palm_out, \(x) x$cal2_post_output[,1])),nrow=3))
disc_palm_250_means<-t(matrix(unlist(lapply(disc_palm_out_0.2_250, \(x) x$cal2_post_output[,1])),nrow=3))
disc_palm_500_means<-t(matrix(unlist(lapply(disc_palm_out_0.2[1:10], \(x) x$cal2_post_output[,1])),nrow=3))
disc_palm_750_means<-t(matrix(unlist(lapply(disc_palm_out_0.2_750, \(x) x$cal2_post_output[,1])),nrow=3))

par(mfrow=c(3,3))
plot(palm_means[,1],disc_palm_250_means[,1])
plot(palm_means[,1],disc_palm_500_means[,1])
plot(palm_means[,1],disc_palm_750_means[,1])

apply(t(matrix(unlist(lapply(1:10,\(x) abs(palm_out[[x]]$palm_post_output[,1]-disc_palm_out_0.2_250[[x]]$palm_post_output[,1]))),nrow=3)),2,mean)
apply(t(matrix(unlist(lapply(1:10,\(x) abs(palm_out[[x]]$palm_post_output[,1]-disc_palm_out_0.2[[x]]$palm_post_output[,1]))),nrow=3)),2,mean)
apply(t(matrix(unlist(lapply(1:10,\(x) abs(palm_out[[x]]$palm_post_output[,1]-disc_palm_out_0.2_750[[x]]$palm_post_output[,1]))),nrow=3)),2,mean)

apply(t(matrix(unlist(lapply(1:10,\(x) abs(palm_out[[x]]$palm_post_output[,4]-disc_palm_out_0.2_250[[x]]$palm_post_output[,4]))),nrow=3)),2,mean)
apply(t(matrix(unlist(lapply(1:10,\(x) abs(palm_out[[x]]$palm_post_output[,4]-disc_palm_out_0.2[[x]]$palm_post_output[,4]))),nrow=3)),2,mean)
apply(t(matrix(unlist(lapply(1:10,\(x) abs(palm_out[[x]]$palm_post_output[,4]-disc_palm_out_0.2_750[[x]]$palm_post_output[,4]))),nrow=3)),2,mean)

apply(t(matrix(unlist(lapply(1:10,\(x) abs(palm_out[[x]]$palm_post_output[,5]-disc_palm_out_0.2_250[[x]]$palm_post_output[,5]))),nrow=3)),2,mean)
apply(t(matrix(unlist(lapply(1:10,\(x) abs(palm_out[[x]]$palm_post_output[,5]-disc_palm_out_0.2[[x]]$palm_post_output[,5]))),nrow=3)),2,mean)
apply(t(matrix(unlist(lapply(1:10,\(x) abs(palm_out[[x]]$palm_post_output[,5]-disc_palm_out_0.2_750[[x]]$palm_post_output[,5]))),nrow=3)),2,mean)


# ------------------------------------------------------------------------------#



