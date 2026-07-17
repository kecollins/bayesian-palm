rho_sd<-1
S_all<-readRDS(paste0('sim_output/lgcp/sim_data_',round(mu),'_',wsize,'.rds'))
palm_out<-readRDS(paste0('sim_output/lgcp/palm_output_',0.2,'_',round(rho_sd),'_',round(mu),'_',wsize,'.rds'))
disc_palm_out_0.2_250<-readRDS(paste0('sim_output/lgcp/disc_palm_output_',0.2,'_',250,'_',round(rho_sd),'_',round(mu),'_',wsize,'.rds'))
disc_palm_out_0.2_750<-readRDS(paste0('sim_output/lgcp/disc_palm_output_',0.2,'_',750,'_',round(rho_sd),'_',round(mu),'_',wsize,'.rds'))
disc_palm_out_0.4<-readRDS(paste0('sim_output/lgcp/disc_palm_output_',0.4,'_',round(rho_sd),'_',round(mu),'_',wsize,'.rds'))
full_out<-readRDS(paste0('sim_output/lgcp/full_output_',round(mu),'_',wsize,'.rds'))

### TEMPORARY REMOVE LATER
ind<-which(unlist(lapply(disc_palm_out_0.4,\(x) class(x)))=="character")
disc_palm_out_0.4[ind]<-NULL
######

### COMPUTATION TIMES IN MINUTES
## MEAN
palm_time<-round(apply(t(matrix(unlist(lapply(palm_out, \(x) c(x$init_time,x$cal2_time,x$init_time+x$cal2_time))),nrow=3))/60,2,mean),2)
disc_palm_time_0.2<-round(apply(t(matrix(unlist(lapply(disc_palm_out_0.2, \(x) c(x$init_time,x$cal2_time,x$init_time+x$cal2_time))),nrow=3))/60,2,mean),2)
disc_palm_time_0.4<-round(apply(t(matrix(unlist(lapply(disc_palm_out_0.4, \(x) c(x$init_time,x$cal2_time,x$init_time+x$cal2_time))),nrow=3))/60,2,mean),2)

round(apply(t(matrix(unlist(lapply(palm_out, \(x) x$palm_post_output[,6]/(x$init_time+x$cal2_time))),nrow=3)),2,mean),2)
round(apply(t(matrix(unlist(lapply(disc_palm_out_0.2, \(x) x$palm_post_output[,6]/(x$init_time+x$cal2_time))),nrow=3)),2,mean),2)
round(apply(t(matrix(unlist(lapply(disc_palm_out_0.4, \(x) x$palm_post_output[,6]/(x$init_time+x$cal2_time))),nrow=3)),2,mean),2)

## SD
round(apply(t(matrix(unlist(lapply(palm_out, \(x) c(x$init_time,x$cal2_time,x$init_time+x$cal2_time))),nrow=3))/60,2,sd),2)
round(apply(t(matrix(unlist(lapply(disc_palm_out_0.2, \(x) c(x$init_time,x$cal2_time,x$init_time+x$cal2_time))),nrow=3))/60,2,sd),2)
round(apply(t(matrix(unlist(lapply(disc_palm_out_0.4, \(x) c(x$init_time,x$cal2_time,x$init_time+x$cal2_time))),nrow=3))/60,2,sd),2)


# BIASES FROM POSTERIOR MEANS
apply(t(matrix(unlist(lapply(palm_out, \(x) x$cal2_post_output[,1]-c(mu,log(sig2),log(phi)))),nrow=3)),2,mean)
apply(t(matrix(unlist(lapply(disc_palm_out_0.2, \(x) x$cal2_post_output[,1]-c(mu,log(sig2),log(phi)))),nrow=3)),2,mean)
apply(t(matrix(unlist(lapply(disc_palm_out_0.4, \(x) x$cal2_post_output[,1]-c(mu,log(sig2),log(phi)))),nrow=3)),2,mean)

# RMSES FROM POSTERIOR MEANS
sqrt(apply(t(matrix(unlist(lapply(palm_out, \(x) (x$cal2_post_output[,1]-c(mu,log(sig2),log(phi)))^2)),nrow=3)),2,mean))
sqrt(apply(t(matrix(unlist(lapply(disc_palm_out_0.2, \(x) (x$cal2_post_output[,1]-c(mu,log(sig2),log(phi)))^2)),nrow=3)),2,mean))
sqrt(apply(t(matrix(unlist(lapply(disc_palm_out_0.4, \(x) (x$cal2_post_output[,1]-c(mu,log(sig2),log(phi)))^2)),nrow=3)),2,mean))

# COVERAGES
apply(t(matrix(unlist(lapply(palm_out, \(x) (x$palm_post_output[,4]<c(mu,log(sig2),log(phi)))*(x$palm_post_output[,5]>c(mu,log(sig2),log(phi))))),nrow=3)),2,mean)
apply(t(matrix(unlist(lapply(disc_palm_out_0.2, \(x) (x$palm_post_output[,4]<c(mu,log(sig2),log(phi)))*(x$palm_post_output[,5]>c(mu,log(sig2),log(phi))))),nrow=3)),2,mean)
apply(t(matrix(unlist(lapply(disc_palm_out_0.4, \(x) (x$palm_post_output[,4]<c(mu,log(sig2),log(phi)))*(x$palm_post_output[,5]>c(mu,log(sig2),log(phi))))),nrow=3)),2,mean)

apply(t(matrix(unlist(lapply(palm_out, \(x) (x$cal1_post_output[,4]<c(mu,log(sig2),log(phi)))*(x$cal1_post_output[,5]>c(mu,log(sig2),log(phi))))),nrow=3)),2,mean)
apply(t(matrix(unlist(lapply(disc_palm_out_0.2, \(x) (x$cal1_post_output[,4]<c(mu,log(sig2),log(phi)))*(x$cal1_post_output[,5]>c(mu,log(sig2),log(phi))))),nrow=3)),2,mean)
apply(t(matrix(unlist(lapply(disc_palm_out_0.4, \(x) (x$cal1_post_output[,4]<c(mu,log(sig2),log(phi)))*(x$cal1_post_output[,5]>c(mu,log(sig2),log(phi))))),nrow=3)),2,mean)

apply(t(matrix(unlist(lapply(palm_out, \(x) (x$cal2_post_output[,2]<c(mu,log(sig2),log(phi)))*(x$cal2_post_output[,3]>c(mu,log(sig2),log(phi))))),nrow=3)),2,mean)
apply(t(matrix(unlist(lapply(disc_palm_out_0.2, \(x) (x$cal2_post_output[,2]<c(mu,log(sig2),log(phi)))*(x$cal2_post_output[,3]>c(mu,log(sig2),log(phi))))),nrow=3)),2,mean)
apply(t(matrix(unlist(lapply(disc_palm_out_0.4, \(x) (x$cal2_post_output[,2]<c(mu,log(sig2),log(phi)))*(x$cal2_post_output[,3]>c(mu,log(sig2),log(phi))))),nrow=3)),2,mean)

# CI LENGTH
apply(t(matrix(unlist(lapply(disc_palm_out_0.2, \(x) x$cal1_post_output[,5]-x$cal1_post_output[,4])),nrow=3)),2,mean)
apply(t(matrix(unlist(lapply(disc_palm_out_0.4, \(x) x$cal2_post_output[,3]-x$cal2_post_output[,2])),nrow=3)),2,mean)



#### COMPARE FOR APPENDIX
apply(t(matrix(unlist(lapply(palm_out, \(x) x$cal2_post_output[,1]-c(mu,log(sig2),log(phi)))),nrow=3)),2,mean)
apply(t(matrix(unlist(lapply(disc_palm_out_0.2_250, \(x) x$cal2_post_output[,1]-c(mu,log(sig2),log(phi)))),nrow=3)),2,mean)
apply(t(matrix(unlist(lapply(disc_palm_out_0.2_750, \(x) x$cal2_post_output[,1]-c(mu,log(sig2),log(phi)))),nrow=3)),2,mean)

