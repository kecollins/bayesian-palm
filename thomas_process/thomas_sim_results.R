

# SIMULATION RESULTS
nsim<-100
mu<-30
nu<-10
sig2<-0.0025
rho<-mu*nu
wsize<-1
S_all<-readRDS(paste0('sim_output/thomas_process/sim_data_',rho,'_',wsize,'.rds'))
disc_palm_out_0.2<-readRDS(paste0('sim_output/thomas_process/thomas_palm_output_',1,'_',rho,'_',wsize,'.rds'))
disc_palm_out_0.2_NE<-readRDS(paste0('sim_output/thomas_process/thomas_palm_output_',100,'_',rho,'_',wsize,'.rds'))


bias_0.2<-round(apply(t(matrix(unlist(lapply(disc_palm_out_0.2,\(x) x$cal2_post_output[,1]-c(log(mu),log(nu),log(sig2)))),nrow=3)),2,mean),2)
bias_0.2_NE<-round(apply(t(matrix(unlist(lapply(disc_palm_out_0.2_NE,\(x) x$cal2_post_output[,1]-c(log(mu),log(nu),log(sig2)))),nrow=3)),2,mean),2)

rmse_0.2<-round(sqrt(apply(t(matrix(unlist(lapply(disc_palm_out_0.2,\(x) (x$cal2_post_output[,1]-c(log(mu),log(nu),log(sig2)))^2)),nrow=3)),2,mean)),2)
rmse_0.2_NE<-round(sqrt(apply(t(matrix(unlist(lapply(disc_palm_out_0.2_NE,\(x) (x$cal2_post_output[,1]-c(log(mu),log(nu),log(sig2)))^2)),nrow=3)),2,mean)),2)

cov_0.2<-round(apply(t(matrix(unlist(lapply(disc_palm_out_0.2,\(x) (x$palm_post_output[,4]<c(log(mu),log(nu),log(sig2)))*(x$palm_post_output[,5]>c(log(mu),log(nu),log(sig2))))),nrow=3)),2,mean),2)
cov_0.2_NE<-round(apply(t(matrix(unlist(lapply(disc_palm_out_0.2_NE,\(x) (x$palm_post_output[,4]<c(log(mu),log(nu),log(sig2)))*(x$palm_post_output[,5]>c(log(mu),log(nu),log(sig2))))),nrow=3)),2,mean),2)

cov_0.2_1<-round(apply(t(matrix(unlist(lapply(disc_palm_out_0.2,\(x) (x$cal1_post_output[,4]<c(log(mu),log(nu),log(sig2)))*(x$cal1_post_output[,5]>c(log(mu),log(nu),log(sig2))))),nrow=3)),2,mean),2)
cov_0.2_1_NE<-round(apply(t(matrix(unlist(lapply(disc_palm_out_0.2_NE,\(x) (x$cal1_post_output[,4]<c(log(mu),log(nu),log(sig2)))*(x$cal1_post_output[,5]>c(log(mu),log(nu),log(sig2))))),nrow=3)),2,mean),2)

cov_0.2_2<-round(apply(t(matrix(unlist(lapply(disc_palm_out_0.2,\(x) (x$cal2_post_output[,2]<c(log(mu),log(nu),log(sig2)))*(x$cal2_post_output[,3]>c(log(mu),log(nu),log(sig2))))),nrow=3)),2,mean),2)
cov_0.2_2_NE<-round(apply(t(matrix(unlist(lapply(disc_palm_out_0.2_NE,\(x) (x$cal2_post_output[,2]<c(log(mu),log(nu),log(sig2)))*(x$cal2_post_output[,3]>c(log(mu),log(nu),log(sig2))))),nrow=3)),2,mean),2)

time_0.2<-round(apply(t(matrix(unlist(lapply(disc_palm_out_0.2,\(x) c(x$init_time,x$cal1_time,x$cal2_time))),nrow=3)),2,mean)/60,2)
time_0.2_NE<-round(apply(t(matrix(unlist(lapply(disc_palm_out_0.2_NE,\(x) c(x$init_time,x$cal1_time,x$cal2_time))),nrow=3)),2,mean)/60,2)

c(bias_0.2[1],rmse_0.2[1],cov_0.2[1],cov_0.2_1[1],cov_0.2_2[1])
bias_0.2_NE
