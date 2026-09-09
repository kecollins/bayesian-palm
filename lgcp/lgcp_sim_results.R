library(kableExtra)

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

# SIMULATION (2) RESULTS
wsize<-1
mu<-log(1200)-1/2
S_all_sim2<-readRDS(paste0('sim_output/lgcp/sim_data_',round(mu),'_',wsize,'.rds'))
disc_palm_out_0.2_sim2<-readRDS(paste0('sim_output/lgcp/disc_palm_output_',0.2,'_',500,'_',1,'_',round(mu),'_',wsize,'.rds'))
disc_palm_out_0.4_sim2<-readRDS(paste0('sim_output/lgcp/disc_palm_output_',0.4,'_',500,'_',1,'_',round(mu),'_',wsize,'.rds'))
full_out_sim2<-readRDS(paste0('sim_output/lgcp/full_output_timing_',round(mu),'_',wsize,'.rds'))

# SIMULATION (3) RESULTS
wsize<-2
mu<-log(300)-1/2
S_all_sim3<-readRDS(paste0('sim_output/lgcp/sim_data_',round(mu),'_',wsize,'.rds'))
disc_palm_out_0.2_sim3<-readRDS(paste0('sim_output/lgcp/disc_palm_output_',0.2,'_',500,'_',1,'_',round(mu),'_',wsize,'.rds'))
disc_palm_out_0.4_sim3<-readRDS(paste0('sim_output/lgcp/disc_palm_output_',0.4,'_',500,'_',1,'_',round(mu),'_',wsize,'.rds'))
full_out_sim3<-readRDS(paste0('sim_output/lgcp/full_output_timing_',round(mu),'_',wsize,'.rds'))


# ------------------------------------------------------------------------------#
### COMPUTATION TIMES IN MINUTES
# RUNTIMES
full_time_sim1<-round(mean(unlist(lapply(full_out[1:10], \(x) x$time/60))),2)
disc_palm_time_0.2_sim1<-round(apply(t(matrix(unlist(lapply(disc_palm_out_0.2[1:10], \(x) c(x$init_time,x$init_time+x$cal1_time,x$init_time+x$cal2_time))),nrow=3))/60,2,mean),2)
disc_palm_time_0.4_sim1<-round(apply(t(matrix(unlist(lapply(disc_palm_out_0.4[1:10], \(x) c(x$init_time,x$init_time+x$cal1_time,x$init_time+x$cal2_time))),nrow=3))/60,2,mean),2)

full_time_sim2<-round(mean(unlist(lapply(full_out_sim2[1:10], \(x) x$time/60))),2)
disc_palm_time_0.2_sim2<-round(apply(t(matrix(unlist(lapply(disc_palm_out_0.2_sim2[1:10], \(x) c(x$init_time,x$init_time+x$cal1_time,x$init_time+x$cal2_time))),nrow=3))/60,2,mean),2)
disc_palm_time_0.4_sim2<-round(apply(t(matrix(unlist(lapply(disc_palm_out_0.4_sim2[1:10], \(x) c(x$init_time,x$init_time+x$cal1_time,x$init_time+x$cal2_time))),nrow=3))/60,2,mean),2)

full_time_sim3<-round(mean(unlist(lapply(full_out_sim3[1:10], \(x) x$time/60*100))),2)
disc_palm_time_0.2_sim3<-round(apply(t(matrix(unlist(lapply(disc_palm_out_0.2_sim3[1:10], \(x) c(x$init_time,x$init_time+x$cal1_time,x$init_time+x$cal2_time))),nrow=3))/60,2,mean),2)
disc_palm_time_0.4_sim3<-round(apply(t(matrix(unlist(lapply(disc_palm_out_0.4_sim3[1:10], \(x) c(x$init_time,x$init_time+x$cal1_time,x$init_time+x$cal2_time))),nrow=3))/60,2,mean),2)


# ESS
full_ESS_sim1<-round(mean(unlist(lapply(full_out[1:10], \(x) x$full_post_output[2,6]/x$time))),2)
disc_0.2_ESS_sim1<-round(apply(t(matrix(unlist(lapply(disc_palm_out_0.2[1:10], \(x) c(x$palm_post_output[2,6]/x$init_time,
                                                                                      x$cal1_post_output[2,6]/(x$init_time+x$cal1_time),
                                                                                      x$palm_post_output[2,6]/(x$init_time+x$cal2_time)))),nrow=3)),2,mean),2)
disc_0.4_ESS_sim1<-round(apply(t(matrix(unlist(lapply(disc_palm_out_0.4[1:10], \(x) c(x$palm_post_output[2,6]/x$init_time,
                                                                                      x$cal1_post_output[2,6]/(x$init_time+x$cal1_time),
                                                                                      x$palm_post_output[2,6]/(x$init_time+x$cal2_time)))),nrow=3)),2,mean),2)

full_ESS_sim2<-round(mean(unlist(lapply(full_out_sim2[1:10], \(x) x$full_post_output[2,6]/x$time))),2)
disc_0.2_ESS_sim2<-round(apply(t(matrix(unlist(lapply(disc_palm_out_0.2_sim2[1:10], \(x) c(x$palm_post_output[2,6]/x$init_time,
                                                                                           x$cal1_post_output[2,6]/(x$init_time+x$cal1_time),
                                                                                           x$palm_post_output[2,6]/(x$init_time+x$cal2_time)))),nrow=3)),2,mean),2)
disc_0.4_ESS_sim2<-round(apply(t(matrix(unlist(lapply(disc_palm_out_0.4_sim2[1:10], \(x) c(x$palm_post_output[2,6]/x$init_time,
                                                                                           x$cal1_post_output[2,6]/(x$init_time+x$cal1_time),
                                                                                           x$palm_post_output[2,6]/(x$init_time+x$cal2_time)))),nrow=3)),2,mean),2)

full_ESS_sim3<-round(mean(unlist(lapply(full_out_sim3[1:10], \(x) x$full_post_output[2,6]/x$time))),2)
disc_0.2_ESS_sim3<-round(apply(t(matrix(unlist(lapply(disc_palm_out_0.2_sim3[1:10], \(x) c(x$palm_post_output[2,6]/x$init_time,
                                                                                           x$cal1_post_output[2,6]/(x$init_time+x$cal1_time),
                                                                                           x$palm_post_output[2,6]/(x$init_time+x$cal2_time)))),nrow=3)),2,mean),2)
disc_0.4_ESS_sim3<-round(apply(t(matrix(unlist(lapply(disc_palm_out_0.4_sim3[1:10], \(x) c(x$palm_post_output[2,6]/x$init_time,
                                                                                           x$cal1_post_output[2,6]/(x$init_time+x$cal1_time),
                                                                                           x$palm_post_output[2,6]/(x$init_time+x$cal2_time)))),nrow=3)),2,mean),2)


## ADJUST FULL TIME SIM 3 TO ACCOUNT FOR ONLY 10 WARMUP, 100 ITERATIONS
full_ESS_sim3<-NA

time_df<-data.frame(rbind(c(full_time_sim1,NA,NA,full_ESS_sim1,NA,NA),
                          c(disc_palm_time_0.2_sim1,disc_0.2_ESS_sim1),
                          c(disc_palm_time_0.4_sim1,disc_0.4_ESS_sim1),
                          c(full_time_sim2,NA,NA,full_ESS_sim2,NA,NA),
                          c(disc_palm_time_0.2_sim2,disc_0.2_ESS_sim2),
                          c(disc_palm_time_0.4_sim2,disc_0.4_ESS_sim2),
                          c(full_time_sim3,NA,NA,full_ESS_sim3,NA,NA),
                          c(disc_palm_time_0.2_sim3,disc_0.2_ESS_sim3),
                          c(disc_palm_time_0.4_sim3,disc_0.4_ESS_sim3)))

time_df<-cbind(c("FL","PLE_0.2","PLE_0.4","FL","PLE_0.2","PLE_0.4","FL","PLE_0.2","PLE_0.4"),time_df)
colnames(time_df)<-NULL
time_df
time_df %>% kable(format="latex")
# ------------------------------------------------------------------------------#



# ------------------------------------------------------------------------------#
# TRUE VALUES
wsize<-1
sig2<-1
phi<-0.1
mu<-log(300)-1/2

### POINT ESTIMATION
# BIASES FROM POSTERIOR MEANS
bias_full<-round(apply(t(matrix(unlist(lapply(full_out, \(x) x$full_post_output[,1]-c(mu,log(sig2),log(phi)))),nrow=3)),2,mean),2)
bias_0.2<-round(apply(t(matrix(unlist(lapply(disc_palm_out_0.2, \(x) x$cal2_post_output[,1]-c(mu,log(sig2),log(phi)))),nrow=3)),2,mean),2)
bias_0.2_NE<-round(apply(t(matrix(unlist(lapply(disc_palm_out_0.2_NE, \(x) x$cal2_post_output[,1]-c(mu,log(sig2),log(phi)))),nrow=3)),2,mean),2)
bias_0.4<-round(apply(t(matrix(unlist(lapply(disc_palm_out_0.4, \(x) x$cal2_post_output[,1]-c(mu,log(sig2),log(phi)))),nrow=3)),2,mean),2)
bias_0.4_NE<-round(apply(t(matrix(unlist(lapply(disc_palm_out_0.4_NE, \(x) x$cal2_post_output[,1]-c(mu,log(sig2),log(phi)))),nrow=3)),2,mean),2)


# RMSES FROM POSTERIOR MEANS
rmse_full<-round(sqrt(apply(t(matrix(unlist(lapply(full_out, \(x) (x$full_post_output[,1]-c(mu,log(sig2),log(phi)))^2)),nrow=3)),2,mean)),2)
rmse_0.2<-round(sqrt(apply(t(matrix(unlist(lapply(disc_palm_out_0.2, \(x) (x$cal2_post_output[,1]-c(mu,log(sig2),log(phi)))^2)),nrow=3)),2,mean)),2)
rmse_0.2_NE<-round(sqrt(apply(t(matrix(unlist(lapply(disc_palm_out_0.2_NE, \(x) (x$cal2_post_output[,1]-c(mu,log(sig2),log(phi)))^2)),nrow=3)),2,mean)),2)
rmse_0.4<-round(sqrt(apply(t(matrix(unlist(lapply(disc_palm_out_0.4, \(x) (x$cal2_post_output[,1]-c(mu,log(sig2),log(phi)))^2)),nrow=3)),2,mean)),2)
rmse_0.4_NE<-round(sqrt(apply(t(matrix(unlist(lapply(disc_palm_out_0.4_NE, \(x) (x$cal2_post_output[,1]-c(mu,log(sig2),log(phi)))^2)),nrow=3)),2,mean)),2)
# ------------------------------------------------------------------------------#

# ------------------------------------------------------------------------------#
### EMPIRICAL COVERAGE RATES
cov_full<-apply(t(matrix(unlist(lapply(full_out, \(x) (x$full_post_output[,4]<c(mu,log(sig2),log(phi)))*(x$full_post_output[,5]>c(mu,log(sig2),log(phi))))),nrow=3)),2,mean)
cov_0.2<-apply(t(matrix(unlist(lapply(disc_palm_out_0.2, \(x) (x$palm_post_output[,4]<c(mu,log(sig2),log(phi)))*(x$palm_post_output[,5]>c(mu,log(sig2),log(phi))))),nrow=3)),2,mean)
cov_0.2_NE<-apply(t(matrix(unlist(lapply(disc_palm_out_0.2_NE, \(x) (x$palm_post_output[,4]<c(mu,log(sig2),log(phi)))*(x$palm_post_output[,5]>c(mu,log(sig2),log(phi))))),nrow=3)),2,mean)
cov_0.4<-apply(t(matrix(unlist(lapply(disc_palm_out_0.4, \(x) (x$palm_post_output[,4]<c(mu,log(sig2),log(phi)))*(x$palm_post_output[,5]>c(mu,log(sig2),log(phi))))),nrow=3)),2,mean)
cov_0.4_NE<-apply(t(matrix(unlist(lapply(disc_palm_out_0.4_NE, \(x) (x$palm_post_output[,4]<c(mu,log(sig2),log(phi)))*(x$palm_post_output[,5]>c(mu,log(sig2),log(phi))))),nrow=3)),2,mean)

cov_0.2_1<-apply(t(matrix(unlist(lapply(disc_palm_out_0.2, \(x) (x$cal1_post_output[,4]<c(mu,log(sig2),log(phi)))*(x$cal1_post_output[,5]>c(mu,log(sig2),log(phi))))),nrow=3)),2,mean)
cov_0.2_1_NE<-apply(t(matrix(unlist(lapply(disc_palm_out_0.2_NE, \(x) (x$cal1_post_output[,4]<c(mu,log(sig2),log(phi)))*(x$cal1_post_output[,5]>c(mu,log(sig2),log(phi))))),nrow=3)),2,mean)
cov_0.4_1<-apply(t(matrix(unlist(lapply(disc_palm_out_0.4, \(x) (x$cal1_post_output[,4]<c(mu,log(sig2),log(phi)))*(x$cal1_post_output[,5]>c(mu,log(sig2),log(phi))))),nrow=3)),2,mean)
cov_0.4_1_NE<-apply(t(matrix(unlist(lapply(disc_palm_out_0.4_NE, \(x) (x$cal1_post_output[,4]<c(mu,log(sig2),log(phi)))*(x$cal1_post_output[,5]>c(mu,log(sig2),log(phi))))),nrow=3)),2,mean)

cov_0.2_2<-apply(t(matrix(unlist(lapply(disc_palm_out_0.2, \(x) (x$cal2_post_output[,2]<c(mu,log(sig2),log(phi)))*(x$cal2_post_output[,3]>c(mu,log(sig2),log(phi))))),nrow=3)),2,mean)
cov_0.2_2_NE<-apply(t(matrix(unlist(lapply(disc_palm_out_0.2_NE, \(x) (x$cal2_post_output[,2]<c(mu,log(sig2),log(phi)))*(x$cal2_post_output[,3]>c(mu,log(sig2),log(phi))))),nrow=3)),2,mean)
cov_0.4_2<-apply(t(matrix(unlist(lapply(disc_palm_out_0.4, \(x) (x$cal2_post_output[,2]<c(mu,log(sig2),log(phi)))*(x$cal2_post_output[,3]>c(mu,log(sig2),log(phi))))),nrow=3)),2,mean)
cov_0.4_2_NE<-apply(t(matrix(unlist(lapply(disc_palm_out_0.4_NE, \(x) (x$cal2_post_output[,2]<c(mu,log(sig2),log(phi)))*(x$cal2_post_output[,3]>c(mu,log(sig2),log(phi))))),nrow=3)),2,mean)

# CI LENGTH
length_full<-round(apply(t(matrix(unlist(lapply(full_out, \(x) (x$full_post_output[,5]-x$full_post_output[,4]))),nrow=3)),2,mean),2)
length_0.2<-round(apply(t(matrix(unlist(lapply(disc_palm_out_0.2, \(x) x$palm_post_output[,5]-x$palm_post_output[,4])),nrow=3)),2,mean),2)
length_0.2_NE<-round(apply(t(matrix(unlist(lapply(disc_palm_out_0.2_NE, \(x) x$palm_post_output[,5]-x$palm_post_output[,4])),nrow=3)),2,mean),2)
length_0.4<-round(apply(t(matrix(unlist(lapply(disc_palm_out_0.4, \(x) x$palm_post_output[,5]-x$palm_post_output[,4])),nrow=3)),2,mean),2)
length_0.4_NE<-round(apply(t(matrix(unlist(lapply(disc_palm_out_0.4_NE, \(x) x$palm_post_output[,5]-x$palm_post_output[,4])),nrow=3)),2,mean),2)

length_0.2_1<-round(apply(t(matrix(unlist(lapply(disc_palm_out_0.2, \(x) x$cal1_post_output[,5]-x$cal1_post_output[,4])),nrow=3)),2,mean),2)
length_0.2_1_NE<-round(apply(t(matrix(unlist(lapply(disc_palm_out_0.2_NE, \(x) x$cal1_post_output[,5]-x$cal1_post_output[,4])),nrow=3)),2,mean),2)
length_0.4_1<-round(apply(t(matrix(unlist(lapply(disc_palm_out_0.4, \(x) x$cal1_post_output[,5]-x$cal1_post_output[,4])),nrow=3)),2,mean),2)
length_0.4_1_NE<-round(apply(t(matrix(unlist(lapply(disc_palm_out_0.4_NE, \(x) x$cal1_post_output[,5]-x$cal1_post_output[,4])),nrow=3)),2,mean),2)

length_0.2_2<-round(apply(t(matrix(unlist(lapply(disc_palm_out_0.2, \(x) x$cal2_post_output[,3]-x$cal2_post_output[,2])),nrow=3)),2,mean),2)
length_0.2_2_NE<-round(apply(t(matrix(unlist(lapply(disc_palm_out_0.2_NE, \(x) x$cal2_post_output[,3]-x$cal2_post_output[,2])),nrow=3)),2,mean),2)
length_0.4_2<-round(apply(t(matrix(unlist(lapply(disc_palm_out_0.4, \(x) x$cal2_post_output[,3]-x$cal2_post_output[,2])),nrow=3)),2,mean),2)
length_0.4_2_NE<-round(apply(t(matrix(unlist(lapply(disc_palm_out_0.4_NE, \(x) x$cal2_post_output[,3]-x$cal2_post_output[,2])),nrow=3)),2,mean),2)
# ------------------------------------------------------------------------------#

# ------------------------------------------------------------------------------#
# INFERENCE TABLES

# MU
mu_df<-as.data.frame(rbind(c(bias_full[1],rmse_full[1],paste0(cov_full[1]," (",length_full[1],")"),"",""),
                           c(bias_0.2[1],rmse_0.2[1],paste0(cov_0.2[1]," (",length_0.2[1],")"),paste0(cov_0.2_1[1]," (",length_0.2_1[1],")"),paste0(cov_0.2_2[1]," (",length_0.2_2[1],")")),
                           c(bias_0.2_NE[1],rmse_0.2_NE[1],paste0(cov_0.2_NE[1]," (",length_0.2_NE[1],")"),paste0(cov_0.2_1_NE[1]," (",length_0.2_1_NE[1],")"),paste0(cov_0.2_2_NE[1]," (",length_0.2_2_NE[1],")")),
                           c(bias_0.4[1],rmse_0.4[1],paste0(cov_0.4[1]," (",length_0.4[1],")"),paste0(cov_0.4_1[1]," (",length_0.4_1[1],")"),paste0(cov_0.4_2[1]," (",length_0.4_2[1],")")),
                           c(bias_0.4_NE[1],rmse_0.4_NE[1],paste0(cov_0.4_NE[1]," (",length_0.4_NE[1],")"),paste0(cov_0.4_1_NE[1]," (",length_0.4_1_NE[1],")"),paste0(cov_0.4_2_NE[1]," (",length_0.4_2_NE[1],")"))))

# LOG(SIG2)
lsig2_df<-as.data.frame(rbind(c(bias_full[2],rmse_full[2],paste0(cov_full[2]," (",length_full[2],")"),"",""),
                           c(bias_0.2[2],rmse_0.2[2],paste0(cov_0.2[2]," (",length_0.2[2],")"),paste0(cov_0.2_1[2]," (",length_0.2_1[2],")"),paste0(cov_0.2_2[2]," (",length_0.2_2[2],")")),
                           c(bias_0.2_NE[2],rmse_0.2_NE[2],paste0(cov_0.2_NE[2]," (",length_0.2_NE[2],")"),paste0(cov_0.2_1_NE[2]," (",length_0.2_1_NE[2],")"),paste0(cov_0.2_2_NE[2]," (",length_0.2_2_NE[2],")")),
                           c(bias_0.4[2],rmse_0.4[2],paste0(cov_0.4[2]," (",length_0.4[2],")"),paste0(cov_0.4_1[2]," (",length_0.4_1[2],")"),paste0(cov_0.4_2[2]," (",length_0.4_2[2],")")),
                           c(bias_0.4_NE[2],rmse_0.4_NE[2],paste0(cov_0.4_NE[2]," (",length_0.4_NE[2],")"),paste0(cov_0.4_1_NE[2]," (",length_0.4_1_NE[2],")"),paste0(cov_0.4_2_NE[2]," (",length_0.4_2_NE[2],")"))))

# LOG(PHI)
lphi_df<-as.data.frame(rbind(c(bias_full[3],rmse_full[3],paste0(cov_full[3]," (",length_full[3],")"),"",""),
                              c(bias_0.2[3],rmse_0.2[3],paste0(cov_0.2[3]," (",length_0.2[3],")"),paste0(cov_0.2_1[3]," (",length_0.2_1[3],")"),paste0(cov_0.2_2[3]," (",length_0.2_2[3],")")),
                              c(bias_0.2_NE[3],rmse_0.2_NE[3],paste0(cov_0.2_NE[3]," (",length_0.2_NE[3],")"),paste0(cov_0.2_1_NE[3]," (",length_0.2_1_NE[3],")"),paste0(cov_0.2_2_NE[3]," (",length_0.2_2_NE[3],")")),
                              c(bias_0.4[3],rmse_0.4[3],paste0(cov_0.4[3]," (",length_0.4[3],")"),paste0(cov_0.4_1[3]," (",length_0.4_1[3],")"),paste0(cov_0.4_2[3]," (",length_0.4_2[3],")")),
                              c(bias_0.4_NE[3],rmse_0.4_NE[3],paste0(cov_0.4_NE[3]," (",length_0.4_NE[3],")"),paste0(cov_0.4_1_NE[3]," (",length_0.4_1_NE[3],")"),paste0(cov_0.4_2_NE[3]," (",length_0.4_2_NE[3],")"))))

point_est_df<-data.frame(cbind(mu_df[,1],lsig2_df[,1],lphi_df[,1],mu_df[,2],lsig2_df[,2],lphi_df[,2]))
rownames(point_est_df)<-c("FL","PLE_0.2","PL_0.2","PLE_0.4","PL_0.4")
colnames(point_est_df)<-NULL
kable(point_est_df,format = "latex",align="c")

cov_df<-data.frame(cbind(mu_df[,3:5],lsig2_df[,3:5],lphi_df[,3:5]))
rownames(cov_df)<-c("FL","PLE_0.2","PL_0.2","PLE_0.4","PL_0.4")
colnames(cov_df)<-NULL
kable(cov_df,format = "latex",align="c")


inference_df<-data.frame(cbind(mu_df,lsig2_df,lphi_df))

kable(inference_df,format = "latex",booktabs = TRUE,align="c")

# ------------------------------------------------------------------------------#







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



