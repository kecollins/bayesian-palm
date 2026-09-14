

# SIMULATION RESULTS
nsim<-100
alpha<-0.05
rho<-75
wsize<-2
S_all<-readRDS(paste0('sim_output/dpp/sim_data_',rho,'_',wsize,'.rds'))
disc_palm_out_0.2<-readRDS(paste0('sim_output/dpp/dpp_palm_output_',1,'_',rho,'_',wsize,'.rds'))
disc_palm_out_0.2_NE<-readRDS(paste0('sim_output/dpp/dpp_palm_output_',100,'_',rho,'_',wsize,'.rds'))


bias_0.2<-round(apply(t(matrix(unlist(lapply(disc_palm_out_0.2,\(x) x$cal2_post_output[,1]-c(rho,log(alpha)))),nrow=2)),2,mean),2)
bias_0.2_NE<-round(apply(t(matrix(unlist(lapply(disc_palm_out_0.2_NE,\(x) x$cal2_post_output[,1]-c(rho,log(alpha)))),nrow=2)),2,mean),2)

rmse_0.2<-round(sqrt(apply(t(matrix(unlist(lapply(disc_palm_out_0.2,\(x) (x$cal2_post_output[,1]-c(rho,log(alpha)))^2)),nrow=2)),2,mean)),2)
rmse_0.2_NE<-round(sqrt(apply(t(matrix(unlist(lapply(disc_palm_out_0.2_NE,\(x) (x$cal2_post_output[,1]-c(rho,log(alpha)))^2)),nrow=2)),2,mean)),2)

cov_0.2<-round(apply(t(matrix(unlist(lapply(disc_palm_out_0.2,\(x) (x$palm_post_output[,4]<c(rho,log(alpha)))*(x$palm_post_output[,5]>c(rho,log(alpha))))),nrow=2)),2,mean),2)
cov_0.2_NE<-round(apply(t(matrix(unlist(lapply(disc_palm_out_0.2_NE,\(x) (x$palm_post_output[,4]<c(rho,log(alpha)))*(x$palm_post_output[,5]>c(rho,log(alpha))))),nrow=2)),2,mean),2)

cov_0.2_1<-round(apply(t(matrix(unlist(lapply(disc_palm_out_0.2,\(x) (x$cal1_post_output[,4]<c(rho,log(alpha)))*(x$cal1_post_output[,5]>c(rho,log(alpha))))),nrow=2)),2,mean),2)
cov_0.2_1_NE<-round(apply(t(matrix(unlist(lapply(disc_palm_out_0.2_NE,\(x) (x$cal1_post_output[,4]<c(rho,log(alpha)))*(x$cal1_post_output[,5]>c(rho,log(alpha))))),nrow=2)),2,mean),2)

cov_0.2_2<-round(apply(t(matrix(unlist(lapply(disc_palm_out_0.2,\(x) (x$cal2_post_output[,2]<c(rho,log(alpha)))*(x$cal2_post_output[,3]>c(rho,log(alpha))))),nrow=2)),2,mean),2)
cov_0.2_2_NE<-round(apply(t(matrix(unlist(lapply(disc_palm_out_0.2_NE,\(x) (x$cal2_post_output[,2]<c(rho,log(alpha)))*(x$cal2_post_output[,3]>c(rho,log(alpha))))),nrow=2)),2,mean),2)

length_0.2<-round(apply(t(matrix(unlist(lapply(disc_palm_out_0.2,\(x) x$palm_post_output[,5]-x$palm_post_output[,4])),nrow=2)),2,mean),2)
length_0.2_NE<-round(apply(t(matrix(unlist(lapply(disc_palm_out_0.2_NE,\(x) x$palm_post_output[,5]-x$palm_post_output[,4])),nrow=2)),2,mean),2)

length_0.2_1<-round(apply(t(matrix(unlist(lapply(disc_palm_out_0.2,\(x) x$cal1_post_output[,5]-x$cal1_post_output[,4])),nrow=2)),2,mean),2)
length_0.2_1_NE<-round(apply(t(matrix(unlist(lapply(disc_palm_out_0.2_NE,\(x) x$cal1_post_output[,5]-x$cal1_post_output[,4])),nrow=2)),2,mean),2)

length_0.2_2<-round(apply(t(matrix(unlist(lapply(disc_palm_out_0.2,\(x) x$cal2_post_output[,3]-x$cal2_post_output[,2])),nrow=2)),2,mean),2)
length_0.2_2_NE<-round(apply(t(matrix(unlist(lapply(disc_palm_out_0.2_NE,\(x) x$cal2_post_output[,3]-x$cal2_post_output[,2])),nrow=2)),2,mean),2)

time_0.2<-round(apply(t(matrix(unlist(lapply(disc_palm_out_0.2,\(x) c(x$init_time,x$init_time+x$cal1_time,x$init_time+x$cal2_time))),nrow=3)),2,mean)/60,2)
time_0.2_NE<-round(apply(t(matrix(unlist(lapply(disc_palm_out_0.2_NE,\(x) c(x$init_time,x$init_time+x$cal1_time,x$init_time+x$cal2_time))),nrow=3)),2,mean)/60,2)

ESS_0.2<-round(c(mean(unlist(lapply(disc_palm_out_0.2,\(x) x$palm_post_output[1,6]/x$init_time))),
                 mean(unlist(lapply(disc_palm_out_0.2,\(x) x$cal1_post_output[1,6]/(x$init_time+x$cal1_time)))),
                 mean(unlist(lapply(disc_palm_out_0.2,\(x) x$palm_post_output[1,6]/(x$init_time+x$cal2_time))))),2)
ESS_0.2_NE<-round(c(mean(unlist(lapply(disc_palm_out_0.2_NE,\(x) x$palm_post_output[1,6]/x$init_time))),
                    mean(unlist(lapply(disc_palm_out_0.2_NE,\(x) x$cal1_post_output[1,6]/(x$init_time+x$cal1_time)))),
                    mean(unlist(lapply(disc_palm_out_0.2_NE,\(x) x$palm_post_output[1,6]/(x$init_time+x$cal2_time))))),2)

c(bias_0.2[1],rmse_0.2[1],cov_0.2[1],cov_0.2_1[1],cov_0.2_2[1])
bias_0.2_NE

bias_df<-data.frame(rbind(c(bias_0.2,rmse_0.2),
                          c(bias_0.2_NE,rmse_0.2_NE)))

bias_df %>% kable(format="latex")

cov_df<-data.frame(rbind(as.vector(t(cbind(cov_0.2,cov_0.2_1,cov_0.2_2))),
                         as.vector(t(cbind(cov_0.2_NE,cov_0.2_1_NE,cov_0.2_2_NE)))))

cov_df<-data.frame(rbind(as.vector(t(cbind(paste0(cov_0.2, " (",length_0.2,")"),paste0(cov_0.2_1, " (",length_0.2_1,")"),
                                           paste0(cov_0.2_2, " (",length_0.2_2,")")))),
                         as.vector(t(cbind(paste0(cov_0.2_NE, " (",length_0.2_NE,")"),paste0(cov_0.2_1_NE, " (",length_0.2_1_NE,")"),
                                           paste0(cov_0.2_2_NE, " (",length_0.2_2_NE,")"))))))

cov_df %>% kable(format="latex")

time_df<-data.frame(rbind(c(time_0.2,ESS_0.2),
                          c(time_0.2_NE,ESS_0.2_NE)))

time_df %>% kable(format="latex")


# ------------------------------------------------------------------------------#
# PRINT OUT ALL TABLE DATA

# Values for DPP section of Table 1
print("Table 1: DPP")
print(time_df)

# Values for Table 6
print("Table 6:")
print(bias_df)

# Values for Table 7
print("Table 7:")
print(cov_df)
