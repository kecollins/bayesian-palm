
simulate_lgcp_data<-function(nsim,mu,sig2,phi,wsize){
  S<-rLGCP(nsim=nsim,model="exponential",mu=mu,param=list(var=sig2,scale=phi),win = owin(c(0,wsize),c(0,wsize)))
  saveRDS(S,paste0('sim_output/lgcp/sim_data_',round(mu),'_',wsize,'.rds'))
}



palm_lgcp_sim_study<-function(nsim,mu,sig2,phi,wsize,R,nboot=100,empirical=FALSE){
  S_all<-readRDS(paste0('sim_output/lgcp/sim_data_',round(mu),'_',wsize,'.rds'))
  if(empirical){rho_sd<-1} else{rho_sd<-30}
  
  all_out<-foreach(k=c(1:nsim))%do%{
    S<-S_all[[k]]
    
    # ORGANIZE DATA
    data<-data_clean(S=S,R=R,wsize)
    
    if(sum(data$window_keep)==0){return("border points")}
    # SETUP STAN DATA
    stan_data<-list(
      rho_mean=S$n/area(owin(c(0,wsize),c(0,wsize))),
      rho_sd=rho_sd,
      N11=length(data$S_dist1R),
      N12=length(data$S_dist2R),
      N2=data$R/data$dr,
      NR=sum(data$window_keep),
      dr=data$dr,
      d1=data$S_dist1R,
      d2=data$S_dist2R,
      dG=seq(data$dr/2,data$R-data$dr/2,by=data$dr)
    )
    
    
    init_time<-system.time({
      # RUN STAN MODEL
      post <- stan(
        file = "lgcp/lgcp_palm_EC.stan",  # Stan program
        data = stan_data,    # named list of data
        chains = 1,             # number of Markov chains
        warmup = 1000,          # number of warmup iterations per chain
        iter = 10000,            # total number of iterations per chain
        cores = 1,              # number of cores (could use one per chain)
        refresh = 0,             # no progress shown
        init=list(list(rho=rho,lsig2=0,lphi=-2.3))#,
        #control = list(adapt_delta = 0.99)
      )
    })
    
    post_mean<-get_posterior_mean(post,pars=c("mu","sig2","phi"))
    
    
    simulate_time<-system.time({
      S_boot<-bootstrap_ppp(nsim=nboot,model="LGCP",pars=list(mu=post_mean[1],pars=list(var=post_mean[2],scale=post_mean[3])),
                            win=S$window)
    })
    
    
    
    cal_time<-system.time({
      boot_out<-foreach(k=c(1:nboot))%dopar%{
        data<-data_clean(S_boot[[k]],R=R,wsize=wsize)
        
        # SETUP STAN DATA
        stan_data_boot<-list(
          rho_mean=S$n/area(owin(c(0,wsize),c(0,wsize))),
          rho_sd=rho_sd,
          N11=length(data$S_dist1R),
          N12=length(data$S_dist2R),
          N2=data$R/data$dr,
          NR=sum(data$window_keep),
          dr=data$dr,
          d1=data$S_dist1R,
          d2=data$S_dist2R,
          dG=seq(data$dr/2,data$R-data$dr/2,by=data$dr)
        )
        
        
        # RUN STAN MODEL
        post_boot <- stan(
          file = "lgcp/lgcp_palm_EC.stan",  # Stan program
          data = stan_data_boot,    # named list of data
          chains = 1,             # number of Markov chains
          warmup = 1000,          # number of warmup iterations per chain
          iter = 10000,            # total number of iterations per chain
          cores = 1,              # number of cores (could use one per chain)
          refresh = 0,             # no progress shown
          init=list(list(rho=rho,lsig2=0,lphi=-2.3))
        )
        
        
        summ<-summary(post_boot,pars=c("mu","lsig2","lphi"),c(0.025,0.975))$summary
        return(summ[,c(1,4,5)])
      }
      
      # CALIBRATE MU
      mu_mean<-unlist(lapply(boot_out,function(x) x[1,1]))
      mu_qs<-t(matrix(unlist(lapply(boot_out,function(x) x[1,2:3])),nrow=2))-mu_mean
      mu_eta<-find_eta(y_mean=mu_mean,y_qs=mu_qs,alpha=0.05,truth=post_mean[1])
      mu_post<-rstan::extract(post,pars="mu")$mu
      calibrated_mu<-mean(mu_post)+mu_eta*(mu_post-mean(mu_post))
      
      # CALIBRATE LSIG2
      lsig2_mean<-unlist(lapply(boot_out,function(x) x[1,1]))
      lsig2_qs<-t(matrix(unlist(lapply(boot_out,function(x) x[1,2:3])),nrow=2))-lsig2_mean
      lsig2_eta<-find_eta(y_mean=lsig2_mean,y_qs=lsig2_qs,alpha=0.05,truth=post_mean[1])
      lsig2_post<-rstan::extract(post,pars="lsig2")$lsig2
      calibrated_lsig2<-mean(lsig2_post)+lsig2_eta*(lsig2_post-mean(lsig2_post))
      
      # CALIBRATE LPHI
      lphi_mean<-unlist(lapply(boot_out,function(x) x[1,1]))
      lphi_qs<-t(matrix(unlist(lapply(boot_out,function(x) x[1,2:3])),nrow=2))-lphi_mean
      lphi_eta<-find_eta(y_mean=lphi_mean,y_qs=lphi_qs,alpha=0.05,truth=post_mean[1])
      lphi_post<-rstan::extract(post,pars="lphi")$lphi
      calibrated_lphi<-mean(lphi_post)+lphi_eta*(lphi_post-mean(lphi_post))
    })
    
    
    
    palm_post_output<-summary(post,pars=c("mu","lsig2","lphi"),c(0.025,0.975))$summary
    calibrated_post_output<-rbind(c(mean(calibrated_mu),quantile(calibrated_mu,c(0.025,0.975))),
                                  c(mean(calibrated_lsig2),quantile(calibrated_lsig2,c(0.025,0.975))),
                                  c(mean(calibrated_lphi),quantile(calibrated_lphi,c(0.025,0.975))))
    
    
    out<-list(palm_post_output=palm_post_output,
              calibrated_post_output=calibrated_post_output,
              init_time=init_time[3],
              simulate_time=simulate_time[3],
              cal_time=cal_time[3])
    
    rm(post,data,calibrated_mu,calibrated_lsig2,calibrated_lphi,S_boot)
    return(out)
  }
  
  saveRDS(all_out,paste0('sim_output/lgcp/sim_output_',R,'_',round(rho_sd),'_',round(mu),'_',wsize,'.rds'))
  
}



full_lgcp_sim_study<-function(nsim,mu,sig2,phi,wsize){
  S_all<-readRDS(paste0('sim_output/lgcp/sim_data_',round(mu),'_',wsize,'.rds'))
  
  all_out<-foreach(k=c(1:nsim))%dopar%{
    S<-S_all[[k]]
    
    # GRID
    nx <- 20*wsize
    grid<-as.matrix(expand.grid(seq(1/nx/2,wsize-1/nx/2,length.out=nx),seq(1/nx/2,wsize-1/nx/2,length.out=nx)))
    
    points<-as.matrix(as.data.frame(S))
    y<-grid_counts(grid,points)
    
    dist_mat<-rdist(grid)
    # SETUP STAN DATA
    stan_data<-list(
      N=nrow(grid),
      dx=1/nx^2,
      y=y,
      d=dist_mat
    )
    
    time<-system.time({
      # RUN STAN MODEL
      post_full <- stan(
        file = "lgcp/lgcp_full.stan",  # Stan program
        data = stan_data,    # named list of data
        chains = 1,             # number of Markov chains
        warmup = 1000,          # number of warmup iterations per chain
        iter = 10000,            # total number of iterations per chain
        cores = 1,              # number of cores (could use one per chain)
        refresh = 0,             # no progress shown
        init=list(list(rho=rho,lsig2=0,lphi=-2.3))
      )
    })
    
    out<-list(full_post_output=full_post_output)
    rm(post_full)
    return(out)
  }
}
