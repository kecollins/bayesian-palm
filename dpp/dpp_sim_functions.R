simulate_dpp_data<-function(nsim,rho,lalpha,wsize){
  mG<-dppGauss(lambda=rho,alpha=exp(lalpha),d=2)
  S<-simulate.dppm(nsim=nsim,mG,W=owin(c(0,wsize),c(0,wsize)))
  saveRDS(S,paste0('sim_output/dpp/sim_data_',rho,'_',wsize,'.rds'))
}


disc_palm_dpp_sim_study<-function(nsim,rho,lalpha,wsize,R,nboot=100,empirical=TRUE,dr=0.1/500,iters=10000,burn=1000){
  S_all<-readRDS(paste0('sim_output/dpp/sim_data_',rho,'_',wsize,'.rds'))
  if(empirical){rho_sd<-1} else{rho_sd<-100}
  
  model<-stan_model('dpp/dpp_palm.stan')
  
  all_out<-foreach(k=c(1:nsim),.errorhandling = "pass")%do%{
    print(k)
    S<-S_all[[k]]
    
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
    
    init_time<-system.time({
      # RUN STAN MODEL
      post <- sampling(
        model,
        data = stan_data,
        chains = 1,
        warmup = burn,
        iter = iters,
        cores = 1,
        refresh = 0,
        init=list(list(rho=rho,lsig2=0,nu=10))
      )
    })
    
    post_mean<-get_posterior_mean(post,pars=c("rho","lalpha","alpha"))
    mG<-dppGauss(lambda=post_mean[1],alpha=exp(post_mean[2]),d=2)
    
    ### ADJUSTMENT 1
    cal1_time<-system.time({
      H_inv<-cov(cbind(rstan::extract(post,"rho")$rho,
                       rstan::extract(post,"lalpha")$lalpha))
      
      score_out<-foreach(k=c(1:1000),.combine='rbind')%dopar%{
        S_temp<-simulate.dppm(nsim=1,mG,W=owin(c(0,wsize),c(0,wsize)))
        if(S_temp$n==0){next}
        data_temp<-data_clean(S_temp,R,wsize)
        
        return(dpp_score_eval(data_temp,rho=post_mean[1],lalpha=post_mean[2],d_lamP=dpp_d_lamP,lamP=dpp_lamP))
      }
      
      J<-cov(score_out)
      
      eta<-1/sum(eigen(H_inv%*%J)$values)
      
      stan_data$eta<-eta
      
      cal1_post<-sampling(
        model,
        data = stan_data,
        chains = 1,
        warmup = burn,
        iter = iters,
        cores = 1,
        refresh = 0,
        init=list(list(rho=rho,lsig2=0,nu=10))
      )
    })
    
    ### ADJUSTEMENT 2
    cal2_time<-system.time({
      
      boot_out<-foreach(k=c(1:nboot))%dopar%{
        S_boot<-simulate.dppm(nsim=1,mG,W=owin(c(0,wsize),c(0,wsize)))
        data<-data_clean(S_boot,R=R,wsize=wsize)
        
        # SETUP STAN DATA
        stan_data_boot<-list(
          rho_mean=S_boot$n/area(owin(c(0,wsize),c(0,wsize))),
          rho_sd=rho_sd,
          dN=length(data$dG),
          dG=data$dG,
          dG_counts=data$dG_counts,
          dS_counts=data$dS_counts,
          dx=data$dx,
          eta=1
        )
        
        
        # RUN STAN MODEL
        post_boot <- sampling(
          model,  # Stan program compiled earlier
          data = stan_data_boot,    # named list of data
          chains = 1,             # number of Markov chains
          warmup = burn,          # number of warmup iterations per chain
          iter = iters,            # total number of iterations per chain
          cores = 1,              # number of cores (could use one per chain)
          refresh = 0,             # no progress shown
          init=list(list(rho=rho,lsig2=0,nu=10))
        )
        
        
        summ<-summary(post_boot,pars=c("rho","lalpha"),c(0.025,0.975))$summary
        return(summ[,c(1,4,5)])
      }
      
      # CALIBRATE RHO
      rho_mean<-unlist(lapply(boot_out,function(x) x[1,1]))
      rho_qs<-t(matrix(unlist(lapply(boot_out,function(x) x[1,2:3])),nrow=2))-rho_mean
      rho_eta<-find_eta(y_mean=rho_mean,y_qs=rho_qs,alpha=0.05,truth=post_mean[1])
      rho_post<-rstan::extract(post,pars="rho")$rho
      calibrated_rho<-mean(rho_post)+rho_eta*(rho_post-mean(rho_post))
      
      # CALIBRATE LALPHA
      lalpha_mean<-unlist(lapply(boot_out,function(x) x[2,1]))
      lalpha_qs<-t(matrix(unlist(lapply(boot_out,function(x) x[2,2:3])),nrow=2))-lalpha_mean
      lalpha_eta<-find_eta(y_mean=lalpha_mean,y_qs=lalpha_qs,alpha=0.05,truth=post_mean[2])
      lalpha_post<-rstan::extract(post,pars="lalpha")$lalpha
      calibrated_lalpha<-mean(lalpha_post)+lalpha_eta*(lalpha_post-mean(lalpha_post))
      
    })
    
    palm_post_output<-summary(post,pars=c("rho","lalpha"),c(0.025,0.975))$summary
    cal1_post_output<-summary(cal1_post,pars=c("rho","lalpha"),c(0.025,0.975))$summary
    cal2_post_output<-rbind(c(mean(calibrated_rho),quantile(calibrated_rho,c(0.025,0.975))),
                            c(mean(calibrated_lalpha),quantile(calibrated_lalpha,c(0.025,0.975))))
    
    
    out<-list(palm_post_output=palm_post_output,
              cal1_post_output=cal1_post_output,
              cal2_post_output=cal2_post_output,
              init_time=init_time[3],
              cal1_time=cal1_time[3],
              cal2_time=cal2_time[3])
    print(init_time[3]+cal1_time[3]+cal2_time[3])
    rm(post,cal1_post,score_out,data,calibrated_rho,calibrated_lalpha)
    return(out)
  }
  
  saveRDS(all_out,paste0('sim_output/dpp/dpp_palm_output_',rho_sd,'_',rho,'_',wsize,'.rds'))
  
}

