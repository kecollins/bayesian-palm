simulate_thomas_data<-function(nsim,mu,nu,si2,wsize){
  S<-rThomas(nsim=nsim,kappa=mu,mu=nu,scale=sqrt(sig2),win = owin(c(0,wsize),c(0,wsize)))
  saveRDS(S,paste0('sim_output/thomas_process/sim_data_',mu*nu,'_',wsize,'.rds'))
}

disc_palm_thomas_sim_study<-function(nsim,mu,sig2,phi,wsize,R,nboot=100,empirical=TRUE,dr=0.1/500,iters=10000,burn=1000){
  S_all<-readRDS(paste0('sim_output/thomas_process/sim_data_',mu*nu,'_',wsize,'.rds'))
  if(empirical){rho_sd<-1} else{rho_sd<-100}
  
  model<-stan_model('thomas_process/thomas_palm.stan')
  
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
        model,  # Stan program compiled earlier
        data = stan_data,    # named list of data
        chains = 1,             # number of Markov chains
        warmup = burn,          # number of warmup iterations per chain
        iter = iters,            # total number of iterations per chain
        cores = 1,              # number of cores (could use one per chain)
        refresh = 0,             # no progress shown
        init=list(list(rho=rho,lsig2=0,nu=10))#,
        #control = list(adapt_delta = 0.99)
      )
    })
    
    post_mean<-get_posterior_mean(post,pars=c("mu","lmu","nu","lnu","sig2","lsig2"))
    
    ### ADJUSTMENT 1
    cal1_time<-system.time({
      H_inv<-cov(cbind(rstan::extract(post,"lmu")$lmu,
                       rstan::extract(post,"lnu")$lnu,
                       rstan::extract(post,"lsig2")$lsig2))
      
      score_out<-foreach(k=c(1:1000),.combine='rbind')%dopar%{
        S_temp<-rThomas(kappa=exp(post_mean[2]),mu=exp(post_mean[4]),scale=sqrt(exp(post_mean[6])))
        if(S_temp$n==0){next}
        data_temp<-data_clean(S_temp,R,wsize)
        
        return(thomas_score_eval(data_temp,lmu=post_mean[2],lnu=post_mean[4],lsig2=post_mean[6],d_lamP=thomas_d_lamP,lamP=thomas_lamP))
      }
      
      J<-cov(score_out)
      
      eta<-1/sum(eigen(H_inv%*%J)$values)
      
      stan_data$eta<-eta
      
      cal1_post<-sampling(
        model,  # Stan program compiled earlier
        data = stan_data,    # named list of data
        chains = 1,             # number of Markov chains
        warmup = burn,          # number of warmup iterations per chain
        iter = iters,            # total number of iterations per chain
        cores = 1,              # number of cores (could use one per chain)
        refresh = 0,             # no progress shown
        init=list(list(rho=rho,lsig2=0,nu=10))#,
        #control = list(adapt_delta = 0.99)
      )
    })
    
    ### ADJUSTEMENT 2
    cal2_time<-system.time({
      S_boot<-bootstrap_ppp(nsim=nboot,model="Thomas",pars=list(kappa=exp(post_mean[2]),mu=exp(post_mean[4]),sig2=exp(post_mean[6])),
                            win=S$window)
      
      boot_out<-foreach(k=c(1:nboot))%dopar%{
        data<-data_clean(S_boot[[k]],R=R,wsize=wsize)
        
        # SETUP STAN DATA
        stan_data_boot<-list(
          rho_mean=S_boot[[k]]$n/area(owin(c(0,wsize),c(0,wsize))),
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
        
        
        summ<-summary(post_boot,pars=c("lmu","lnu","lsig2"),c(0.025,0.975))$summary
        return(summ[,c(1,4,5)])
      }
      
      # CALIBRATE LMU
      lmu_mean<-unlist(lapply(boot_out,function(x) x[1,1]))
      lmu_qs<-t(matrix(unlist(lapply(boot_out,function(x) x[1,2:3])),nrow=2))-lmu_mean
      lmu_eta<-find_eta(y_mean=lmu_mean,y_qs=lmu_qs,alpha=0.05,truth=post_mean[2])
      lmu_post<-rstan::extract(post,pars="lmu")$lmu
      calibrated_lmu<-mean(lmu_post)+lmu_eta*(lmu_post-mean(lmu_post))
      
      # CALIBRATE LNU
      lnu_mean<-unlist(lapply(boot_out,function(x) x[2,1]))
      lnu_qs<-t(matrix(unlist(lapply(boot_out,function(x) x[2,2:3])),nrow=2))-lnu_mean
      lnu_eta<-find_eta(y_mean=lnu_mean,y_qs=lnu_qs,alpha=0.05,truth=post_mean[4])
      lnu_post<-rstan::extract(post,pars="lnu")$lnu
      calibrated_lnu<-mean(lnu_post)+lnu_eta*(lnu_post-mean(lnu_post))
      
      # CALIBRATE LSIG2
      lsig2_mean<-unlist(lapply(boot_out,function(x) x[3,1]))
      lsig2_qs<-t(matrix(unlist(lapply(boot_out,function(x) x[3,2:3])),nrow=2))-lsig2_mean
      lsig2_eta<-find_eta(y_mean=lsig2_mean,y_qs=lsig2_qs,alpha=0.05,truth=post_mean[6])
      lsig2_post<-rstan::extract(post,pars="lsig2")$lsig2
      calibrated_lsig2<-mean(lsig2_post)+lsig2_eta*(lsig2_post-mean(lsig2_post))

    })
    
    palm_post_output<-summary(post,pars=c("lmu","lnu","lsig2"),c(0.025,0.975))$summary
    cal1_post_output<-summary(cal1_post,pars=c("lmu","lnu","lsig2"),c(0.025,0.975))$summary
    cal2_post_output<-rbind(c(mean(calibrated_lmu),quantile(calibrated_lmu,c(0.025,0.975))),
                            c(mean(calibrated_lnu),quantile(calibrated_lnu,c(0.025,0.975))),
                            c(mean(calibrated_lsig2),quantile(calibrated_lsig2,c(0.025,0.975))))
    
    
    out<-list(palm_post_output=palm_post_output,
              cal1_post_output=cal1_post_output,
              cal2_post_output=cal2_post_output,
              init_time=init_time[3],
              cal1_time=cal1_time[3],
              cal2_time=cal2_time[3])
    print(init_time[3]+cal1_time[3]+cal2_time[3])
    rm(post,cal1_post,score_out,data,calibrated_lmu,calibrated_lnu,calibrated_lsig2,S_boot)
    return(out)
  }
  
  saveRDS(all_out,paste0('sim_output/thomas_process/thomas_palm_output_',rho_sd,'_',mu*nu,'_',wsize,'.rds'))
  
}


