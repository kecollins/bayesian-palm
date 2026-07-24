simulate_thomas_data<-function(nsim,mu,nu,si2,wsize){
  S<-rThomas(nsim=nsim,kappa=mu,mu=nu,scale=sqrt(sig2),win = owin(c(0,wsize),c(0,wsize)))
  saveRDS(S,paste0('sim_output/thomas_process/sim_data_',mu*nu,'_',wsize,'.rds'))
}

disc_palm_thomas_sim_study<-function(nsim,mu,sig2,phi,wsize,R,nboot=100,empirical=TRUE,dr=0.1/500,iters=10000,burn=1000){
  S_all<-readRDS(paste0('sim_output/thomas_process/sim_data_',mu*nu,'_',wsize,'.rds'))
  if(empirical){rho_sd<-1} else{rho_sd<-100}
  
  model<-stan_model('thomas_process/thomas_palm.stan')
  
  all_out<-foreach(k=c(1:nsim))%do%{
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
    
    post_mean<-get_posterior_mean(post,pars=c("mu","nu","lnu","sig2","lsig2"))
    
    ### ADJUSTMENT 1
    cal1_time<-system.time({
      H_inv<-cov(cbind(rstan::extract(post,"mu")$mu,
                       rstan::extract(post,"nu")$nu,
                       rstan::extract(post,"lsig2")$lsig2))
      
      score_out<-foreach(k=c(1:1000),.combine='rbind')%dopar%{
        S_temp<-rThomas(kappa=post_mean[1],mu=exp(post_mean[3]),scale=sqrt(exp(post_mean[5])))
        data_temp<-data_clean(S_temp,R,wsize)
        
        return(thomas_score_eval(data_temp,rho=post_mean[1],lnu=post_mean[3],lsig2=post_mean[5],d_lamP=thomas_d_lamP,lamP=thomas_lamP))
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
      S_boot<-bootstrap_ppp(nsim=nboot,model="Thomas",pars=list(kappa=post_mean[1],mu=post_mean[2],sig2=post_mean[4]),
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
        
        
        summ<-summary(post_boot,pars=c("mu","nu","lsig2"),c(0.025,0.975))$summary
        return(summ[,c(1,4,5)])
      }
      
      # CALIBRATE MU
      mu_mean<-unlist(lapply(boot_out,function(x) x[1,1]))
      mu_qs<-t(matrix(unlist(lapply(boot_out,function(x) x[1,2:3])),nrow=2))-mu_mean
      mu_eta<-find_eta(y_mean=mu_mean,y_qs=mu_qs,alpha=0.05,truth=post_mean[1])
      mu_post<-rstan::extract(post,pars="mu")$mu
      calibrated_mu<-mean(mu_post)+mu_eta*(mu_post-mean(mu_post))
      
      # CALIBRATE NU
      nu_mean<-unlist(lapply(boot_out,function(x) x[1,1]))
      nu_qs<-t(matrix(unlist(lapply(boot_out,function(x) x[1,2:3])),nrow=2))-nu_mean
      nu_eta<-find_eta(y_mean=nu_mean,y_qs=nu_qs,alpha=0.05,truth=post_mean[1])
      nu_post<-rstan::extract(post,pars="nu")$nu
      calibrated_nu<-mean(nu_post)+nu_eta*(nu_post-mean(nu_post))
      
      # CALIBRATE LSIG2
      lsig2_mean<-unlist(lapply(boot_out,function(x) x[1,1]))
      lsig2_qs<-t(matrix(unlist(lapply(boot_out,function(x) x[1,2:3])),nrow=2))-lsig2_mean
      lsig2_eta<-find_eta(y_mean=lsig2_mean,y_qs=lsig2_qs,alpha=0.05,truth=post_mean[1])
      lsig2_post<-rstan::extract(post,pars="lsig2")$lsig2
      calibrated_lsig2<-mean(lsig2_post)+lsig2_eta*(lsig2_post-mean(lsig2_post))

    })
    
    palm_post_output<-summary(post,pars=c("mu","nu","lsig2"),c(0.025,0.975))$summary
    cal1_post_output<-summary(cal1_post,pars=c("mu","nu","lsig2"),c(0.025,0.975))$summary
    cal2_post_output<-rbind(c(mean(calibrated_mu),quantile(calibrated_mu,c(0.025,0.975))),
                            c(mean(calibrated_nu),quantile(calibrated_nu,c(0.025,0.975))),
                            c(mean(calibrated_lsig2),quantile(calibrated_lsig2,c(0.025,0.975))))
    
    
    out<-list(palm_post_output=palm_post_output,
              cal1_post_output=cal1_post_output,
              cal2_post_output=cal2_post_output,
              init_time=init_time[3],
              cal1_time=cal1_time[3],
              cal2_time=cal2_time[3])
    print(init_time[3]+cal1_time[3]+cal2_time[3])
    rm(post,cal1_post,score_out,data,calibrated_mu,calibrated_nu,calibrated_lsig2,S_boot)
    return(out)
  }
  
  saveRDS(all_out,paste0('sim_output/thomas-process/thomas_palm_output_',mu*nu,'_',wsize,'.rds'))
  
}



### ONE OFF MODEL FIT WRAPPER
fit_thomas_lgcp<-function(S,R,wsize,empirical=TRUE,nboot=100,iters=10000,burn=1000){
  model<-stan_model('lgcp/lgcp_palm_disc.stan')
  if(empirical){rho_sd<-1} else{rho_sd<-100}
  # ORGANIZE DATA
  data<-data_clean(S=S,R=R,wsize=wsize,dr=dr)
  
  cat("\n Sampling initial MCMC...")
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
      init=list(list(rho=rho,lsig2=0,lphi=-2.3))#,
      #control = list(adapt_delta = 0.99)
    )
  })
  cat(paste0(round(init_time[3],2)," seconds"))
  
  post_mean<-get_posterior_mean(post,pars=c("mu","sig2","phi","rho","lsig2","lphi"))
  
  cat("\n Computing first posterior adjustment...")
  ### ADJUSTMENT 1
  cal1_time<-system.time({
    H_inv<-cov(cbind(rstan::extract(post,"mu")$mu,
                     rstan::extract(post,"lsig2")$lsig2,
                     rstan::extract(post,"lphi")$lphi))
    
    score_out<-foreach(k=c(1:1000),.combine='rbind')%dopar%{
      S_temp<-rLGCP(model="exponential",mu=log(post_mean[4])-exp(post_mean[5])/2,param=list(var=exp(post_mean[5]),scale=exp(post_mean[6])),
                    win=S$window)
      data_temp<-data_clean(S_temp,R,wsize)
      
      return(disc_score_eval(data_temp,rho=post_mean[4],lsig2=post_mean[5],lphi=post_mean[6],d_lamP = lgcp_d_lamP,lamP=lgcp_lamP))
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
      init=list(list(rho=rho,lsig2=0,lphi=-2.3))#,
      #control = list(adapt_delta = 0.99)
    )
  })
  cat(paste0(round(cal1_time[3],2)," seconds"))
  
  cat("\n Computing second posterior adjustment...")
  ### ADJUSTEMENT 2
  cal2_time<-system.time({
    S_boot<-bootstrap_ppp(nsim=nboot,model="LGCP",pars=list(mu=post_mean[1],pars=list(var=post_mean[2],scale=post_mean[3])),
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
  cat(paste0(round(cal2_time[3],2)," seconds"))
  
  return(list(palm_post=matrix(unlist(extract(post,pars=c("mu","lsig2","lphi"))),nrow=9000),
              cal1_post=matrix(unlist(extract(cal1_post,pars=c("mu","lsig2","lphi"))),nrow=9000),
              cal2_post=cbind(calibrated_mu,calibrated_lsig2,calibrated_lphi)))
}

