### LGCP
fit_palm_lgcp<-function(S,R,wsize,empirical=TRUE,nboot=100,iters=10000,burn=1000,dr=0.1/500){
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
      model,
      data = stan_data,
      chains = 1,
      warmup = burn,
      iter = iters,
      cores = 1,
      refresh = 0,
      init=list(list(rho=rho,lsig2=0,lphi=-2.3))
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
    
    score_out<-foreach(k=c(1:1000),.combine='rbind',.packages = "spatstat")%dopar%{
      S_temp<-rLGCP(model="exponential",mu=log(post_mean[4])-exp(post_mean[5])/2,param=list(var=exp(post_mean[5]),scale=exp(post_mean[6])))
      if(S_temp$n==0){
        return(c(NA,NA,NA))
      }
      
      data_temp<-data_clean(S_temp,R,wsize)
      
      return(disc_score_eval(data_temp,mu=post_mean[1],lsig2=post_mean[5],lphi=post_mean[6],d_lamP = lgcp_d_lamP,lamP=lgcp_lamP))
    }
    
    J<-cov(score_out,use="na.or.complete")
    
    eta<-1/sum(eigen(H_inv%*%J)$values)
    
    stan_data$eta<-eta
    
    suppressWarnings({
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
    lsig2_mean<-unlist(lapply(boot_out,function(x) x[2,1]))
    lsig2_qs<-t(matrix(unlist(lapply(boot_out,function(x) x[2,2:3])),nrow=2))-lsig2_mean
    lsig2_eta<-find_eta(y_mean=lsig2_mean,y_qs=lsig2_qs,alpha=0.05,truth=post_mean[5])
    lsig2_post<-rstan::extract(post,pars="lsig2")$lsig2
    calibrated_lsig2<-mean(lsig2_post)+lsig2_eta*(lsig2_post-mean(lsig2_post))
    
    # CALIBRATE LPHI
    lphi_mean<-unlist(lapply(boot_out,function(x) x[3,1]))
    lphi_qs<-t(matrix(unlist(lapply(boot_out,function(x) x[3,2:3])),nrow=2))-lphi_mean
    lphi_eta<-find_eta(y_mean=lphi_mean,y_qs=lphi_qs,alpha=0.05,truth=post_mean[6])
    lphi_post<-rstan::extract(post,pars="lphi")$lphi
    calibrated_lphi<-mean(lphi_post)+lphi_eta*(lphi_post-mean(lphi_post))
  })

  cat(paste0(round(cal2_time[3],2)," seconds"))
  
  return(list(palm_post=matrix(unlist(extract(post,pars=c("mu","lsig2","lphi"))),nrow=9000),
              cal1_post=matrix(unlist(extract(cal1_post,pars=c("mu","lsig2","lphi"))),nrow=9000),
              cal2_post=cbind(calibrated_mu,calibrated_lsig2,calibrated_lphi)))
}

### THOMAS PROCESS
fit_palm_thomas<-function(S,R,wsize,empirical=TRUE,nboot=100,iters=10000,burn=1000,dr=0.1/500){
  model<-stan_model('thomas_process/thomas_palm.stan')
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
      model,
      data = stan_data,
      chains = 1,
      warmup = burn,
      iter = iters,
      cores = 1,
      refresh = 0,
      init=list(list(rho=rho,lsig2=0,lphi=-2.3))
    )
  })
  cat(paste0(round(init_time[3],2)," seconds"))
  
  post_mean<-get_posterior_mean(post,pars=c("mu","lmu","nu","lnu","sig2","lsig2"))
  
  cat("\n Computing first posterior adjustment...")
  
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
    
    suppressWarnings({
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
  })
  cat(paste0(round(cal1_time[3],2)," seconds"))
  
  cat("\n Computing second posterior adjustment...")
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
  
  cat(paste0(round(cal2_time[3],2)," seconds"))
  
  return(list(palm_post=matrix(unlist(extract(post,pars=c("lmu","lnu","lsig2"))),nrow=9000),
              cal1_post=matrix(unlist(extract(cal1_post,pars=c("lmu","lnu","lsig2"))),nrow=9000),
              cal2_post=cbind(calibrated_lmu,calibrated_lnu,calibrated_lsig2)))
}

### THOMAS PROCESS
fit_palm_dpp<-function(S,R,wsize,empirical=TRUE,nboot=100,iters=10000,burn=1000,dr=0.1/500){
  model<-stan_model('dpp/dpp_palm.stan')
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
      model,
      data = stan_data,
      chains = 1,
      warmup = burn,
      iter = iters,
      cores = 1,
      refresh = 0,
      init=list(list(rho=rho,lsig2=0,lphi=-2.3))
    )
  })
  cat(paste0(round(init_time[3],2)," seconds"))
  
  post_mean<-get_posterior_mean(post,pars=c("mu","lmu","nu","lnu","sig2","lsig2"))
  
  cat("\n Computing first posterior adjustment...")
  
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
    
    suppressWarnings({
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
  })
  cat(paste0(round(cal1_time[3],2)," seconds"))
  
  cat("\n Computing second posterior adjustment...")
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
  
  cat(paste0(round(cal2_time[3],2)," seconds"))
  
  return(list(palm_post=matrix(unlist(extract(post,pars=c("lmu","lnu","lsig2"))),nrow=9000),
              cal1_post=matrix(unlist(extract(cal1_post,pars=c("lmu","lnu","lsig2"))),nrow=9000),
              cal2_post=cbind(calibrated_lmu,calibrated_lnu,calibrated_lsig2)))
}

