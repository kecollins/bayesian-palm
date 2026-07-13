

### SIMULATE B BOOTSTRAPPED DATASETS
bootstrap_ppp<-function(nsim,model,pars=NULL,dpp_object=NULL,win){
  if(model=="LGCP"){
    S_boot<-rLGCP(model="exponential",mu=pars$mu,param=pars$pars,win=win,nsim=nsim)
  } else if(model=="Thomas"){
    S_boot<-rThomas(kappa=pars$kappa,scale=pars$scale,mu=pars$mu,win=win,nsim=nsim)
  } else if(model=="dppGauss"){
    simulate.dppm(object=dpp_object,W=win,nsim=nsim)
  }
  
  return(S_boot)
}

### ORGANIZE DATA FOR MODEL
data_clean<-function(S,R,wsize,dr=R/1000,nx=20){
  
  # GET POINTS WITHIN R OF BORDER OF SQUARE REGION
  square<-st_cast(st_polygon(list((matrix(c(0,0, wsize,0, wsize,wsize, 0,wsize, 0,0), ncol=2, byrow=TRUE)))),"LINESTRING")
  S_df <- data.frame(S)
  S_sf <- st_as_sf(S_df, coords = c("x", "y"))
  window_keep<-lengths(st_is_within_distance(S_sf,square,dist=R))==0 # edge correction
  
  # GET PAIRWISE DISTANCE OF POINTS
  S_dist1<-as.vector(as.matrix(rdist(S_df[window_keep,],S_df[window_keep,])))
  S_dist2<-as.vector(as.matrix(rdist(S_df[window_keep,],S_df[!window_keep,])))
  temp<-S_dist1[which((S_dist1<R)&(S_dist1>0))]
  S_dist1R<-unique(temp) # unique because palm intensity is symmetric
  S_dist2R<-S_dist2[which((S_dist2<R)&(S_dist2>0))]
  
  # BIN PAIRWISE DISTANCES
  dG<-seq(dr/2,R-dr/2,by=dr)
  
  S_dist<-rdist(S_df,S_df)
  S_dist<-S_dist[upper.tri(S_dist)]
  S_distR<-S_dist[which(S_dist<R)]
  
  G_dist<-as.vector(rdist(S_df,grid))
  G_distR<-G_dist[which(G_dist<R)]

  dS_counts<-grid_counts_1d(dG,S_distR)
  dG_counts<-grid_counts_1d(dG,G_distR)
  
  return(list(window_keep=window_keep,S_dist1R=S_dist1R,S_dist2R=S_dist2R,R=R,dr=dr,dG=dG,
              dG_counts=dG_counts,dS_counts=dS_counts,S_distR=S_distR,G_distR=G_distR,dx=1/nx^2))
}

### BISECTION SEARCH TO FIND ETA
find_eta<-function(y_mean,y_qs,init_interval=c(0,5),alpha,truth){
  interval<-init_interval
  stop<-FALSE
  while(stop!=TRUE){
    eta<-mean(interval)
    eta_qs<-y_mean+eta*y_qs
    cov<-mean((eta_qs[,1]<truth)*(eta_qs[,2]>truth))
    if(cov<(1-alpha-0.02)){
      interval[1]<-eta
    } else if(cov>(1-alpha+0.02)){
      interval[2]<-eta
    } else{
      stop<-TRUE
    }
    
    if(abs(interval[1]-interval[2])<1e-5){interval<-init_interval+c(0,100)} # fail safe if initial upper bound is too small
  }
  
  return(eta)
}

### COUNT POINTS IN GRID CELL
grid_counts<-function(grid,x){
  cells<-nn2(grid,x,k=1)$nn.id
  cells<-data.frame(v1=cells)%>%group_by(v1)%>%summarize(n=n())
  cells<-rbind(cells,data.frame(v1=1:nrow(grid),n=0))
  return(as.vector((cells%>%group_by(v1)%>%summarize(n=sum(n)))[,2])$n)
}

grid_counts_1d<-function(grid,x){
  cells<-nn2(grid,x,k=1)$nn.id
  cells<-data.frame(v1=cells)%>%group_by(v1)%>%summarize(n=n())
  cells<-rbind(cells,data.frame(v1=1:length(grid),n=0))
  return(as.vector((cells%>%group_by(v1)%>%summarize(n=sum(n)))[,2])$n)
}

### CALCULATING SCORES
lgcp_lamP<-function(d,rho,lsig2,lphi){
  rho*exp(exp(lsig2)*exp(-d/exp(lphi)))
}

lgcp_d_lamP<-function(d,rho,lsig2,lphi){
  d1<-exp(exp(lsig2)*exp(-d/exp(lphi)))
  d2<-rho*exp(lsig2)*exp(-d/exp(lphi))*exp(exp(lsig2)*exp(-d/exp(lphi)))
  d3<-rho*d*exp(exp(lsig2-d*exp(-lphi))+lsig2-lphi-d*exp(-lphi))
  
  return(cbind(d1,d2,d3))
}

score_eval<-function(data,rho,lsig2,lphi,d_lamP,lamP){
  r<-seq(0,data$R,by=data$dr)
  apply(2*d_lamP(data$S_dist1R,rho,lsig2,lphi)/lamP(data$S_dist1R,rho,lsig2,lphi),2,sum)+
    apply(d_lamP(data$S_dist2R,rho,lsig2,lphi)/lamP(data$S_dist2R,rho,lsig2,lphi),2,sum)-
    sum(data$window_keep)*2*pi*apply(d_lamP(r,rho,lsig2,lphi)*r*data$dr,2,sum)
}



