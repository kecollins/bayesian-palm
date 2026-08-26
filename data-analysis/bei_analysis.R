source('helper_functions.R')
library(RANN)
library(spatstat)
library(fields)
library(tidyverse)
library(rstan)
library(doParallel)
registerDoParallel(cores=10)


## DESIGN MATRIX, CENTER AND SCALE
X<-cbind(1,as.vector(t(bei.extra$elev$v)),as.vector(t(bei.extra$grad$v)))
X[,2]<-(X[,2]-mean(X[,2]))/sd(X[,2])
X[,3]<-(X[,3]-mean(X[,3]))/sd(X[,3])


S_df<-cbind(bei$x,bei$y)
dx<-bei.extra$elev$xstep
x_grid<-seq(bei.extra$elev$xrange[1]+dx/2,bei.extra$elev$xrange[2]-dx/2,by=dx)
y_grid<-seq(bei.extra$elev$yrange[1]+dx/2,bei.extra$elev$yrange[2]-dx/2,by=dx)

grid<-as.matrix(expand.grid(x_grid,y_grid))

R<-200
dr<-200/500

real_data_clean<-function(S_df,X,grid,R,dr,dx){
  dG<-seq(dr/2,R-dr/2,by=dr)
  
  S_dist_mat<-as.matrix(rdist(S_df))
  S_distR<-S_dist_mat[upper.tri(S_dist_mat)]
  S_distR<-S_distR[which(S_distR<R)]
  
  dS_counts<-grid_counts_1d(dG,S_distR)
  rm(S_distR)
  
  dS_ind<-which((S_dist_mat<R)&(S_dist_mat>0),arr.ind=TRUE)
  rm(S_dist_mat)
  
  #### IDENTIFY NEAREST GRID CELL TO EACH POINT
  cells<-nn2(grid,S_df,k=1)$nn.id
  cell_pairs<-data.frame(cbind(cells,as.vector(table(dS_ind[,1]))))
  rm(cells,dS_ind)
  cell_pairs<-rbind(cell_pairs,data.frame(X1=1:nrow(grid),X2=0))
  X_counts<-(cell_pairs %>% group_by(X1) %>% summarize(ct=sum(X2)))$ct
  rm(cell_pairs)
  
  t_X<-t(X)
  
  ## GET BIN COUNTS FOR COMPENSATOR
  G_dist_mat<-as.matrix(rdist(S_df,grid))
  keep<-G_dist_mat<R
  temp<-cbind(which(keep,arr.ind=TRUE)[,2],G_dist_mat[keep])
  rm(G_dist_mat,keep)
  bins<-factor(floor((temp[,2]-0)/dr)+1,1:length(dG))
  mat_counts<-matrix(table(temp[,1],bins),ncol=length(dG))
  rm(temp,bins)
  
  list(
    dN = length(dG),
    p = nrow(t_X),
    xN = ncol(t_X),
    dG = dG,
    dS_counts = dS_counts,
    t_X = t_X,
    X_counts = X_counts,
    mat_counts=mat_counts,
    dx = dx
  )

}



stan_data<-real_data_clean(S_df,X,grid,R=R,dr=dr,dx=dx)

model<-stan_model('data-analysis/lgcp_w_covariates.stan')

iters<-5000
burn<-1000
post <- sampling(
  model,
  data = stan_data,
  chains = 1,
  warmup = burn,
  iter = iters,
  cores = 1,
  refresh = 100,
  init=list(list(beta=c(-10,0,0),lsig2=0,lphi=4.14))
)
saveRDS(post,'sim_output/bei/post.rds')
post<-readRDS('sim_output/bei/post.rds')

post_mean<-get_posterior_mean(post,pars=c("beta","lsig2","lphi"))

bei.extra_scaled<-bei.extra
bei.extra_scaled$elev<-(bei.extra$elev-mean(bei.extra$elev))/sd(bei.extra$elev)
bei.extra_scaled$grad<-(bei.extra$grad-mean(bei.extra$grad))/sd(bei.extra$grad)

mu<-post_mean[1]+post_mean[2]*bei.extra_scaled$elev+post_mean[3]*bei.extra_scaled$grad


all_out<-list()
times<-NULL
all_out<-readRDS('sim_output/bei/all_out.rds')
times<-readRDS('sim_output/bei/times.rds')
for(j in c(9:20)){
  print(j)
  t<-system.time({
    S<-rLGCP(nsim=5,model="exponential",mu=mu,param=list(var=exp(post_mean[4]),scale=exp(post_mean[5])))
    stan_data_sub <- foreach(i=c(1:5))%dopar%{
      S_df<-cbind(S[[i]]$x,S[[i]]$y)
      real_data_clean(S_df,X,grid,R=R,dr=dr,dx=dx)
    }
    
    
    out<-foreach(stan_data_sub_k=stan_data_sub)%dopar%{
  
      post_boot <- sampling(
        model,
        data = stan_data_sub_k,
        chains = 1,
        warmup = burn,
        iter = iters,
        cores = 1,
        refresh = 0,
        init=list(list(beta=c(-10,0,0),lsig2=0,lphi=4.14))
      )
      
      summ<-summary(post_boot,pars=c("beta","lsig2","lphi"),prob=c(0.025,0.975))$summary
      t<-sum(get_elapsed_time(post_boot))
      rm(post_boot)
      return(list(t,summ[,c(1,4,5)]))
    }
    all_out<-c(all_out,out)
  })
  times[j]<-t[3]
  saveRDS(times,'sim_output/bei/times.rds')
  saveRDS(all_out,'sim_output/bei/all_out.rds')
}
