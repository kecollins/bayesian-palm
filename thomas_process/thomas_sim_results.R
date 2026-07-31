
apply(t(matrix(unlist(lapply(temp,\(x) x$cal2_post_output[,1]-c(log(mu),log(nu),log(sig2)))),nrow=3)),2,mean)

apply(t(matrix(unlist(lapply(temp,\(x) (x$palm_post_output[,4]<c(log(mu),log(nu),log(sig2)))*(x$palm_post_output[,5]>c(log(mu),log(nu),log(sig2))))),nrow=3)),2,mean)
apply(t(matrix(unlist(lapply(temp,\(x) (x$cal1_post_output[,4]<c(log(mu),log(nu),log(sig2)))*(x$cal1_post_output[,5]>c(log(mu),log(nu),log(sig2))))),nrow=3)),2,mean)
apply(t(matrix(unlist(lapply(temp,\(x) (x$cal2_post_output[,2]<c(log(mu),log(nu),log(sig2)))*(x$cal2_post_output[,3]>c(log(mu),log(nu),log(sig2))))),nrow=3)),2,mean)
      
lapply(temp,\(x) x$cal1_post_output)
