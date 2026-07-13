functions {

  vector palm_intensity_calc(vector d, real rho, real sig2, real phi){
    int N = num_elements(d);
    vector[N] lambda;
    lambda = rho*exp(sig2*exp(-d/phi));
    return lambda;
  }

  real palm_loglik(vector dS_counts, vector dG, vector dG_counts, real dx, real rho, real sig2, real phi){
    real pll;
    real loglam1;
    real lam2;
    loglam1 = sum(2*dS_counts.*log(palm_intensity_calc(dG, rho, sig2, phi)));
    lam2 = sum(dG_counts.*palm_intensity_calc(dG, rho, sig2, phi)*dx);
    pll = loglam1 - lam2;
    return pll;
  }
}
data {
  real<lower=0> rho_mean;          // points per unit
  real<lower=0> rho_sd;     
  int<lower=0> dN;          // number of bins
  vector<lower=0>[dN] dG;          // bin centroids
  vector<lower=0>[dN] dS_counts;   
  vector<lower=0>[dN] dG_counts;  
  real<lower=0> dx; 
}
parameters {
  real<lower=0> rho;
  real lsig2;
  real lphi;
}
transformed parameters {
  real<lower=0> sig2;
  real<lower=0> phi;
  real<lower=0> mu;
  sig2 = exp(lsig2);
  phi = exp(lphi);
  mu = log(rho)-sig2/2;
}
model {
  target += palm_loglik(dS_counts, dG, dG_counts, dx, rho, sig2, phi);
  target += normal_lpdf(rho | rho_mean, rho_sd);
  target += normal_lpdf(lsig2 | 0, 3);
  target += normal_lpdf(lphi | -2.3, 0.3);
}
