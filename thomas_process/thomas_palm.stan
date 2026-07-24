functions {

  vector palm_intensity_calc(vector d, real rho, real nu, real sig2){
    int N = num_elements(d);
    vector[N] lambda;
    lambda = rho+nu/(4*pi()*sig2)*exp(-square(d)/(4*sig2));
    return lambda;
  }

  real palm_loglik(vector dS_counts, vector dG, vector dG_counts, real dx, real rho, real nu, real sig2, int dN, real eta){
    real pll;
    real loglam1;
    real lam2;
    vector[dN] lam_dG;
    lam_dG = palm_intensity_calc(dG, rho, nu, sig2);
    loglam1 = sum(2*dS_counts.*log(lam_dG));
    lam2 = sum(dG_counts.*lam_dG*dx);
    pll = eta*(loglam1 - lam2);
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
  real<lower=0> eta;
}
parameters {
  real rho;
  real lnu;
  real lsig2;
}
transformed parameters {
  real<lower=0> sig2;
  real<lower=0> nu;
  real mu;
  sig2 = exp(lsig2);
  nu = exp(lnu);
  mu = rho/nu;
}
model {
  target += palm_loglik(dS_counts, dG, dG_counts, dx, rho, nu, sig2, dN, eta);
  target += normal_lpdf(rho | rho_mean, rho_sd);
  target += normal_lpdf(lnu | 0, 10);
  target += normal_lpdf(lsig2 | 0, 3);
}
