functions {
  vector palm_intensity_calc(vector d, real rho, real alpha){
    int N = num_elements(d);
    vector[N] lambda;
    lambda = rho*(1 - square(exp(-square(d/alpha))));
    return lambda;
  }

  real palm_loglik(vector dS_counts, vector dG, vector dG_counts, real dx, real rho, real alpha, int dN, real eta){
    real pll;
    real loglam1;
    real lam2;
    vector[dN] lam_dG;
    lam_dG = palm_intensity_calc(dG, rho, alpha);
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
  real lalpha;
}
transformed parameters {
  real<lower=0> alpha;
  alpha = exp(lalpha);
}
model {
  target += palm_loglik(dS_counts, dG, dG_counts, dx, rho, alpha, dN, eta);
  target += normal_lpdf(rho | rho_mean, rho_sd);
  target += normal_lpdf(lalpha | -3, 3);
}
