functions {

  vector palm_intensity_calc(vector d, real rho, real sig2, real phi){
    int N = num_elements(d);
    vector[N] lambda;
    lambda = rho*exp(sig2*exp(-d/phi));
    return lambda;
  }

  real palm_loglik(vector dS, vector dG, real dx, real rho, real sig2, real phi, real eta){
    real pll;
    real loglam1;
    real lam2;
    loglam1 = sum(2*log(palm_intensity_calc(dS, rho, sig2, phi)));
    lam2 = sum(palm_intensity_calc(dG, rho, sig2, phi)*dx);
    pll = eta*(loglam1 - lam2);
    return pll;
  }
}
data {
  real<lower=0> rho_mean;          // points per unit
  real<lower=0> rho_sd;     
  int<lower=0> N1;          // number of point pairs
  int<lower=0> N2;          // number of quadrature pairs
  vector<lower=0>[N1] dS;               // point pair distances
  vector<lower=0>[N2] dG;            // quadrature distances
  real<lower=0> dx; 
  real<lower=0> eta;
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
  target += palm_loglik(dS, dG, dx, rho, sig2, phi, eta);
  target += normal_lpdf(rho | rho_mean, rho_sd);
  target += normal_lpdf(lsig2 | 0, 3);
  target += normal_lpdf(lphi | -2.3, 0.3);
}
