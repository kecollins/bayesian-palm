functions {

  vector palm_intensity_calc(vector d, real rho, real nu, real sig2){
    int N = num_elements(d);
    vector[N] lambda;
    lambda = rho+nu/(4*square(pi())*sig2)*exp(-square(d)/(4*sig2));
    return lambda;
  }

  real palm_loglik(vector d1, vector d2, vector dG, real rho, real nu, real sig2, real dr, real NR){
    real pll;
    real loglam1;
    real lam2;
    loglam1 = sum(2*log(palm_intensity_calc(d1, rho, sig2, phi)))+sum(log(palm_intensity_calc(d2, rho, sig2, phi)));
    lam2 = NR*2*pi()*sum(palm_intensity_calc(dG, rho, sig2, phi).*dG*dr);
    pll = loglam1 - lam2;
    return pll;
  }
}
data {
  real<lower=0> rho_mean;          // points per unit
  real rho_sd;     
  int<lower=0> N11;          // number of point pairs (symmetric)
  int<lower=0> N12;          // number of point pairs
  int<lower=0> N2;          // number of quadrature distances
  int<lower=0> NR;          // number of points inside window
  real<lower=0> dr;          // stepsize for compensator approx
  vector<lower=0>[N11] d1;               // point pair distances(symmetric)
  vector<lower=0>[N12] d2;               // point pair distances
  vector<lower=0>[N2] dG;            // quadrature distances
}
parameters {
  real rho;
  real nu;
  real lsig2;
}
transformed parameters {
  real<lower=0> sig2;
  real<lower=0> phi;
  real mu;
  sig2 = exp(lsig2);
  phi = exp(lphi);
  mu = rho/nu;
}
model {
  target += palm_loglik(d1, d2, dG, rho, nu, sig2, dr, NR);
  target += normal_lpdf(rho | rho_mean, rho_sd);
  target += normal_lpdf(lsig2 | 0, 3);
  target += normal_lpdf(lphi | -3, -1.6);
}
