functions {

  vector palm_intensity_calc(vector d, real rho, real sig2, real phi){
    int N = num_elements(d);
    vector[N] lambda;
    lambda = rho*exp(sig2*exp(-d/phi));
    return lambda;
  }

  real palm_loglik(vector d1, vector dG, real rho, real sig2, real phi, real dr, real NR){
    real pll;
    real loglam1;
    real lam2;
    loglam1 = sum(2*log(palm_intensity_calc(d1, rho, sig2, phi)));
    lam2 = NR*2*pi()*sum(palm_intensity_calc(dG, rho, sig2, phi).*dG*dr);
    pll = loglam1 - lam2;
    return pll;
  }
}
data {
  int<lower=0> N;          // number of points
  real rho_sd;          
  int<lower=0> N11;          // number of point pairs (symmetric)
  int<lower=0> N2;          // number of quadrature distances
  int<lower=0> NR;          // number of points inside window
  real<lower=0> dr;          // stepsize for compensator approx
  vector<lower=0>[N11] d1;               // point pair distances(symmetric)
  vector<lower=0>[N2] dG;            // quadrature distances
}
parameters {
  real rho;
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
  target += palm_loglik(d1, dG, rho, sig2, phi, dr, NR);
  target += normal_lpdf(rho | N, rho_sd);
  target += normal_lpdf(lsig2 | 0, 5);
  target += uniform_lpdf(lphi | -3, -1.6);
}
