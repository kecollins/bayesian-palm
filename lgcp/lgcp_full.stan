functions {
  real loglik(vector y, vector W, real dx, real mu){
    real pll;
    pll = sum(y.*(mu+W) - dx*exp(mu+W));
    return pll;
  }
}
data {
  int<lower=0> N;          // number of grid cells
  real<lower=0> dx;        // grid cell size
  vector[N] y;             // cell counts
  array[N] vector[N] d;           // grid centroid distances
}
parameters {
  real mu;
  real lsig2;
  real lphi;
  vector[N] W;
}
transformed parameters {
  real<lower=0> sig2;
  real<lower=0> phi;
  matrix[N,N] Sigma;
  sig2 = exp(lsig2);
  phi = exp(lphi);
  Sigma = gp_exponential_cov(d, sig2, phi);
}
model {
  target += loglik(y, W, dx, mu);
  target += normal_lpdf(mu | 0, 100);
  target += multi_normal_lpdf(W | rep_vector(0,N), Sigma);
  target += normal_lpdf(lsig2 | 0, 3);
  target += uniform_lpdf(lphi | -2.3, 0.1);
}
