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
  // array[N] vector[N] d;           // grid centroid distances
  matrix[N,N] d;
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
  vector[N] W0;
  sig2 = exp(lsig2);
  phi = exp(lphi);
  //W0 = W - mean(W);
}
model {
  matrix[N,N] Sigma;
  //Sigma = gp_exponential_cov(d,sqrt(sig2),phi);
  //for (i in 1:N) {
  //  for (j in 1:N) {
  //    Sigma[i,j] = sig2 * exp(-d[i][j] / phi);
  //  }
  //}
  Sigma = sig2*exp(-d/phi);
  target += loglik(y, W, dx, mu);
  target += normal_lpdf(mu | 0, 100);
  target += multi_normal_lpdf(W | rep_vector(0,N), Sigma);
  target += normal_lpdf(lsig2 | 0, 3);
  target += normal_lpdf(lphi | -2.3, 0.3);
}
