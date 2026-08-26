functions {
  real palm_loglik(vector dS_counts, vector dG, real dx, matrix t_X, matrix mat_counts, vector X_counts, row_vector beta, real sig2, real phi){
    int N = num_elements(dG);
    int xN = cols(t_X);
    real pll;
    real eq1;
    real comp;
    vector[N] g;
    row_vector[xN] lam1;
    g = exp(sig2*exp(-dG/phi));
    lam1 = exp(beta*t_X + sig2/2);
    eq1 = dot_product(X_counts,log(lam1)) + 2*dot_product(dS_counts,log(g));
    comp = dot_product(lam1,mat_counts*g)*square(dx);
    pll = eq1 - comp;
    return pll;
  }
}
data {
  int<lower=0> dN;          // number of bins
  int<lower=0> p; // number of covariates
  int<lower=0> xN; // number of covariate tiles
  vector<lower=0>[dN] dG;          // bin centroids
  vector<lower=0>[dN] dS_counts;   
  matrix[p, xN] t_X;
  vector[xN] X_counts;
  matrix[xN, dN] mat_counts;
  real<lower=0> dx; 
}
parameters {
  row_vector[p] beta;
  real lsig2;
  real lphi;
}
transformed parameters {
  real<lower=0> sig2;
  real<lower=0> phi;
  sig2 = exp(lsig2);
  phi = exp(lphi);
}
model {
  target += palm_loglik(dS_counts, dG, dx, t_X, mat_counts, X_counts, beta, sig2, phi);
  target += normal_lpdf(beta | 0, 31.62);
  target += normal_lpdf(lsig2 | 0, 3);
  target += normal_lpdf(lphi | 4.6, 0.5);
}
