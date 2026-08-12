data {
  int<lower=1> N;
  int<lower=1> P;

  matrix[N, P] X;
  vector[N] Y;
  vector[N] D;
}

parameters {

  // nuisance regressions
  vector[P] delta;
  vector[P] gamma;

  // covariance structure
  vector<lower=0>[2] sigma;
  cholesky_factor_corr[2] Lcorr;
}

transformed parameters {
  matrix[2,2] Sigma;

  Sigma =
    diag_pre_multiply(sigma, Lcorr) *
    diag_pre_multiply(sigma, Lcorr)';
}

model {

  vector[N] mu_y;
  vector[N] mu_d;

  // -------------------------
  // Priors (BDML-Basic style)
  // -------------------------
  delta ~ normal(0, 5);
  gamma ~ normal(0, 5);

  sigma ~ cauchy(0, 2.5);
  Lcorr ~ lkj_corr_cholesky(4);

  // -------------------------
  // Means
  // -------------------------
  mu_y = X * delta;
  mu_d = X * gamma;

  // -------------------------
  // Likelihood
  // -------------------------
  for (n in 1:N) {
    vector[2] z;
    vector[2] mu;

    z[1] = Y[n];
    z[2] = D[n];

    mu[1] = mu_y[n];
    mu[2] = mu_d[n];

    target += multi_normal_lpdf(z | mu, Sigma);
  }
}

generated quantities {
  real alpha;

  // causal parameter from covariance
  alpha = Sigma[1,2] / Sigma[2,2];
}