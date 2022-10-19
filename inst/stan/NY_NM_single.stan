data {
  // number of observations
  int<lower=0> N;
  // number of columns in design matrix excluding A (and M)
  int<lower=0> P;
  // design matrix, excluding treatment A
  matrix[N, P] X;
  // observed treatment
  vector[N] A;
  // observed mediator
  vector [N] M;
  // outcome
  vector [N] Y;
  // mean of regression priors
  vector[P + 1] location_m;
  vector[P + 2] location_y;
  // variance-covariance of regression priors
  cov_matrix[P + 2] scale_y;
  cov_matrix[P + 1] scale_m;
  // scale parameters for standard devation for the mediator-outcome model
  real<lower = 0> scale_sd_m;
  real<lower = 0> scale_sd_y;
}
transformed data {
  // make vector of 1/N for (classical) bootstrapping
  vector[N] boot_probs = rep_vector(1.0/N, N);
  // make vector version of M
  vector[N] Mv = to_vector(M);
}
parameters {
  // regression coefficients (outcome model)
  vector[P + 2] alpha;
  // regression coefficients (mediator model)
  vector[P + 1] beta;
  // residual standard devation for the mediator model
  real<lower = 0> sigma_m;
  real<lower = 0> sigma_y;

}

transformed parameters {
  // partial M coefficient parameters
  vector[P] betaZ = head(beta, P);
  real betaA = beta[P + 1];
  // partial Y coefficient parameters
  vector[P] alphaZ = head(alpha, P);
  real alphaA = alpha[P + 1];
  real alphaM = alpha[P + 2];
}
model {
  // priors on causal coefficients weakly informative for binary exposure
  alpha ~ multi_normal(location_m, scale_m);
  beta ~ multi_normal(location_y, scale_y);
  // prior for the residual standrd devation of the mediator model
  target += normal_lpdf(sigma_y | 0, scale_sd_y) - normal_lcdf(0 | 0, scale_sd_y);
  target += normal_lpdf(sigma_m | 0, scale_sd_m) - normal_lcdf(0 | 0, scale_sd_m);

  // likelihoods
  M ~ normal (X * betaZ + A * betaA, sigma_m);
  Y ~ normal (X * alphaZ + A * alphaA + Mv * alphaM, sigma_y);
}
generated quantities {
  // row index to be sampled for bootstrap
  int row_i;
  // calculate NDE in the bootstrapped sample
  real NDE_control = 0;
  real NDE_treated = 0;
  real NIE_control = 0;
  real NIE_treated = 0;
  real TE = 0;
  real ANDE = 0;
  real ANIE = 0;
  vector[N] M_a0;
  vector[N] M_a1;
  vector[N] Y_a1Ma0;
  vector[N] Y_a0Ma0;
  vector[N] Y_a1Ma1;
  vector[N] Y_a0Ma1;
  for (n in 1:N) {
    // sample baseline covariates
    row_i = categorical_rng(boot_probs);
    // sample Ma where a = 0
    M_a0[n] = normal_rng (X[row_i] * betaZ, sigma_m);
    // sample Ma where a = 1
    M_a1[n] = normal_rng (X[row_i] * betaZ + betaA, sigma_m);
    // sample Y_(a=1, M=M_0) and Y_(a=0, M=M_0)
    Y_a1Ma0[n] = normal_rng(X[row_i] * alphaZ + M_a0[n] * alphaM + alphaA, sigma_y);
    Y_a0Ma0[n] = normal_rng(X[row_i] * alphaZ + M_a0[n] * alphaM, scale_sd_y);
    // sample Y_(a=1, M=M_1) and Y_(a=0, M=M_1)
    Y_a1Ma1[n] = normal_rng(X[row_i] * alphaZ + M_a1[n] * alphaM + alphaA, sigma_y);
    Y_a0Ma1[n] = normal_rng(X[row_i] * alphaZ + M_a1[n] * alphaM, sigma_y);
    // add contribution of this observation to the bootstrapped NDE
    NDE_control = NDE_control + (Y_a1Ma0[n] - Y_a0Ma0[n])/N;//control
    NDE_treated = NDE_treated + (Y_a1Ma1[n] - Y_a0Ma1[n])/N;//treated
    //natural indirect effect Yi(t, Mi(1)) - Yi(t, Mi(0))
    NIE_treated = NIE_treated + (Y_a1Ma1[n] - Y_a1Ma0[n])/N;//treated
    NIE_control = NIE_control + (Y_a0Ma1[n] - Y_a0Ma0[n])/N;//control
    TE = NDE_control + NIE_treated;
    ANDE = (NDE_control + NDE_treated)/2;
    ANIE = (NIE_control + NIE_treated)/2;
  }
}
