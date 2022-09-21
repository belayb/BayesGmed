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
  int<lower=0,upper=1> M[N];
  // outcome
  int<lower=0,upper=1> Y[N];
  // mean of regression priors
  vector[P + 2] location_y;
  vector[P + 1] location_m;
  // variance-covariance of regression priors
  cov_matrix[P + 2] scale_y;
  cov_matrix[P + 1] scale_m;
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
  alpha ~ multi_normal(location_y, scale_y);
  beta ~ multi_normal(location_m, scale_m);
  // likelihoods
  M ~ bernoulli_logit (X * betaZ + A * betaA);
  Y ~ bernoulli_logit(X * alphaZ + A * alphaA + Mv * alphaM);
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
    M_a0[n] = bernoulli_logit_rng (X[row_i] * betaZ);
    // sample Ma where a = 1
    M_a1[n] = bernoulli_logit_rng (X[row_i] * betaZ + betaA);
    // sample Y_(a=1, M=M_0) and Y_(a=0, M=M_0)
    Y_a1Ma0[n] = bernoulli_logit_rng(X[row_i] * alphaZ + M_a0[n] * alphaM + alphaA);
    Y_a0Ma0[n] = bernoulli_logit_rng(X[row_i] * alphaZ + M_a0[n] * alphaM);
    // sample Y_(a=1, M=M_1) and Y_(a=0, M=M_1)
    Y_a1Ma1[n] = bernoulli_logit_rng(X[row_i] * alphaZ + M_a1[n] * alphaM + alphaA);
    Y_a0Ma1[n] = bernoulli_logit_rng(X[row_i] * alphaZ + M_a1[n] * alphaM);
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
