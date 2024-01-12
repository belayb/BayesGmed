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
  array[N] int<lower=0,upper=1> M;
  // outcome
  vector [N] Y;
  // mean of regression priors 
  vector[P + 2] location_y; 
  vector[P + 1] location_m; 
  vector[4] location_gamma;
  // variance-covariance of regression priors 
  cov_matrix[P + 2] scale_y; 
  cov_matrix[P + 1] scale_m;
  cov_matrix[4] scale_gamma;
  // scale parameter for residual error 
  real<lower=0> scale_sd_y;
}
transformed data { 
  // make vector of 1/N for (classical) bootstrapping 
  vector[N] boot_probs = rep_vector(1.0/N, N);
  // make vector version of M 
  vector[N] Mv = to_vector(M); //Take a note here - it not nessary 
  // design matrix for the unmeasured confounder model
    matrix[N, 2] Xu = append_col(rep_vector(1.0, N), A); 
}
parameters { 
  // regression coefficients (confounder model)
  vector[4] gamma;
  // regression coefficients (outcome model) 
  vector[P + 2] alpha;
  // regression coefficients (mediator model) 
  vector[P + 1] beta;
   // residual standard devation for the outcome model
  real<lower = 0> sigma_y;
}

transformed parameters { 
  // you may not need all this need
  // partial M coefficient parameters 
  vector[P] betaZ = head(beta, P); 
  real betaU = gamma[3];
  real betaA = beta[P + 1];
  // partial Y coefficient parameters 
  vector[P] alphaZ = head(alpha, P); 
  real alphaU = gamma[4]; 
  real alphaA = alpha[P + 1]; 
  real alphaM = alpha[P + 2];
  vector[2] gammaU = head(gamma, 2);
}
model { 
  // linear predictors 
  // U regression 
  vector[N] eta_u;
  // M regression, if U = 0
  vector[N] eta_mu0;
  // Y regression, if U = 0
  vector[N] eta_yu0;
  //log-likelihood contribution for U = 0 and U = 1 cases
  real ll_0;
  real ll_1;
  // calculate linear predictor for U = 0 case
  eta_u = Xu * gammaU;
  eta_mu0 = X * betaZ + A * betaA;
  eta_yu0 = X * alphaZ + A * alphaA + Mv * alphaM;

  // informative prior 
  alpha ~ multi_normal(location_y, scale_y); 
  beta ~ multi_normal(location_m, scale_m);
  gamma ~ multi_normal(location_gamma, scale_gamma); 
  // prior for the residual standrd devation of the outcome model 
   target += normal_lpdf(sigma_y | 0, scale_sd_y) - normal_lcdf(0 | 0, scale_sd_y);
  // likelihoods
  for (n in 1:N){
    //contribution if U = 0
    ll_0 = log_inv_logit(eta_mu0[n]) * M[n]+log1m_inv_logit(eta_mu0[n]) * (1 - M[n])+
           normal_lpdf(Y[n]|eta_yu0[n],sigma_y)+log1m_inv_logit(eta_u[n]);
    
    //contribution if U = 0

    ll_1 = log_inv_logit(eta_mu0[n] + betaU) * M[n]+log1m_inv_logit(eta_mu0[n] + betaU) * (1 - M[n])+
           normal_lpdf(Y[n]|eta_yu0[n] + alphaU,sigma_y)+log1m_inv_logit(eta_u[n]);
    
    // contribution is summation over U possibilities
    
    target += log_sum_exp(ll_0, ll_1);
  }

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
  vector[N] U;
  for (n in 1:N) { 
    // sample baseline covariates 
    row_i = categorical_rng(boot_probs);
    // sample U
    U[n] = bernoulli_logit_rng(Xu[row_i] * gammaU);
    // sample Ma where a = 0 
    M_a0[n] = bernoulli_logit_rng(X[row_i] * betaZ + U[n] * betaU);
    // sample Ma where a = 1 
    M_a1[n] = bernoulli_logit_rng(X[row_i] * betaZ + betaA + U[n] * betaU);
    // sample Y_(a=1, M=M_0) and Y_(a=0, M=M_0) 
    Y_a1Ma0[n] = normal_rng(X[row_i] * alphaZ + M_a0[n] * alphaM + alphaA + U[n] * alphaU, sigma_y); 
    Y_a0Ma0[n] = normal_rng(X[row_i] * alphaZ + M_a0[n] * alphaM + U[n] * alphaU, sigma_y);
    // sample Y_(a=1, M=M_1) and Y_(a=0, M=M_1) 
    Y_a1Ma1[n] = normal_rng(X[row_i] * alphaZ + M_a1[n] * alphaM + alphaA + U[n] * alphaU, sigma_y); 
    Y_a0Ma1[n] = normal_rng(X[row_i] * alphaZ + M_a1[n] * alphaM + U[n] * alphaU, sigma_y);
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
