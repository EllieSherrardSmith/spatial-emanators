// bernoulli_logistic transformed data function
data {
  
  int<lower=1> N;                  // rows of data
  
  int<lower=0> n_t[N];             // Total number of mosquitoes entering treated huts
  int<lower=0> duf_t[N];             // Number mosquites dead unfed sprayed hut
  int<lower=0> df_t[N];             // Number mosquites dead fed sprayed hut
  int<lower=0> fed_t[N];           // Number of mosquitoes feeding in treated huts assuming equal feeding for dead and alive ones
  int<lower=0> deterrence[N];  // Number of mosquitoes in sprayed huts
  int<lower=0> deterrence_total[N]; //Total number of mosquitoes in both sprayed and control huts
  
  // int<lower=1> N_treatments;           // treatments (control and trt)
  int<lower=0> treatments[N];
  vector<lower=0>[N] time;       // predictor

}

parameters {
  //Consider death. This is the proportion of mosquitoes dying (d_t) in treated huts (n_t)
  real alpha1;
  real alpha2;
  real alpha3;
  
  //Consider death and fed. This is the proportion of mosquitoes dying (d_t) in treated huts (n_t)
  real alphaf1;
  real alphaf2;
  real alphaf3;
  
  //Consider feeding. This is the proportion of mosquitoes that successfully fed in treatment (f_t)
  real beta1;
  real beta2;
  real beta3;
  //real beta4;
  
  //Consider deterrence. This is the proportion of mosquitoes deterred in treatment (d_t)
  real omega1;
  real omega2;
  
  // Overdispersion (concentration) parameters
  real<lower=0> phi_duf;     // for death
  real<lower=0> phi_df;     // for death
  real<lower=0> phi_f;     // for feeding
  real<lower=0> phi_det;   // for deterrence
  
}

model {
  vector[N] logit_p_duf;
  vector[N] logit_p_df;
  vector[N] logit_p_f;
  vector[N] logit_p_det;

  // Priors
  alpha1 ~ normal(0, 10);
  alpha2 ~ normal(0, 10);
  alpha3 ~ normal(0, 10);

  alphaf1 ~ normal(0, 10);
  alphaf2 ~ normal(0, 10);
  alphaf3 ~ normal(0, 10);

  beta1 ~ normal(0, 1);
  beta2 ~ normal(0, 10);
  beta3 ~ normal(0, 1);
  // beta4 ~ normal(0, 10);

  omega1 ~ normal(0, 10);
  omega2 ~ normal(0, 10);

  phi_duf ~ exponential(1);
  phi_df ~ exponential(1);
  phi_f ~ exponential(1);
  phi_det ~ exponential(1);

  for (n in 1:N) {
    logit_p_duf[n] = alpha1 + alpha3 * treatments[n] + alpha2 * time[n];
    logit_p_df[n] = alphaf1 + alphaf3 * treatments[n] + alphaf2 * time[n];
    logit_p_f[n] = beta1 + beta3 * treatments[n] + beta2 * time[n]; //+ beta4 * time[n] * treatments[n];
    logit_p_det[n] = omega1 + omega2 * time[n];
  }

  for (n in 1:N) {
    // Convert logit_p to mean p, clipped to avoid 0 or 1
      real p_duf = fmin(fmax(inv_logit(logit_p_duf[n]), 1e-6), 1 - 1e-6);
      real p_df = fmin(fmax(inv_logit(logit_p_df[n]), 1e-6), 1 - 1e-6);
      real p_f = fmin(fmax(inv_logit(logit_p_f[n]), 1e-6), 1 - 1e-6);
      real p_det = fmin(fmax(inv_logit(logit_p_det[n]), 1e-6), 1 - 1e-6);

    // Get alpha/beta shape parameters
    real alpha_d = p_duf * phi_duf;
    real beta_d  = (1 - p_duf) * phi_duf;

    real alpha_df = p_df * phi_df;
    real beta_df  = (1 - p_df) * phi_df;

    real alpha_f = p_f * phi_f;
    real beta_f  = (1 - p_f) * phi_f;

    real alpha_det = p_det * phi_det;
    real beta_det  = (1 - p_det) * phi_det;

    // Likelihoods using beta-binomial
    duf_t[n] ~ beta_binomial(n_t[n], alpha_d, beta_d);
    df_t[n] ~ beta_binomial(n_t[n], alpha_df, beta_df);
    fed_t[n] ~ beta_binomial(n_t[n], alpha_f, beta_f);
    deterrence[n] ~ beta_binomial(deterrence_total[n], alpha_det, beta_det);
  }
}

generated quantities {
  int duf_t_rep[N];       // Replicated death counts
  int df_t_rep[N];       // Replicated death counts
  int fed_t_rep[N];     // Replicated fed counts
  int det_rep[N];       // Replicated deterrence counts

  real log_lik_duf[N];    // Log-likelihood for deaths
  real log_lik_df[N];    // Log-likelihood for deaths
  real log_lik_f[N];    // Log-likelihood for feeding
  real log_lik_det[N];  // Log-likelihood for deterrence

  for (n in 1:N) {
    // Stabilized p values
    real p_duf = fmin(fmax(inv_logit(alpha1 + alpha3 * treatments[n] + alpha2 * time[n]), 1e-6), 1 - 1e-6);
    real p_df = fmin(fmax(inv_logit(alphaf1 + alphaf3 * treatments[n] + alphaf2 * time[n]), 1e-6), 1 - 1e-6);
    real p_f = fmin(fmax(inv_logit(beta1 + beta3 * treatments[n] + beta2 * time[n]), 1e-6), 1 - 1e-6);// + beta4 * time[n] * treatments[n]), 1e-6), 1 - 1e-6);
    real p_det = fmin(fmax(inv_logit(omega1 + omega2 * time[n]), 1e-6), 1 - 1e-6);

    // Beta parameters
    real alpha_d = p_duf * phi_duf;
    real beta_d  = (1 - p_duf) * phi_duf;

    real alpha_df = p_df * phi_df;
    real beta_df  = (1 - p_df) * phi_df;

    real alpha_f = p_f * phi_f;
    real beta_f  = (1 - p_f) * phi_f;

    real alpha_det = p_det * phi_det;
    real beta_det  = (1 - p_det) * phi_det;

    // Posterior predictive draws
    duf_t_rep[n] = beta_binomial_rng(n_t[n], alpha_d, beta_d);
    df_t_rep[n] = beta_binomial_rng(n_t[n], alpha_df, beta_df);
    fed_t_rep[n] = beta_binomial_rng(n_t[n], alpha_f, beta_f);
    det_rep[n] = beta_binomial_rng(deterrence_total[n], alpha_det, beta_det);

    // Log-likelihoods
    log_lik_duf[n] = beta_binomial_lpmf(duf_t[n] | n_t[n], alpha_d, beta_d);
    log_lik_df[n] = beta_binomial_lpmf(df_t[n] | n_t[n], alpha_df, beta_df);
    log_lik_f[n] = beta_binomial_lpmf(fed_t[n] | n_t[n], alpha_f, beta_f);
    log_lik_det[n] = beta_binomial_lpmf(deterrence[n] | deterrence_total[n], alpha_det, beta_det);
  }
}

