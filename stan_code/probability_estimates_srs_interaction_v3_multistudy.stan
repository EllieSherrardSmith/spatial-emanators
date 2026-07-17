data {

  // DATASET 1
  int<lower=1> N1;

  int<lower=0> n_t1[N1];
  int<lower=0> duf_t1[N1];
  int<lower=0> df_t1[N1];
  int<lower=0> fed_t1[N1];
  int<lower=0> deterrence1[N1];
  int<lower=0> deterrence_total1[N1];

  int<lower=0> treatments1[N1];
  vector<lower=0>[N1] time1;

  // DATASET 2
  int<lower=1> N2;

  int<lower=0> n_t2[N2];
  int<lower=0> duf_t2[N2];
  int<lower=0> df_t2[N2];
  int<lower=0> fed_t2[N2];
  int<lower=0> deterrence2[N2];
  int<lower=0> deterrence_total2[N2];

  int<lower=0> treatments2[N2];
  vector<lower=0>[N2] time2;
}

parameters {

  // dead unfed
  real alpha1;
  real alpha2;
  real alpha3;

  // dead fed
  real alphaf1;
  real alphaf2;
  real alphaf3;

  // feeding
  real beta1;
  real beta2;
  real beta3;

  // deterrence
  real omega1;
  real omega2;

  // beta-binomial concentrations
  real<lower=1e-3> phi_duf;
  real<lower=1e-3> phi_df;
  real<lower=1e-3> phi_f;
  real<lower=1e-3> phi_det;
}

model {

  // Priors
  alpha1 ~ normal(0,10);
  alpha2 ~ normal(0,10);
  alpha3 ~ normal(0,10);

  alphaf1 ~ normal(0,10);
  alphaf2 ~ normal(0,10);
  alphaf3 ~ normal(0,10);

  beta1 ~ normal(0,1);
  beta2 ~ normal(0,10);
  beta3 ~ normal(0,1);

  omega1 ~ normal(0,10);
  omega2 ~ normal(0,10);

  phi_duf ~ exponential(1);
  phi_df  ~ exponential(1);
  phi_f   ~ exponential(1);
  phi_det ~ exponential(1);

  // ---------------------------
  // DATASET 1
  // ---------------------------

  for (n in 1:N1) {

    real p_duf =
      fmin(fmax(inv_logit(alpha1 +
                          alpha3 * treatments1[n] +
                          alpha2 * time1[n]),
                1e-6), 1 - 1e-6);

    real p_df =
      fmin(fmax(inv_logit(alphaf1 +
                          alphaf3 * treatments1[n] +
                          alphaf2 * time1[n]),
                1e-6), 1 - 1e-6);

    real p_f =
      fmin(fmax(inv_logit(beta1 +
                          beta3 * treatments1[n] +
                          beta2 * time1[n]),
                1e-6), 1 - 1e-6);

    real p_det =
      fmin(fmax(inv_logit(omega1 +
                          omega2 * time1[n]),
                1e-6), 1 - 1e-6);

    real alpha_d  = p_duf * phi_duf;
    real beta_d   = (1 - p_duf) * phi_duf;

    real alpha_df_ = p_df * phi_df;
    real beta_df_  = (1 - p_df) * phi_df;

    real alpha_f_ = p_f * phi_f;
    real beta_f_  = (1 - p_f) * phi_f;

    real alpha_det = p_det * phi_det;
    real beta_det  = (1 - p_det) * phi_det;

    duf_t1[n] ~ beta_binomial(n_t1[n], alpha_d, beta_d);
    df_t1[n] ~ beta_binomial(n_t1[n], alpha_df_, beta_df_);
    fed_t1[n] ~ beta_binomial(n_t1[n], alpha_f_, beta_f_);
    deterrence1[n] ~ beta_binomial(deterrence_total1[n],
                                   alpha_det,
                                   beta_det);
  }

  // ---------------------------
  // DATASET 2
  // ---------------------------

  for (n in 1:N2) {

    real p_duf =
      fmin(fmax(inv_logit(alpha1 +
                          alpha3 * treatments2[n] +
                          alpha2 * time2[n]),
                1e-6), 1 - 1e-6);

    real p_df =
      fmin(fmax(inv_logit(alphaf1 +
                          alphaf3 * treatments2[n] +
                          alphaf2 * time2[n]),
                1e-6), 1 - 1e-6);

    real p_f =
      fmin(fmax(inv_logit(beta1 +
                          beta3 * treatments2[n] +
                          beta2 * time2[n]),
                1e-6), 1 - 1e-6);

    real p_det =
      fmin(fmax(inv_logit(omega1 +
                          omega2 * time2[n]),
                1e-6), 1 - 1e-6);

    real alpha_d  = p_duf * phi_duf;
    real beta_d   = (1 - p_duf) * phi_duf;

    real alpha_df_ = p_df * phi_df;
    real beta_df_  = (1 - p_df) * phi_df;

    real alpha_f_ = p_f * phi_f;
    real beta_f_  = (1 - p_f) * phi_f;

    real alpha_det = p_det * phi_det;
    real beta_det  = (1 - p_det) * phi_det;

    duf_t2[n] ~ beta_binomial(n_t2[n], alpha_d, beta_d);
    df_t2[n] ~ beta_binomial(n_t2[n], alpha_df_, beta_df_);
    fed_t2[n] ~ beta_binomial(n_t2[n], alpha_f_, beta_f_);
    deterrence2[n] ~ beta_binomial(deterrence_total2[n],
                                   alpha_det,
                                   beta_det);
  }
}

generated quantities {

  vector[N1] log_lik1;
  vector[N2] log_lik2;

  for (n in 1:N1) {

    real p_duf =
      fmin(fmax(inv_logit(alpha1 +
                          alpha3 * treatments1[n] +
                          alpha2 * time1[n]),
                1e-6), 1 - 1e-6);

    real alpha_d = p_duf * phi_duf;
    real beta_d  = (1 - p_duf) * phi_duf;

    log_lik1[n] =
      beta_binomial_lpmf(duf_t1[n] |
                         n_t1[n],
                         alpha_d,
                         beta_d);
  }

  for (n in 1:N2) {

    real p_duf =
      fmin(fmax(inv_logit(alpha1 +
                          alpha3 * treatments2[n] +
                          alpha2 * time2[n]),
                1e-6), 1 - 1e-6);

    real alpha_d = p_duf * phi_duf;
    real beta_d  = (1 - p_duf) * phi_duf;

    log_lik2[n] =
      beta_binomial_lpmf(duf_t2[n] |
                         n_t2[n],
                         alpha_d,
                         beta_d);
  }
}

