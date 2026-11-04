data {
  
  // ---- within-hut compositional data (one row per hut-night) ----
  int<lower=1> J;
  array[J, 4] int<lower=0> y;            // [ ,1]=unfed-alive [ ,2]=unfed-dead [ ,3]=fed-dead [ ,4]=fed-alive
  array[J] int<lower=0, upper=1> arm;    // 0 = control, 1 = treatment
  vector[J] time;                        // time since deployment, per hut-row

  // ---- deterrence (entry) data, aggregated to trial-day level ----
  int<lower=1> K;                        // number of trial-day units with both arms present
  array[K] int<lower=0> Tcatch;          // total catch summed across treatment huts, trial-day k
  array[K] int<lower=0> Ccatch;          // total catch summed across control huts, trial-day k
  vector[K] time_det;                    // time for trial-day k (same scale as `time` above)
  vector[K] log_hut_ratio;               // log(n_treatment_huts / n_control_huts) that day; 0 if balanced

  // ---- link from each hut-row to its trial-day deterrence group (for generated quantities) ----
  array[J] int<lower=1, upper=K> det_idx;

  real<lower=0> pi_prior_a;
  real<lower=0> pi_prior_b;
}
parameters {
  real<lower=0, upper=1> pi;

  real b0_R; real b1_R; real b2_R; real b3_R;
  real b0_K; real b1_K; real b2_K; real b3_K;
  real b0_M; real b1_M; real b2_M; real b3_M;

  real<lower=0> phi;                     // shared Dirichlet-multinomial precision (within-hut composition)

  // Deterrence: treatment-vs-control contrast in entry, control arm as the reference (delta_control == 0)
  real g0;                               // deterrence intercept (log-odds scale, at time 0)
  real g1;                               // deterrence time trend
  real<lower=0, upper=1> rho_delta;      // overdispersion, ICC scale, for the deterrence beta-binomial
}
transformed parameters {
  real<lower=0> phi_delta = (1 - rho_delta) / rho_delta;
}
model {
  // ---- priors ----
  pi ~ beta(pi_prior_a, pi_prior_b);

  b0_R ~ normal(0, 2); b1_R ~ normal(0, 2); b2_R ~ normal(0, 1);
  b0_K ~ normal(0, 2); b1_K ~ normal(0, 2); b2_K ~ normal(0, 1);
  b0_M ~ normal(0, 2); b1_M ~ normal(0, 2); b2_M ~ normal(0, 1);
  phi ~ gamma(2, 0.1);   // weakly informative; review against expected catch sizes

  g0 ~ normal(0, 2);
  g1 ~ normal(0, 1);
  rho_delta ~ beta(2, 5);

  // ---- within-hut compositional likelihood (entry-conditional, Dirichlet-multinomial) ----
  for (j in 1:J) {
    real muR = inv_logit(b0_R + b1_R * arm[j] + b2_R * arm[j] * time[j]);
    real muK = inv_logit(b0_K + b1_K * arm[j] + b2_K * arm[j] * time[j]);
    real muM = inv_logit(b0_M + b1_M * arm[j] + b2_M * arm[j] * time[j]);

    real kappa = (arm[j] == 0) ? 0.0 : muK;

    vector[4] p;
    p[1] = (1 - pi) * (1 - kappa) + pi * muR;
    p[2] = (1 - pi) * kappa + pi * (1 - muR) * muK;
    p[3] = pi * (1 - muR) * (1 - muK) * muM;
    p[4] = pi * (1 - muR) * (1 - muK) * (1 - muM);

    y[j] ~ dirichlet_multinomial(p * phi);
  }

  // ---- deterrence (entry) likelihood, on matched trial-day totals ----
  // Conditioning on the combined catch (Tcatch[k] + Ccatch[k]) each trial-day cancels out the
  // (highly variable, unmodelled) nightly mosquito abundance, leaving only the treatment-vs-control
  // share to explain -- the same logic as a matched-pairs / conditional logistic analysis.
  for (k in 1:K) {
    real eta_k = log_hut_ratio[k] + g0 + g1 * time_det[k];
    real p_k = inv_logit(eta_k);          // expected share of the combined catch entering treatment huts
    int total_k = Tcatch[k] + Ccatch[k];
    Tcatch[k] ~ beta_binomial(total_k, p_k * phi_delta, (1 - p_k) * phi_delta);
  }
}
generated quantities {
  // Deterrence at trial-day k, relative to the control-arm reference (delta_control == 0 by convention,
  // mirroring the kappa_control == 0 constraint). Negative values indicate net attraction.
  vector[K] delta = 1 - exp(g0 + g1 * time_det);

  // Fully unconditional-on-entry outcome probabilities per hut-row: control rows equal the
  // within-hut (fixed-effects) probabilities directly (delta==0 by definition); treatment rows are
  // scaled by (1-delta) for that row's trial-day, with the deterred mass reported separately.
  vector[J] p_unfed_alive_uncond;
  vector[J] p_unfed_dead_uncond;
  vector[J] p_fed_dead_uncond;
  vector[J] p_fed_alive_uncond;
  vector[J] p_deterred_uncond;

  for (j in 1:J) {
    real muR = inv_logit(b0_R + b1_R * arm[j] + b2_R * arm[j] * time[j]); // Remove one of these terms.
    real muK = inv_logit(b0_K + b1_K * arm[j] + b2_K * arm[j] * time[j]);
    real muM = inv_logit(b0_M + b1_M * arm[j] + b2_M * arm[j] * time[j]);

    real kappa = (arm[j] == 0) ? 0.0 : muK;
    real one_minus_delta = (arm[j] == 0) ? 1.0 : (1 - delta[det_idx[j]]);

    real p1 = (1 - pi) * (1 - kappa) + pi * muR;
    real p2 = (1 - pi) * kappa + pi * (1 - muR) * muK;
    real p3 = pi * (1 - muR) * (1 - muK) * muM;
    real p4 = pi * (1 - muR) * (1 - muK) * (1 - muM);

    p_unfed_alive_uncond[j] = one_minus_delta * p1;
    p_unfed_dead_uncond[j]  = one_minus_delta * p2;
    p_fed_dead_uncond[j]    = one_minus_delta * p3;
    p_fed_alive_uncond[j]   = one_minus_delta * p4;
    p_deterred_uncond[j]    = 1 - one_minus_delta;
  }
}
