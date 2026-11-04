// Seasonal version of multinomial_hill_eps_v2.stan -- same tree, same eps-reparametrised
// Hill decays for mu_R(t), mu_K(t), mu_M(t) and deterrence, same Dirichlet-multinomial
// composition and beta-binomial deterrence likelihoods. The only change is that the
// constant biting probability pi is replaced by a seasonal pi(t), a Fourier series on
// the logit scale in hut-trial day (no calendar dates needed -- day 1 of the trial is
// simply the phase origin of the cycle):
//
//   logit pi(t) = alpha0 + sum_{h=1}^{H} [ beta_cos[h]*cos(2*pi*h*t/P) + beta_sin[h]*sin(2*pi*h*t/P) ]
//
// P = `period` (days per seasonal cycle, e.g. 365) and H = `n_harmonics` are data.
// H = 0 recovers the constant-pi model with pi = inv_logit(alpha0).
//
//   f(t) = f0 / (1 + rho*(t/T_ref)^n),  rho = 1/eps - 1   =>   f(T_ref)/f0 = eps exactly
//
//   Outcome probabilities (as v2, with pi -> pi(t)):
//   A  = 1 - (1-nu_u)*(1-mu_K),  hF = 1 - (1-nu_f)*(1-mu_M)
//   p[alive,unfed] = (1-A)*(1-hF) * ((1-pi(t)) + pi(t)*mu_R)
//   p[dead, unfed] = A + (1-A)*hF * ((1-pi(t)) + pi(t)*mu_R)
//   p[dead, fed]   = pi(t)*(1-mu_R)*(1-A)*hF
//   p[alive,fed]   = pi(t)*(1-mu_R)*(1-A)*(1-hF)

functions {
  // Hill / log-logistic decay, parametrised by the residual fraction eps remaining at
  // t = T_ref: f0 at t=0, f0*eps at t=T_ref, shape n.
  real hill_decay_eps(real t, real f0, real eps, real n, real T_ref) {
    real rho = 1 / eps - 1;
    return f0 / (1 + rho * pow(fmax(t, 0) / T_ref, n));
  }

  // Fourier design matrix: columns [cos_1, sin_1, cos_2, sin_2, ...].
  // Mirrored by fourier_basis() in R/prepare_stan_data_seasonal.R -- keep them in sync.
  matrix fourier_basis(vector t, real period, int H) {
    int N = num_elements(t);
    matrix[N, 2 * H] X;
    for (h in 1:H) {
      X[, 2 * h - 1] = cos(2 * pi() * h * t / period);
      X[, 2 * h]     = sin(2 * pi() * h * t / period);
    }
    return X;
  }
}
data {
  // ---- within-hut compositional data (one row per hut-night) ----
  int<lower=1> J;
  array[J, 4] int<lower=0> y;            // [ ,1]=alive&unfed [ ,2]=dead&unfed [ ,3]=dead&fed [ ,4]=alive&fed
  array[J] int<lower=0, upper=1> arm;    // 0 = control, 1 = treatment
  vector<lower=0>[J] time;               // hut-trial day, per hut-row (drives both decay and season)

  // ---- deterrence (entry) data, aggregated to trial-day level ----
  int<lower=1> K;                        // number of trial-day units with both arms present
  array[K] int<lower=0> Tcatch;          // total catch summed across treatment huts, trial-day k
  array[K] int<lower=0> Ccatch;          // total catch summed across control huts, trial-day k
  vector<lower=0>[K] time_det;           // time for trial-day k (same scale as `time` above)
  vector[K] log_hut_ratio;               // log(n_treatment_huts / n_control_huts) that day; 0 if balanced

  // ---- link from each hut-row to its trial-day deterrence group (for generated quantities) ----
  array[J] int<lower=1, upper=K> det_idx;

  real<lower=0> T_ref;                   // reference time (days) at which evidence indicates
                                          // effectively no residual effect -- e.g. 40.

  // ---- seasonal biting probability ----
  real<lower=0> period;                  // days per seasonal cycle (e.g. 365)
  int<lower=0> n_harmonics;              // number of Fourier harmonics H
  int<lower=1> N_grid;                   // grid of trial days on which to report pi(t)
  vector[N_grid] time_grid;
}
transformed data {
  matrix[J, 2 * n_harmonics] X_pi           = fourier_basis(time, period, n_harmonics);
  matrix[N_grid, 2 * n_harmonics] X_pi_grid = fourier_basis(time_grid, period, n_harmonics);
}
parameters {
  real alpha0;                           // logit-scale mean biting probability across the cycle
  vector[2 * n_harmonics] beta_pi;       // Fourier coefficients [cos_1, sin_1, cos_2, sin_2, ...]

  real<lower=0, upper=1> nu_u;           // background hazard, before average feeding time (both arms)
  real<lower=0, upper=1> nu_f;           // background hazard, after average feeding time (both arms)

  real<lower=0, upper=1> R0; real<lower=0, upper=1> eps_R; real<lower=0> nR;   // mu_R: repelled before feeding
  real<lower=0, upper=1> K0; real<lower=0, upper=1> eps_K; real<lower=1> nK;   // mu_K: killed before feeding
  real<lower=0, upper=1> M0; real<lower=0, upper=1> eps_M; real<lower=1> nM;   // mu_M: killed after feeding

  real<lower=0> phi;                     // shared Dirichlet-multinomial precision (within-hut composition)

  real<lower=0, upper=1> delta0; real<lower=0, upper=1> eps_d; real<lower=0> nd;
  real<lower=0, upper=1> rho_delta;      // overdispersion, ICC scale, for the deterrence beta-binomial
}
transformed parameters {
  real<lower=0> phi_delta = (1 - rho_delta) / rho_delta;
}
model {
  // ---- priors ----
  alpha0  ~ normal(0, 1.5);   // roughly flat on the probability scale, like the old pi ~ beta(1,1)
  // alpha0 ~ normal(logit(0.40), 0.2); // centre around 30% consistent with data
  beta_pi ~ normal(0, 1);     // seasonal amplitude: +/-1 on the logit scale is already a strong season
  // beta_pi ~ normal(0, 0.15);     // seasonal amplitude: if the baseline probability is 30%, adding a seasonal effect of +1 increases the probability to approximately 52%
  nu_u ~ beta(1, 1);
  nu_f ~ beta(1, 1);
  // nu_u ~ beta(1, 49);
  // nu_f ~ beta(1, 49);

  R0 ~ beta(2,4); nR ~ lognormal(log(2), 0.5);       // 1,1
  K0 ~ beta(2,4); nK ~ lognormal(log(2), 0.5);       // 3, 17 initial effect ~ 0.2 [beta(4, 16)] BUT [beta(3,17) gives 0.15] and gradual decay: nK > 1
  M0 ~ beta(2,4); nM ~ lognormal(log(2), 0.5);     // 2, 4 initial effect small ish and gradual decay: nK > 1
  delta0 ~ beta(2,4); nd ~ lognormal(log(2), 0.5); // 1,1

  // Strong priors: evidence indicates effectively no residual effect by T_ref.
  eps_R ~ beta(1, 49);
  eps_K ~ beta(1, 49);
  eps_M ~ beta(1, 49);
  eps_d ~ beta(1, 49);

  phi ~ gamma(2, 0.1);   // weakly informative; review against expected catch sizes
  rho_delta ~ beta(2, 5);

  // ---- within-hut compositional likelihood (Dirichlet-multinomial) ----
  vector[J] pi_t = inv_logit(alpha0 + X_pi * beta_pi);

  for (j in 1:J) {
    real muR_j; real muK_j; real muM_j;

    if (arm[j] == 1) {
      muR_j = hill_decay_eps(time[j], R0, eps_R, nR, T_ref);
      muK_j = hill_decay_eps(time[j], K0, eps_K, nK, T_ref);
      muM_j = hill_decay_eps(time[j], M0, eps_M, nM, T_ref);
    } else {
      muR_j = 0; muK_j = 0; muM_j = 0;    // structural constraint: mu_R = mu_K = mu_M = 0 in control
    }

    real A  = 1 - (1 - nu_u) * (1 - muK_j);   // combined hazard of dying before the average feeding time
    real hF = 1 - (1 - nu_f) * (1 - muM_j);   // combined hazard of dying after the average feeding time

    vector[4] p;
    p[1] = (1 - A) * (1 - hF) * ((1 - pi_t[j]) + pi_t[j] * muR_j);   // alive, unfed
    p[2] = A + (1 - A) * hF * ((1 - pi_t[j]) + pi_t[j] * muR_j);     // dead,  unfed
    p[3] = pi_t[j] * (1 - muR_j) * (1 - A) * hF;                      // dead,  fed
    p[4] = pi_t[j] * (1 - muR_j) * (1 - A) * (1 - hF);                // alive, fed

    y[j] ~ dirichlet_multinomial(p * phi);
  }

  // ---- deterrence (entry) likelihood, on matched trial-day totals ----
  for (k in 1:K) {
    real delta_k = hill_decay_eps(time_det[k], delta0, eps_d, nd, T_ref);
    real ratio_k = (1 - delta_k) * exp(log_hut_ratio[k]);       // expected treatment:control catch ratio
    real p_k     = ratio_k / (1 + ratio_k);                     // expected share of combined catch in treatment
    int total_k  = Tcatch[k] + Ccatch[k];
    Tcatch[k] ~ beta_binomial(total_k, p_k * phi_delta, (1 - p_k) * phi_delta);
  }
}
generated quantities {
  // seasonal biting probability on the reporting grid, and its cycle average
  vector[N_grid] pi_grid = inv_logit(alpha0 + X_pi_grid * beta_pi);
  real pi_mean = mean(pi_grid);

  vector[K] delta;
  for (k in 1:K) delta[k] = hill_decay_eps(time_det[k], delta0, eps_d, nd, T_ref);

  vector[J] p_unfed_alive_uncond;
  vector[J] p_unfed_dead_uncond;
  vector[J] p_fed_dead_uncond;
  vector[J] p_fed_alive_uncond;
  vector[J] p_deterred_uncond;

  {
    vector[J] pi_t = inv_logit(alpha0 + X_pi * beta_pi);

    for (j in 1:J) {
      real muR_j; real muK_j; real muM_j;
      if (arm[j] == 1) {
        muR_j = hill_decay_eps(time[j], R0, eps_R, nR, T_ref);
        muK_j = hill_decay_eps(time[j], K0, eps_K, nK, T_ref);
        muM_j = hill_decay_eps(time[j], M0, eps_M, nM, T_ref);
      } else {
        muR_j = 0; muK_j = 0; muM_j = 0;
      }

      real A  = 1 - (1 - nu_u) * (1 - muK_j);
      real hF = 1 - (1 - nu_f) * (1 - muM_j);
      real one_minus_delta = (arm[j] == 0) ? 1.0 : (1 - delta[det_idx[j]]);

      real p1 = (1 - A) * (1 - hF) * ((1 - pi_t[j]) + pi_t[j] * muR_j);
      real p2 = A + (1 - A) * hF * ((1 - pi_t[j]) + pi_t[j] * muR_j);
      real p3 = pi_t[j] * (1 - muR_j) * (1 - A) * hF;
      real p4 = pi_t[j] * (1 - muR_j) * (1 - A) * (1 - hF);

      p_unfed_alive_uncond[j] = one_minus_delta * p1;
      p_unfed_dead_uncond[j]  = one_minus_delta * p2;
      p_fed_dead_uncond[j]    = one_minus_delta * p3;
      p_fed_alive_uncond[j]   = one_minus_delta * p4;
      p_deterred_uncond[j]    = 1 - one_minus_delta;
    }
  }
}
