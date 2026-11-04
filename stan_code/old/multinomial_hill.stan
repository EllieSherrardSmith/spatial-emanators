// Dirichlet-multinomial version of nested_hill_model.stan -- same tree, same
// Hill-decay hazards mu_R(t), mu_K(t), mu_M(t) and Hill-decay deterrence, but
// compositional overdispersion is carried by a single shared Dirichlet-multinomial
// concentration phi instead of per-row Beta-mixing on mu_R/mu_K/mu_M. Use this as
// a robustness check against nested_hill_model.stan: if the two give materially
// different answers, that itself is informative about which overdispersion
// structure the data actually support.
//
//   nu_u = background mortality hazard, BEFORE the average feeding time (both arms)
//   nu_f = background mortality hazard, AFTER the average feeding time (both arms)
//   mu_R, mu_K, mu_M = intervention-attributable hazards (repelled before feeding; killed
//           before feeding; killed after feeding), each a Hill/log-logistic decay in time,
//           structurally zero in control (arm == 0). Here mu_R(t)/mu_K(t)/mu_M(t) are the
//           fixed-effects Hill-decay values directly -- no per-row latent Beta-mixing,
//           since the Dirichlet-multinomial's own phi already absorbs compositional
//           overdispersion at the hut-row level.
//
// Outcome probabilities (identical derivation to the nested model):
//   A  = 1 - (1-nu_u)*(1-mu_K)          combined hazard of dying before the average feeding time
//   hF = 1 - (1-nu_f)*(1-mu_M)          combined hazard of dying after the average feeding time
//   p[alive,unfed] = (1-A)*(1-hF) * ((1-pi) + pi*mu_R)
//   p[dead, unfed] = A + (1-A)*hF * ((1-pi) + pi*mu_R)
//   p[dead, fed]   = pi*(1-mu_R)*(1-A)*hF
//   p[alive,fed]   = pi*(1-mu_R)*(1-A)*(1-hF)
// Setting mu_R=mu_K=mu_M=0 (control) recovers the control-only equations exactly.

functions {
  // Hill / log-logistic decay: f0 at t=0, half-life thalf, shape n, -> 0 as t -> infinity.
  real hill_decay(real t, real f0, real thalf, real n) {
    real tn    = pow(fmax(t, 0), n);
    real t50n  = pow(thalf, n);
    return f0 * t50n / (t50n + tn);
  }
}
data {
  // ---- within-hut compositional data (one row per hut-night) ----
  int<lower=1> J;
  array[J, 4] int<lower=0> y;            // [ ,1]=alive&unfed [ ,2]=dead&unfed [ ,3]=dead&fed [ ,4]=alive&fed
  array[J] int<lower=0, upper=1> arm;    // 0 = control, 1 = treatment
  vector<lower=0>[J] time;               // time since deployment, per hut-row

  // ---- deterrence (entry) data, aggregated to trial-day level ----
  int<lower=1> K;                        // number of trial-day units with both arms present
  array[K] int<lower=0> Tcatch;          // total catch summed across treatment huts, trial-day k
  array[K] int<lower=0> Ccatch;          // total catch summed across control huts, trial-day k
  vector<lower=0>[K] time_det;           // time for trial-day k (same scale as `time` above)
  vector[K] log_hut_ratio;               // log(n_treatment_huts / n_control_huts) that day; 0 if balanced

  // ---- link from each hut-row to its trial-day deterrence group (for generated quantities) ----
  array[J] int<lower=1, upper=K> det_idx;

  real<lower=0> pi_prior_a;              // Beta prior shape params for pi
  real<lower=0> pi_prior_b;
}
parameters {
  real<lower=0, upper=1> pi;

  real<lower=0, upper=1> nu_u;           // background hazard, before average feeding time (both arms)
  real<lower=0, upper=1> nu_f;           // background hazard, after average feeding time (both arms)

  // Treatment-arm-only intervention hazards, each a Hill decay from a peak at t=0 towards 0.
  real<lower=0, upper=1> R0; real<lower=0> thalf_R; real<lower=0> nR;   // mu_R: repelled before feeding
  real<lower=0, upper=1> K0; real<lower=0> thalf_K; real<lower=0> nK;   // mu_K: killed before feeding
  real<lower=0, upper=1> M0; real<lower=0> thalf_M; real<lower=0> nM;   // mu_M: killed after feeding

  real<lower=0> phi;                     // shared Dirichlet-multinomial precision (within-hut composition)

  // Deterrence: Hill decay, control arm is the delta == 0 reference by construction
  real<lower=0, upper=1> delta0; real<lower=0> thalf_d; real<lower=0> nd;
  real<lower=0, upper=1> rho_delta;      // overdispersion, ICC scale, for the deterrence beta-binomial
}
transformed parameters {
  real<lower=0> phi_delta = (1 - rho_delta) / rho_delta;
}
model {
  // ---- priors ----
  pi   ~ beta(pi_prior_a, pi_prior_b);
  nu_u ~ beta(1, 1);
  nu_f ~ beta(1, 1);

  R0 ~ beta(1, 1); thalf_R ~ normal(20, 20); nR ~ lognormal(log(2), 0.7);
  K0 ~ beta(1, 1); thalf_K ~ normal(20, 20); nK ~ lognormal(log(2), 0.7);
  M0 ~ beta(1, 1); thalf_M ~ normal(20, 20); nM ~ lognormal(log(2), 0.7);

  phi ~ gamma(2, 0.1);   // weakly informative; review against expected catch sizes

  delta0 ~ beta(1, 1); thalf_d ~ normal(20, 20); nd ~ lognormal(log(2), 0.7);
  rho_delta ~ beta(2, 5);

  // ---- within-hut compositional likelihood (Dirichlet-multinomial) ----
  for (j in 1:J) {
    real muR_j; real muK_j; real muM_j;

    if (arm[j] == 1) {
      muR_j = hill_decay(time[j], R0, thalf_R, nR);
      muK_j = hill_decay(time[j], K0, thalf_K, nK);
      muM_j = hill_decay(time[j], M0, thalf_M, nM);
    } else {
      muR_j = 0; muK_j = 0; muM_j = 0;    // structural constraint: mu_R = mu_K = mu_M = 0 in control
    }

    real A  = 1 - (1 - nu_u) * (1 - muK_j);   // combined hazard of dying before the average feeding time
    real hF = 1 - (1 - nu_f) * (1 - muM_j);   // combined hazard of dying after the average feeding time

    vector[4] p;
    p[1] = (1 - A) * (1 - hF) * ((1 - pi) + pi * muR_j);          // alive, unfed
    p[2] = A + (1 - A) * hF * ((1 - pi) + pi * muR_j);            // dead,  unfed
    p[3] = pi * (1 - muR_j) * (1 - A) * hF;                       // dead,  fed
    p[4] = pi * (1 - muR_j) * (1 - A) * (1 - hF);                 // alive, fed

    y[j] ~ dirichlet_multinomial(p * phi);
  }

  // ---- deterrence (entry) likelihood, on matched trial-day totals ----
  for (k in 1:K) {
    real delta_k = hill_decay(time_det[k], delta0, thalf_d, nd);
    real ratio_k = (1 - delta_k) * exp(log_hut_ratio[k]);       // expected treatment:control catch ratio
    real p_k     = ratio_k / (1 + ratio_k);                     // expected share of combined catch in treatment
    int total_k  = Tcatch[k] + Ccatch[k];
    Tcatch[k] ~ beta_binomial(total_k, p_k * phi_delta, (1 - p_k) * phi_delta);
  }
}
generated quantities {
  vector[K] delta;
  for (k in 1:K) delta[k] = hill_decay(time_det[k], delta0, thalf_d, nd);

  // Fully unconditional-on-entry outcome probabilities per hut-row, folding in deterrence.
  vector[J] p_unfed_alive_uncond;
  vector[J] p_unfed_dead_uncond;
  vector[J] p_fed_dead_uncond;
  vector[J] p_fed_alive_uncond;
  vector[J] p_deterred_uncond;

  for (j in 1:J) {
    real muR_j; real muK_j; real muM_j;
    if (arm[j] == 1) {
      muR_j = hill_decay(time[j], R0, thalf_R, nR);
      muK_j = hill_decay(time[j], K0, thalf_K, nK);
      muM_j = hill_decay(time[j], M0, thalf_M, nM);
    } else {
      muR_j = 0; muK_j = 0; muM_j = 0;
    }

    real A  = 1 - (1 - nu_u) * (1 - muK_j);
    real hF = 1 - (1 - nu_f) * (1 - muM_j);
    real one_minus_delta = (arm[j] == 0) ? 1.0 : (1 - delta[det_idx[j]]);

    real p1 = (1 - A) * (1 - hF) * ((1 - pi) + pi * muR_j);
    real p2 = A + (1 - A) * hF * ((1 - pi) + pi * muR_j);
    real p3 = pi * (1 - muR_j) * (1 - A) * hF;
    real p4 = pi * (1 - muR_j) * (1 - A) * (1 - hF);

    p_unfed_alive_uncond[j] = one_minus_delta * p1;
    p_unfed_dead_uncond[j]  = one_minus_delta * p2;
    p_fed_dead_uncond[j]    = one_minus_delta * p3;
    p_fed_alive_uncond[j]   = one_minus_delta * p4;
    p_deterred_uncond[j]    = 1 - one_minus_delta;
  }
}
