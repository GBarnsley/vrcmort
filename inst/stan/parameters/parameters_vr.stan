// ----------------------------
  // Mortality process: log lambda
  // ----------------------------
  vector[G] alpha0;                    // cause-specific intercept
  matrix[A, G] alpha_age;              // age effects by cause
  matrix[S, G] alpha_sex;              // sex effects by cause

  matrix[R, G] u_lambda_raw;           // region RE (raw), will be centred and scaled
  array[G] vector[T] v_lambda_eps;     // time RW1 innovations (raw)

  vector<lower=0>[G] sigma_u_lambda;
  vector<lower=0>[G] sigma_v_lambda;

  vector<lower=0>[G] beta_conf;        // constrain >= 0 to avoid "conflict reduces mortality" pathology
  // Optional region-varying conflict effect (random slope by region)
  matrix[use_beta_conf_re == 1 ? R : 0, G] beta_conf_re_raw;
  vector<lower=0>[G] sigma_beta_conf;
  matrix[G, K_mort] beta_mort;         // additional mortality covariate effects

  // Healthcare facility effects (mortality)
  // beta_best < 0 (better facilities decrease mortality)
  vector<upper=0>[K_fac_mort > 0 ? G : 0] beta_fac_best_mort; 
  array[K_fac_mort > 0 ? G : 0] simplex[K_fac_mort > 0 ? K_fac_mort : 1] gap_ratios_fac_mort;

  // ----------------------------
  // Reporting process: logit rho
  // ----------------------------
  vector[G] kappa0;                    // baseline reporting (pre-conflict)
  vector[G] kappa_post;                // post-conflict shift

  matrix[R, G] u_rho_raw;              // region RE (raw)
  array[G] vector[T] v_rho_eps;        // time RW1 innovations (raw)

  vector<lower=0>[G] sigma_u_rho;
  vector<lower=0>[G] sigma_v_rho;

  vector[G] gamma_conf;                // conflict effect on reporting
  // Optional region-varying conflict effect (random slope by region)
  matrix[use_gamma_conf_re == 1 ? R : 0, G] gamma_conf_re_raw;
  vector<lower=0>[G] sigma_gamma_conf;
  matrix[G, K_rep] gamma_rep;          // additional reporting covariate effects

  // Healthcare facility effects (reporting)
  // gamma_best > 0 (better facilities increase reporting)
  vector<lower=0>[K_fac_rep > 0 ? G : 0] gamma_fac_best_rep;
  array[K_fac_rep > 0 ? G : 0] simplex[K_fac_rep > 0 ? K_fac_rep : 1] gap_ratios_fac_rep;

  // Optional region-specific time random walks (deviations around the national trend)
  array[G] matrix[use_rw_region_lambda == 1 ? T : 0, use_rw_region_lambda == 1 ? R : 0] v_lambda_region_eps;
  vector<lower=0>[G] sigma_v_lambda_region;

  array[G] matrix[use_rw_region_rho == 1 ? T : 0, use_rw_region_rho == 1 ? R : 0] v_rho_region_eps;
  vector<lower=0>[G] sigma_v_rho_region;

  // Age-selective post-conflict reporting drop for non-trauma:
  // build a monotone increasing "age penalty" then apply it negatively post-conflict.
  vector<lower=0>[A-1] delta_age_incr; // increments between age groups (positive)
  real<lower=0> delta_age_scale;

  // Observation dispersion (NB2)
  vector<lower=0>[G] phi;

  // Labeling probability (MAR)
  real<lower=0, upper=1> omega;
