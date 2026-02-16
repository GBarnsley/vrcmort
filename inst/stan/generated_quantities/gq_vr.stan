vector[N * R] mu_rep;
  vector[N * R] rho_rep;
  vector[N * R] lambda_rep;
  vector[N * R] log_lik;
  array[N * R] int y_rep;

  vector[N] log_lik_miss;

  // Final coefficients for extraction
  matrix[K_fac_mort, G] beta_fac_mort_pars = beta_fac_mort;
  matrix[K_fac_rep, G] gamma_fac_rep_pars = gamma_fac_rep;

  for (j in 1:N) {
    int g = cause[j];
    int t = time[j];
    int a = age[j];
    int s = sex[j];

    vector[R] mus;
    for (r in 1:R) {
      int i = (j - 1) * R + r;

      real mort_x = 0;
      real rep_x = 0;
      if (K_mort > 0) mort_x = X_mort[i] * (beta_mort[g]');
      if (K_rep > 0)  rep_x  = X_rep[i]  * (gamma_rep[g]');

      // Facility contributions (monotonic)
      if (K_fac_mort > 0) mort_x += X_fac_mort[i] * beta_fac_mort[, g];
      if (K_fac_rep > 0)  rep_x  += X_fac_rep[i]  * gamma_fac_rep[, g];

      real log_lambda_r = alpha0[g]
                        + alpha_age[a, g]
                        + alpha_sex[s, g]
                        + u_lambda[r, g]
                        + v_lambda[g][t]
                        + beta_conf_rg[r, g] * conflict[j, r]
                        + mort_x;
      
      if (use_rw_region_lambda == 1) log_lambda_r += v_lambda_region[g][t, r];

      real logit_rho_r = kappa0[g]
                       + kappa_post[g] * post[t]
                       + u_rho[r, g]
                       + v_rho[g][t]
                       + gamma_conf_rg[r, g] * conflict[j, r]
                       + rep_x;
      
      if (use_rw_region_rho == 1) logit_rho_r += v_rho_region[g][t, r];

      if (g == 2) {
        logit_rho_r += - delta_age[a] * post[t];
      }

      lambda_rep[i] = exp(log_lambda_r);
      rho_rep[i] = inv_logit(logit_rho_r);
      
      real log_mu = log_exposure[j, r] + log_lambda_r + log_inv_logit(logit_rho_r);
      real log_mu_rep = (use_mar_labels == 1) ? (log_mu + log(omega)) : log_mu;
      
      mu_rep[i] = exp(log_mu_rep);

      // Posterior predictive draw
      y_rep[i] = safe_neg_binomial_2_log_rng(log_mu_rep, phi[g]);

      // Pointwise log-likelihood (for LOO / WAIC)
      log_lik[i] = neg_binomial_2_log_lpmf(y[j, r] | log_mu_rep, phi[g]);
      
      mus[r] = log_mu;
    }

    if (use_mar_labels == 1) {
      log_lik_miss[j] = neg_binomial_2_sum_fast_log_lpmf(y_miss[j] | mus + log1m(omega), phi[g]);
    } else {
      log_lik_miss[j] = 0;
    }
  }
