// ----------------------------
  // Likelihood
  // ----------------------------
  if (prior_PD == 0) {
    // 1. Calculate Log Lambda (Vectorized)
    vector[M] log_lambda;
    for (m in 1:M) {
      int g = idx_g[m];
      int t = idx_t[m];
      int a = idx_a[m];
      int s = idx_s[m];
      int r = idx_r[m];

      real mort_x = 0;
      if (K_mort > 0) mort_x = X_mort[m] * (beta_mort[g]');
      if (K_fac_mort > 0) mort_x += X_fac_mort[m] * beta_fac_mort[, g];

      log_lambda[m] = alpha0[g]
                    + alpha_age[a, g]
                    + alpha_sex[s, g]
                    + u_lambda[r, g]
                    + v_lambda[g][t]
                    + beta_conf_rg[r, g] * conflict_flat[m]
                    + mort_x;
      
      if (use_rw_region_lambda == 1) log_lambda[m] += v_lambda_region[g][t, r];
    }

    // 2. Calculate Logit Rho (Vectorized)
    vector[M] logit_rho;
    for (m in 1:M) {
      int g = idx_g[m];
      int t = idx_t[m];
      int a = idx_a[m];
      int r = idx_r[m];

      real rep_x = 0;
      if (K_rep > 0) rep_x = X_rep[m] * (gamma_rep[g]');
      if (K_fac_rep > 0) rep_x += X_fac_rep[m] * gamma_fac_rep[, g];

      logit_rho[m] = kappa0[g]
                   + kappa_post[g] * post[t]
                   + u_rho[r, g]
                   + v_rho[g][t]
                   + gamma_conf_rg[r, g] * conflict_flat[m]
                   + rep_x;
      
      if (use_rw_region_rho == 1) logit_rho[m] += v_rho_region[g][t, r];

      if (g == cause_idx_age_penalty) {
        logit_rho[m] += - delta_age[a] * post[t];
      }
    }

    // 3. Combine to log_mus
    vector[M] log_mus_flat = log_exposure_flat + log_lambda + log_inv_logit(logit_rho);
    
    // 4. Sampling Statement (Vectorized)
    // We must still handle cause-specific phi and labeling
    for (j in 1:N) {
      int g = cause[j];
      int i_start = (j - 1) * R + 1;
      int i_end = j * R;
      
      if (use_mar_labels == 1) {
        y[j] ~ neg_binomial_2_log(to_array_1d(log_mus_flat[i_start:i_end] + log(omega)), phi[g]);
      } else {
        y[j] ~ neg_binomial_2_log(to_array_1d(log_mus_flat[i_start:i_end]), phi[g]);
      }

      if (use_mar_labels == 1) {
        target += neg_binomial_2_sum_fast_log_lpmf(y_miss[j] | log_mus_flat[i_start:i_end] + log1m(omega), phi[g]);
      }
    }
  }
