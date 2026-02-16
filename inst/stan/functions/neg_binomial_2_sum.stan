  /**
   * Log PMF of the sum of independent Negative Binomial 2 variables
   * using an exact recurrence relation.
   *
   * @param y Observed sum
   * @param mus Vector of means for each component
   * @param phi Dispersion parameter (shared)
   * @return Log probability mass
   */
  real neg_binomial_2_sum_lpmf(int y, vector mus, real phi) {
    int R = num_elements(mus);
    
    // Base case: y = 0 (Probability all components are 0)
    if (y == 0) {
      real log_p0 = 0;
      for (r in 1:R) {
        log_p0 += neg_binomial_2_lpmf(0 | mus[r], phi);
      }
      return log_p0;
    }

    // Pre-calculate log(p) and log(1-p) for each component
    // In Stan's neg_binomial_2: p = phi / (phi + mu)
    vector[R] log_one_minus_p; 
    for (r in 1:R) {
      log_one_minus_p[r] = log(mus[r]) - log(mus[r] + phi);
    }

    // Recurrence array
    vector[y + 1] lp; 
    lp[1] = 0;
    for (r in 1:R) {
      lp[1] += neg_binomial_2_lpmf(0 | mus[r], phi);
    }

    // C[j] = phi * sum( (1-p_r)^j )
    // We vectorize the calculation of log_C to improve gradient speed
    vector[y] log_C;
    real log_phi = log(phi);
    for (j in 1:y) {
      log_C[j] = log_phi + log_sum_exp(j * log_one_minus_p);
    }

    // Main Recurrence: lp[k+1] is log(P(Sum = k))
    for (k in 1:y) {
      // Log-space convolution of previous probabilities and C coefficients
      vector[k] components;
      for (j in 1:k) {
        components[j] = lp[k - j + 1] + log_C[j];
      }
      lp[k + 1] = log_sum_exp(components) - log(k);
    }

    return lp[y + 1];
  }

  /**
   * Log PMF of the sum of independent Negative Binomial 2 variables
   * using an exact recurrence relation, with means provided on log scale.
   */
  real neg_binomial_2_sum_log_lpmf(int y, vector log_mus, real phi) {
    int R = num_elements(log_mus);
    
    if (y == 0) {
      real log_p0 = 0;
      for (r in 1:R) {
        log_p0 += neg_binomial_2_log_lpmf(0 | log_mus[r], phi);
      }
      return log_p0;
    }

    // p = phi / (phi + exp(log_mu))
    // 1-p = exp(log_mu) / (phi + exp(log_mu))
    // log(1-p) = log_mu - log(phi + exp(log_mu))
    vector[R] log_one_minus_p; 
    for (r in 1:R) {
      log_one_minus_p[r] = log_mus[r] - log_sum_exp(log(phi), log_mus[r]);
    }

    vector[y + 1] lp; 
    lp[1] = 0;
    for (r in 1:R) {
      lp[1] += neg_binomial_2_log_lpmf(0 | log_mus[r], phi);
    }

    vector[y] log_C;
    real log_phi = log(phi);
    for (j in 1:y) {
      log_C[j] = log_phi + log_sum_exp(j * log_one_minus_p);
    }

    for (k in 1:y) {
      vector[k] components;
      for (j in 1:k) {
        components[j] = lp[k - j + 1] + log_C[j];
      }
      lp[k + 1] = log_sum_exp(components) - log(k);
    }

    return lp[y + 1];
  }

  /**
   * Log PMF of the sum of independent Negative Binomial 2 variables
   * using an moment matching, should be faster than the full method
   *
   * @param y Observed sum
   * @param mus Vector of means for each component
   * @param phi Dispersion parameter (shared)
   * @return Log probability mass
   */
  real neg_binomial_2_sum_fast_lpmf(int y, vector mus, real phi) {

    real mus_sum = 0;
    real mus_var_diff = 0;

    int R = num_elements(mus);

    for (r in 1:R) {
      mus_sum += mus[r];
      mus_var_diff += square(mus[r]) / phi;
    }

    real phi_eff = square(mus_sum) / mus_var_diff;

    return neg_binomial_2_lpmf(y | mus_sum, phi_eff);
  }

  /**
   * Log PMF of the sum of independent Negative Binomial 2 variables
   * using moment matching, with means provided on log scale.
   */
  real neg_binomial_2_sum_fast_log_lpmf(int y, vector log_mus, real phi) {
    real log_mus_sum = log_sum_exp(log_mus);
    real log_mus_var_diff = log_sum_exp(2 * log_mus) - log(phi);
    real phi_eff = exp(2 * log_mus_sum - log_mus_var_diff);

    return neg_binomial_2_log_lpmf(y | log_mus_sum, phi_eff);
  }

  /**
   * Safe RNG for Negative Binomial 2.
   * Caps the mean to prevent Gamma draw overflow (Stan limit ~1e9).
   */
  int safe_neg_binomial_2_rng(real mu, real phi) {
    if (mu > 1e8) {
       // If mu is huge, the exact draw doesn't matter for validation, 
       // but we must not crash. Return a large safe value.
       return 100000000; 
    }
    if (mu <= 1e-15) return 0;
    
    // Also guard against extremely small phi which makes variance explode
    real safe_phi = fmax(phi, 1e-5);
    return neg_binomial_2_rng(mu, safe_phi);
  }

  /**
   * Safe RNG for Negative Binomial 2, mean on log scale.
   */
  int safe_neg_binomial_2_log_rng(real log_mu, real phi) {
    if (log_mu > 18.42) { // log(1e8) approx 18.42
       return 100000000; 
    }
    if (log_mu < -34.5) return 0; // log(1e-15) approx -34.5
    
    real safe_phi = fmax(phi, 1e-5);
    return neg_binomial_2_rng(exp(log_mu), safe_phi);
  }
