/**
 * Compute monotonic coefficients for healthcare facility levels.
 * 
 * Level 1 is most functional, Level K is least functional (provided).
 * Baseline (non-functioning) is not provided and has effect 0.
 * 
 * Result: weight_1 = 1.0 > weight_2 > ... > weight_K > 0.0.
 *
 * For mortality: beta_fac = B_best * weights (B_best < 0)
 * beta_fac_1 < beta_fac_2 < ... < beta_fac_K < 0.
 * 
 * For reporting: gamma_fac = B_best * weights (B_best > 0)
 * gamma_fac_1 > gamma_fac_2 > ... > gamma_fac_K > 0.
 *
 * gap_ratios is a simplex of size K.
 */
vector compute_fac_beta(int K, real B_best, vector gap_ratios) {
  vector[K] beta;
  if (K > 0) {
    // weight_K = gap_ratios_K (gap between level K and baseline 0)
    // weight_K-1 = gap_ratios_K + gap_ratios_K-1
    // ...
    // weight_1 = sum(gap_ratios) = 1.0
    
    vector[K] cs = cumulative_sum(gap_ratios); 
    // weights are 1 - sum(gaps before)
    // weight_1 = 1.0
    // weight_2 = 1.0 - gaps_1
    // ...
    beta[1] = B_best;
    for (k in 2:K) {
      beta[k] = B_best * (1.0 - cs[k-1]);
    }
  }
  return beta;
}
