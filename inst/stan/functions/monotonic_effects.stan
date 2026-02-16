/**
 * Compute monotonic coefficients using the "Gaps" structure.
 * 
 * beta_K: The baseline (smallest coefficient).
 * B: The total positive range (beta_1 - beta_K).
 * gap_ratios: A simplex of size K-1 describing the distribution of gaps.
 * 
 * Result: beta_1 > beta_2 > ... > beta_K.
 * beta_i = beta_K + B * sum(gap_ratios[i : K-1])
 */
vector compute_mono_beta(int K, real beta_K, real B, vector gap_ratios) {
  vector[K] beta;
  beta[K] = beta_K;
  if (K > 1) {
    // Cumulative sum approach:
    // sum(gap_ratios[i:K-1]) = 1 - sum(gap_ratios[1:i-1])
    vector[K-1] cs = cumulative_sum(gap_ratios);
    beta[1] = beta_K + B;
    for (i in 2:(K-1)) {
      beta[i] = beta_K + B * (1.0 - cs[i-1]);
    }
  }
  return beta;
}
