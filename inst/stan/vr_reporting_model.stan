// Base VR mortality + reporting model (assembled from include files)
//
// This file mirrors the modular Stan layout used by the epidemia R package.
// Individual pieces live under inst/stan/{functions,data,parameters,tparameters,model,generated_quantities}.

#include "include/license.stan"

functions {
#include "functions/rw1_centered.stan"
#include "functions/neg_binomial_2_sum.stan"
#include "functions/facility_effects.stan"
}

data {
#include "data/data_vr.stan"
}

transformed data {
  matrix[N, R] log_exposure = log(exposure);
  
  int M = N * R;
  array[M] int y_flat;
  array[M] int idx_g;
  array[M] int idx_t;
  array[M] int idx_a;
  array[M] int idx_s;
  array[M] int idx_r;
  vector[M] log_exposure_flat;
  vector[M] conflict_flat;

  for (j in 1:N) {
    for (r in 1:R) {
      int pos = (j - 1) * R + r;
      y_flat[pos] = y[j, r];
      idx_g[pos] = cause[j];
      idx_t[pos] = time[j];
      idx_a[pos] = age[j];
      idx_s[pos] = sex[j];
      idx_r[pos] = r;
      log_exposure_flat[pos] = log_exposure[j, r];
      conflict_flat[pos] = conflict[j, r];
    }
  }
}

parameters {
#include "parameters/parameters_vr.stan"
}

transformed parameters {
#include "tparameters/tparameters_vr.stan"
}

model {
#include "model/priors_vr.stan"
#include "model/likelihood_vr_nb2.stan"
}

generated quantities {
#include "generated_quantities/gq_vr.stan"
}
