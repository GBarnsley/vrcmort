#' Build Stan data for the VR reporting model
#'
#' @description
#' Convert a long-format VR dataset into the list structure expected by the
#' Stan model shipped with the package.
#'
#' @param data A data.frame in VR long format.
#' @param t0 Conflict start time.
#' @param mortality_covariates Optional formula for additional mortality covariates.
#' @param reporting_covariates Optional formula for additional reporting covariates.
#' @param mortality_facilities Optional character vector of level column names for mortality.
#'   Provided in decreasing order of functionality (perfect -> minimal).
#' @param reporting_facilities Optional character vector of level column names for reporting.
#'   Provided in decreasing order of functionality (perfect -> minimal).
#' @param mortality_conflict How to model conflict in mortality ("fixed" or "region").
#' @param reporting_conflict How to model conflict in reporting ("fixed" or "region").
#' @param mortality_time How to model time in mortality ("national" or "region").
#' @param reporting_time How to model time in reporting ("national" or "region").
#' @param standardise Logical. If TRUE, standardise non-facility covariates.
#' @param scale_binary Logical.
#' @param drop_na_y Logical.
#' @param duplicates Handling of duplicates.
#' @param use_mar_labels Logical.
#' @param priors Optional vrc_priors object.
#' @param prior_PD Logical.
#'
#' @export
vrc_standata <- function(
  data,
  t0,
  mortality_covariates = NULL,
  reporting_covariates = NULL,
  mortality_facilities = NULL,
  reporting_facilities = NULL,
  mortality_conflict = c("fixed", "region"),
  reporting_conflict = c("fixed", "region"),
  mortality_time = c("national", "region"),
  reporting_time = c("national", "region"),
  standardise = TRUE,
  scale_binary = FALSE,
  drop_na_y = TRUE,
  duplicates = c("error", "sum"),
  use_mar_labels = FALSE,
  priors = NULL,
  prior_PD = FALSE
) {
  stopifnot(is.data.frame(data))
  duplicates <- match.arg(duplicates)
  mortality_conflict <- match.arg(mortality_conflict)
  reporting_conflict <- match.arg(reporting_conflict)
  mortality_time <- match.arg(mortality_time)
  reporting_time <- match.arg(reporting_time)

  idx <- vrc_index(
    data = data,
    t0 = t0,
    duplicates = duplicates,
    sort_time = TRUE
  )
  df <- idx$data
  df_miss <- idx$data_miss
  meta <- idx$meta

  meta$model_options <- list(
    mortality_conflict = mortality_conflict,
    reporting_conflict = reporting_conflict,
    mortality_time = mortality_time,
    reporting_time = reporting_time,
    use_mar_labels = use_mar_labels,
    mortality_facilities = mortality_facilities,
    reporting_facilities = reporting_facilities
  )

  if (isTRUE(drop_na_y)) {
    df <- df[!is.na(df$y), , drop = FALSE]
  }

  check_positive(df$exposure, "exposure")

  if (standardise) {
    sc <- standardise_vector(df$conflict)
    df$conflict_z <- sc$x
    conflict_scaling <- sc
  } else {
    df$conflict_z <- as.numeric(df$conflict)
    conflict_scaling <- list(x = NULL, centre = 0, scale = 1)
  }

  R <- meta$R
  T <- meta$T
  G <- meta$G
  join_cols <- c("time_id", "age_id", "sex_id", "cause_id")

  groups_in_df <- df |> dplyr::distinct(dplyr::across(dplyr::all_of(join_cols)))
  groups_in_miss <- df_miss |>
    dplyr::distinct(dplyr::across(dplyr::all_of(join_cols)))

  if (use_mar_labels && nrow(groups_in_miss) < nrow(groups_in_df)) {
    warning("Some cells with labeled deaths have no corresponding missing-region entry. They will be treated as having zero missing-region deaths.", call. = FALSE)
  }

  all_groups <- dplyr::union(groups_in_df, groups_in_miss) |>
    dplyr::arrange(.data$time_id, .data$age_id, .data$sex_id, .data$cause_id)
  N_groups <- nrow(all_groups)

  if (!use_mar_labels || nrow(df_miss) == 0L) {
    df_miss_final <- all_groups
    df_miss_final$y <- 0
  } else {
    df_miss_final <- all_groups |> dplyr::left_join(df_miss, by = join_cols)
    df_miss_final$y[is.na(df_miss_final$y)] <- 0
  }

  grid <- expand.grid(region_id = seq_len(R), group_id = seq_len(N_groups))
  df_complete <- grid |>
    dplyr::inner_join(
      all_groups |> dplyr::mutate(group_id = dplyr::row_number()),
      by = "group_id"
    ) |>
    dplyr::left_join(df, by = c(join_cols, "region_id")) |>
    dplyr::arrange(.data$group_id, .data$region_id)

  df_complete$y[is.na(df_complete$y)] <- 0
  df_complete$exposure[is.na(df_complete$exposure)] <- 1e-9
  df_complete$conflict_z[is.na(df_complete$conflict_z)] <- 0

  X_mort <- model_matrix_no_intercept(mortality_covariates, df_complete)
  X_rep <- model_matrix_no_intercept(reporting_covariates, df_complete)

  if (standardise) {
    X_mort <- standardise_matrix(X_mort, scale_binary = scale_binary)
    X_rep <- standardise_matrix(X_rep, scale_binary = scale_binary)
  }

  get_fac <- function(cols, data) {
    if (is.null(cols)) {
      return(matrix(0, nrow(data), 0))
    }
    stop_if_missing_cols(data, cols)
    out <- as.matrix(data[, cols, drop = FALSE])
    out[is.na(out)] <- 0
    out
  }

  X_fac_mort <- get_fac(mortality_facilities, df_complete)
  X_fac_rep <- get_fac(reporting_facilities, df_complete)

  post <- as.integer(seq_len(T) >= meta$t0)

  standata <- list(
    N = N_groups,
    R = R,
    T = T,
    A = meta$A,
    S = meta$S,
    G = G,
    time = as_stan_array_int(all_groups$time_id),
    age = as_stan_array_int(all_groups$age_id),
    sex = as_stan_array_int(all_groups$sex_id),
    cause = as_stan_array_int(all_groups$cause_id),
    y = matrix(df_complete$y, nrow = N_groups, ncol = R, byrow = TRUE),
    exposure = matrix(
      df_complete$exposure,
      nrow = N_groups,
      ncol = R,
      byrow = TRUE
    ),
    conflict = matrix(
      df_complete$conflict_z,
      nrow = N_groups,
      ncol = R,
      byrow = TRUE
    ),
    use_beta_conf_re = as.integer(mortality_conflict == "region"),
    use_gamma_conf_re = as.integer(reporting_conflict == "region"),
    use_rw_region_lambda = as.integer(mortality_time == "region"),
    use_rw_region_rho = as.integer(reporting_time == "region"),
    K_mort = ncol(X_mort),
    X_mort = X_mort,
    K_rep = ncol(X_rep),
    X_rep = X_rep,
    K_fac_mort = ncol(X_fac_mort),
    X_fac_mort = X_fac_mort,
    K_fac_rep = ncol(X_fac_rep),
    X_fac_rep = X_fac_rep,
    post = as_stan_array_int(post),
    t0 = as.integer(meta$t0),
    y_miss = as_stan_array_int(df_miss_final$y),
    use_mar_labels = as.integer(use_mar_labels),
    prior_PD = as.integer(prior_PD)
  )

  priors_resolved <- vrc_resolve_priors(
    priors,
    G,
    ncol(X_mort),
    ncol(X_rep),
    X_mort,
    X_rep
  )
  standata <- c(standata, priors_resolved)

  list(
    standata = standata,
    df = df_complete,
    meta = meta,
    scaling = list(
      conflict = conflict_scaling,
      X_mort = attr(X_mort, "scaling"),
      X_rep = attr(X_rep, "scaling")
    ),
    priors = if (is.null(priors)) vrc_priors() else priors,
    priors_resolved = priors_resolved
  )
}
