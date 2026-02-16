# High-level interface

#' Mortality model component specification
#'
#' @description
#' Create a mortality component specification for [vrcm()].
#'
#' @param formula A model formula for additional mortality covariates.
#' @param monotonic Optional character vector of level column names to treat as monotonic.
#' @param conflict How to model conflict in mortality ("fixed" or "region").
#' @param time How to model time in mortality ("national" or "region").
#'
#' @export
vrc_mortality <- function(
  formula = ~1,
  monotonic = NULL,
  conflict = c("fixed", "region"),
  time = c("national", "region")
) {
  if (!inherits(formula, "formula")) stop("formula must be a formula")
  structure(
    list(formula = formula, monotonic = monotonic, conflict = match.arg(conflict), time = match.arg(time)),
    class = "vrc_mortality"
  )
}

#' Reporting model component specification
#'
#' @description
#' Create a reporting component specification for [vrcm()].
#'
#' @param formula A model formula for additional reporting covariates.
#' @param monotonic Optional character vector of level column names to treat as monotonic.
#' @param conflict How to model conflict in reporting ("fixed" or "region").
#' @param time How to model time in reporting ("national" or "region").
#'
#' @export
vrc_reporting <- function(
  formula = ~1,
  monotonic = NULL,
  conflict = c("fixed", "region"),
  time = c("national", "region")
) {
  if (!inherits(formula, "formula")) stop("formula must be a formula")
  structure(
    list(formula = formula, monotonic = monotonic, conflict = match.arg(conflict), time = match.arg(time)),
    class = "vrc_reporting"
  )
}

#' Fit a VR mortality + reporting model
#'
#' @param mortality A [vrc_mortality()] specification.
#' @param reporting A [vrc_reporting()] specification.
#' @param data A data.frame in canonical VR long format.
#' @param t0 Conflict start time.
#' @param use_mar_labels Logical.
#' @param priors Optional vrc_priors object.
#' @param ... Passed to [vrc_fit()].
#'
#' @export
vrcm <- function(
  mortality = vrc_mortality(~1),
  reporting = vrc_reporting(~1),
  data,
  t0,
  use_mar_labels = FALSE,
  priors = NULL,
  ...
) {
  if (inherits(mortality, "formula")) mortality <- vrc_mortality(mortality)
  if (inherits(reporting, "formula")) reporting <- vrc_reporting(reporting)

  out <- vrc_fit(
    data = data, t0 = t0,
    mortality_covariates = mortality$formula,
    reporting_covariates = reporting$formula,
    mortality_monotonic = mortality$monotonic,
    reporting_monotonic = reporting$monotonic,
    mortality_conflict = mortality$conflict,
    reporting_conflict = reporting$conflict,
    mortality_time = mortality$time,
    reporting_time = reporting$time,
    use_mar_labels = use_mar_labels,
    priors = priors,
    ...
  )
  out$mortality <- mortality
  out$reporting <- reporting
  out
}

#' Summarise covariate effects
#'
#' @param x A `vrcfit` object.
#' @param probs Quantiles to include.
#' @param original_scale Logical.
#'
#' @export
vrc_coef_summary <- function(x, probs = c(0.1, 0.5, 0.9), original_scale = FALSE) {
  if (!inherits(x, "vrcfit")) stop("x must be a vrcfit object")
  
  K_mort <- x$standata$K_mort
  K_rep <- x$standata$K_rep
  K_mono_mort <- x$standata$K_mono_mort
  K_mono_rep <- x$standata$K_mono_rep

  pars <- c("beta_conf", "gamma_conf")
  if (K_mort > 0) pars <- c(pars, "beta_mort")
  if (K_rep > 0) pars <- c(pars, "gamma_rep")
  if (K_mono_mort > 0) pars <- c(pars, "beta_mono_mort_pars")
  if (K_mono_rep > 0) pars <- c(pars, "beta_mono_rep_pars")

  e <- rstan::extract(x$stanfit, pars = pars, permuted = TRUE)
  cause_levels <- x$meta$cause_levels

  summarise_mat <- function(mat) {
    if (is.null(dim(mat))) mat <- matrix(mat, ncol = 1)
    mean_ <- colMeans(mat)
    sd_ <- apply(mat, 2, stats::sd)
    qs <- t(apply(mat, 2, stats::quantile, probs = probs))
    out <- data.frame(mean = mean_, sd = sd_, qs, check.names = FALSE)
    colnames(out)[-(1:2)] <- paste0("q", gsub("\\.", "", format(probs, trim = TRUE)))
    out
  }

  rows <- list()

  if (!is.null(e$beta_conf)) {
    mat <- if (length(dim(e$beta_conf)) == 1) matrix(e$beta_conf, ncol = 1) else e$beta_conf
    sm <- summarise_mat(mat)
    rows[[length(rows) + 1]] <- data.frame(component = "mortality", parameter = "beta_conf", cause = cause_levels[seq_len(ncol(mat))], term = "conflict", sm)
  }

  if (K_mort > 0 && !is.null(e$beta_mort)) {
    terms <- x$scaling$X_mort$colnames
    for (g in seq_len(dim(e$beta_mort)[2])) {
      mat <- matrix(e$beta_mort[, g, ], ncol = K_mort)
      sm <- summarise_mat(mat)
      rows[[length(rows) + 1]] <- data.frame(component = "mortality", parameter = "beta_mort", cause = rep(cause_levels[g], K_mort), term = terms, sm)
    }
  }

  if (K_mono_mort > 0 && !is.null(e$beta_mono_mort_pars)) {
    terms <- x$mortality$monotonic
    for (g in seq_len(dim(e$beta_mono_mort_pars)[3])) {
      mat <- matrix(e$beta_mono_mort_pars[, , g], ncol = K_mono_mort)
      sm <- summarise_mat(mat)
      rows[[length(rows) + 1]] <- data.frame(component = "mortality", parameter = "beta_mono_mort", cause = rep(cause_levels[g], K_mono_mort), term = terms, sm)
    }
  }

  if (!is.null(e$gamma_conf)) {
    mat <- if (length(dim(e$gamma_conf)) == 1) matrix(e$gamma_conf, ncol = 1) else e$gamma_conf
    sm <- summarise_mat(mat)
    rows[[length(rows) + 1]] <- data.frame(component = "reporting", parameter = "gamma_conf", cause = cause_levels[seq_len(ncol(mat))], term = "conflict", sm)
  }

  if (K_rep > 0 && !is.null(e$gamma_rep)) {
    terms <- x$scaling$X_rep$colnames
    for (g in seq_len(dim(e$gamma_rep)[2])) {
      mat <- matrix(e$gamma_rep[, g, ], ncol = K_rep)
      sm <- summarise_mat(mat)
      rows[[length(rows) + 1]] <- data.frame(component = "reporting", parameter = "gamma_rep", cause = rep(cause_levels[g], K_rep), term = terms, sm)
    }
  }

  if (K_mono_rep > 0 && !is.null(e$beta_mono_rep_pars)) {
    terms <- x$reporting$monotonic
    for (g in seq_len(dim(e$beta_mono_rep_pars)[3])) {
      mat <- matrix(e$beta_mono_rep_pars[, , g], ncol = K_mono_rep)
      sm <- summarise_mat(mat)
      rows[[length(rows) + 1]] <- data.frame(component = "reporting", parameter = "beta_mono_rep", cause = rep(cause_levels[g], K_mono_rep), term = terms, sm)
    }
  }

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' @export
plot.vrcfit <- function(x, type = c("reporting", "mortality"), ...) {
  type <- match.arg(type)
  if (type == "reporting") return(plot_reporting(x, ...))
  plot_mortality(x, ...)
}

#' @export
fitted.vrcfit <- function(object, ...) {
  df <- posterior_expected_counts(object, draws = FALSE)
  df$mu_mean
}

#' @export
residuals.vrcfit <- function(object, ...) {
  mu <- stats::fitted(object); y <- object$data$y; y - mu
}

#' @export
coef.vrcfit <- function(object, ...) {
  sm <- vrc_coef_summary(object, probs = 0.5, original_scale = FALSE)
  nm <- paste0(sm$component, ":", sm$cause, ":", sm$term)
  stats::setNames(sm$mean, nm)
}
