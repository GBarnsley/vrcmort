#' Fit the VR reporting model
#'
#' @description
#' Fit the hierarchical VR mortality model with facility functionality effects.
#'
#' @param data A data.frame in long format.
#' @param t0 Conflict start time.
#' @param mortality_covariates Optional formula for additional mortality covariates.
#' @param reporting_covariates Optional formula for additional reporting covariates.
#' @param mortality_facilities Optional character vector of facility level column names for mortality.
#' @param reporting_facilities Optional character vector of facility level column names for reporting.
#' @param mortality_conflict How to model conflict in mortality ("fixed" or "region").
#' @param reporting_conflict How to model conflict in reporting ("fixed" or "region").
#' @param mortality_time How to model time in mortality ("national" or "region").
#' @param reporting_time How to model time in reporting ("national" or "region").
#' @param standardise Logical.
#' @param scale_binary Logical.
#' @param drop_na_y Logical.
#' @param duplicates Duplicate handling.
#' @param use_mar_labels Logical.
#' @param algorithm Inference algorithm.
#' @param priors Optional prior bundle.
#' @param prior_PD Logical.
#' @param backend Backend.
#' @param stan_model Model name.
#' @param ... Passed to sampling.
#'
#' @export
vrc_fit <- function(
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
  algorithm = c("sampling", "meanfield", "fullrank"),
  priors = NULL,
  prior_PD = FALSE,
  backend = c("rstan"),
  stan_model = "vr_reporting_model",
  ...
) {
  call <- match.call(expand.dots = TRUE)

  sdat_obj <- vrc_standata(
    data = data,
    t0 = t0,
    mortality_covariates = mortality_covariates,
    reporting_covariates = reporting_covariates,
    mortality_facilities = mortality_facilities,
    reporting_facilities = reporting_facilities,
    mortality_conflict = match.arg(mortality_conflict),
    reporting_conflict = match.arg(reporting_conflict),
    mortality_time = match.arg(mortality_time),
    reporting_time = match.arg(reporting_time),
    standardise = standardise,
    scale_binary = scale_binary,
    drop_na_y = drop_na_y,
    duplicates = match.arg(duplicates),
    use_mar_labels = use_mar_labels,
    priors = priors,
    prior_PD = prior_PD
  )

  dots <- list(...)
  if (!is.null(dots$chains) && identical(dots$chains, 0)) {
    return(sdat_obj)
  }

  spec <- vrc_model_spec(model = stan_model, backend = match.arg(backend))
  sm <- vrc_model(spec)

  if (is.null(dots$init_r)) {
    dots$init_r <- 1e-6
  }

  args <- c(dots, list(object = sm, data = sdat_obj$standata))

  fit <- if (match.arg(algorithm) == "sampling") {
    do.call(rstan::sampling, args)
  } else {
    args$algorithm <- match.arg(algorithm)
    do.call(rstan::vb, args)
  }

  out <- list(
    call = call,
    stanfit = fit,
    standata = sdat_obj$standata,
    data = sdat_obj$df,
    meta = sdat_obj$meta,
    scaling = sdat_obj$scaling,
    priors = sdat_obj$priors,
    priors_resolved = sdat_obj$priors_resolved,
    algorithm = match.arg(algorithm),
    backend = match.arg(backend),
    stan_model = spec$name,
    stan_file = spec$file,
    model_spec = spec
  )
  class(out) <- "vrcfit"
  out
}

#' @export
print.vrcfit <- function(x, ...) {
  cat("vrcmort model fit\n")
  cat("- Stan model:", x$stan_model, "\n")
  cat(
    "- Dimensions: R=",
    x$standata$R,
    ", T=",
    x$standata$T,
    ", G=",
    x$standata$G,
    "\n"
  )
  invisible(x)
}

#' @export
summary.vrcfit <- function(
  object,
  pars = c("beta_conf", "kappa0", "kappa_post", "gamma_conf", "phi"),
  ...
) {
  rstan::summary(object$stanfit, pars = pars, ...)$summary
}
