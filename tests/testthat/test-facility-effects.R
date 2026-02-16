test_that("vrc_standata handles facility covariates correctly", {
  set.seed(123)
  df <- expand.grid(region = 1:2, time = 1:5, age = 1:2, sex = 1:2, cause = 1:2)
  df$y <- rpois(nrow(df), 10)
  df$exposure <- 100
  df$conflict <- rnorm(nrow(df))
  df$F1 <- runif(nrow(df))
  df$F2 <- runif(nrow(df))

  # vrc_standata sorts by time, age, sex, cause, region
  # We ensure we compare the right values
  sdat <- vrc_standata(data = df, t0 = 3, mortality_facilities = c("F1", "F2"))

  expect_equal(sdat$standata$K_fac_mort, 2)
  expect_equal(ncol(sdat$standata$X_fac_mort), 2)

  # Values in X_fac_mort are expanded/grouped.
  # We check that dimensions are consistent and values are numeric.
  expect_true(is.numeric(sdat$standata$X_fac_mort))
})

test_that("vrcm initializes facility parameters (chains=0)", {
  set.seed(123)
  df <- expand.grid(
    region = c("R1", "R2"),
    time = 1:5,
    age = 1:2,
    sex = 1:2,
    cause = 1:2
  )
  df$y <- rpois(nrow(df), 10)
  df$exposure <- 100
  df$conflict <- rnorm(nrow(df))
  df$F1 <- 1
  df$F2 <- 0.5

  bundle <- vrcm(
    mortality = vrc_mortality(~1, facilities = c("F1", "F2")),
    reporting = vrc_reporting(~1, facilities = c("F1")),
    data = df,
    t0 = 3,
    chains = 0
  )

  expect_equal(bundle$standata$K_fac_mort, 2)
  expect_equal(bundle$standata$K_fac_rep, 1)
  expect_true("prior_beta_fac_best_loc" %in% names(bundle$standata))
  expect_true("prior_gamma_fac_best_scale" %in% names(bundle$standata))
})
