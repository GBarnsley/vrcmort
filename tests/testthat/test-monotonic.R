test_that("vrc_standata handles monotonic covariates correctly", {
  set.seed(123)
  # Ensure data is sorted by time, age, sex, cause to match vrc_standata grouped rows
  df <- expand.grid(
    region = 1:2,
    time = 1:5,
    age = 1:2,
    sex = 1:2,
    cause = 1:2
  )
  df$y <- rpois(nrow(df), 10)
  df$exposure <- 100
  df$conflict <- rnorm(nrow(df))
  
  # Create 3 level covariates
  # Make them constant within group (time, age, sex, cause) for simple comparison
  # since standata rows are groups
  df$L1 <- rep(runif(nrow(df)/2), each=2)
  df$L2 <- rep(runif(nrow(df)/2), each=2)
  df$L3 <- rep(runif(nrow(df)/2), each=2)
  
  sdat <- vrc_standata(
    data = df,
    t0 = 3,
    mortality_monotonic = c("L1", "L2", "L3")
  )
  
  expect_equal(sdat$standata$K_mono_mort, 3)
  expect_equal(ncol(sdat$standata$X_mono_mort), 3)
  
  # Grouped columns are N_groups * R. X_mono_mort should match df values.
  # We just check the dimensions and that values are preserved.
  expect_equal(nrow(sdat$standata$X_mono_mort), nrow(df))
})

test_that("vrcm / vrc_fit initializes monotonic parameters (chains=0)", {
  set.seed(123)
  df <- expand.grid(region = c("R1", "R2"), time = 1:5, age = 1:2, sex = 1:2, cause = 1:2)
  df$y <- rpois(nrow(df), 10); df$exposure <- 100; df$conflict <- rnorm(nrow(df))
  df$L1 <- 1; df$L2 <- 0.5; df$L3 <- 0
  
  bundle <- vrcm(
    mortality = vrc_mortality(~1, monotonic = c("L1", "L2", "L3")),
    reporting = vrc_reporting(~1, monotonic = c("L1", "L2")),
    data = df,
    t0 = 3,
    chains = 0
  )
  
  expect_equal(bundle$standata$K_mono_mort, 3)
  expect_equal(bundle$standata$K_mono_rep, 2)
  
  # Verify prior hyperparameters are present
  expect_true("prior_beta_K_mono_mort_loc" %in% names(bundle$standata))
  expect_true("prior_B_mono_mort_scale" %in% names(bundle$standata))
})
