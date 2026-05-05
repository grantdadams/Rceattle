testthat::test_that("prior_* constructors build valid Rceattle_prior", {
  p1 <- Rceattle::prior_normal(0, 1)
  testthat::expect_s3_class(p1, "Rceattle_prior")
  testthat::expect_equal(p1$family, "normal")
  testthat::expect_equal(p1$p1, 0)
  testthat::expect_equal(p1$p2, 1)

  p2 <- Rceattle::prior_lognormal(-1, 0.5)
  testthat::expect_equal(p2$family, "lognormal")

  p3 <- Rceattle::prior_gamma(2, 3)
  testthat::expect_equal(p3$family, "gamma")
  testthat::expect_equal(c(p3$p1, p3$p2), c(2, 3))

  p4 <- Rceattle::prior_beta(2, 5)
  testthat::expect_equal(p4$family, "beta")
})


testthat::test_that("prior constructors reject bad parameters", {
  testthat::expect_error(Rceattle::prior_normal(0, 0),    "sd")
  testthat::expect_error(Rceattle::prior_normal(0, -1),   "sd")
  testthat::expect_error(Rceattle::prior_lognormal(0, 0), "sdlog")
  testthat::expect_error(Rceattle::prior_gamma(0, 1),     "shape")
  testthat::expect_error(Rceattle::prior_gamma(1, 0),     "shape")
  testthat::expect_error(Rceattle::prior_beta(0, 1),      "shape")
  testthat::expect_error(Rceattle::prior_normal("a", 1),  "finite numeric")
  testthat::expect_error(Rceattle::prior_normal(c(0, 0), 1), "finite numeric")
})


testthat::test_that("linkage_spec resolves shorthand prior names via NSE", {
  spec <- Rceattle:::linkage_spec(
    formula = ~ temp + PDO,
    param   = "log_alpha",
    by      = ~species,
    link    = "log",
    priors  = list(
      temp = normal(0, 1),
      PDO  = lognormal(0, 0.5)
    )
  )
  testthat::expect_s3_class(spec, "Rceattle_linkage_spec")
  testthat::expect_s3_class(spec$priors$temp, "Rceattle_prior")
  testthat::expect_equal(spec$priors$temp$family, "normal")
  testthat::expect_equal(spec$priors$PDO$family, "lognormal")
})


testthat::test_that("NSE in priors does not mask base::gamma/base::beta", {
  # Outside the priors arg, base::gamma must still be the math function.
  testthat::expect_equal(gamma(5), 24)
  testthat::expect_equal(beta(2, 3), 1 / 12)

  # Inside priors, gamma() and beta() resolve to prior constructors.
  spec <- Rceattle:::linkage_spec(
    formula = ~ temp,
    param   = "log_M",
    priors  = list(`(Intercept)` = gamma(2, 3), temp = beta(2, 5))
  )
  testthat::expect_equal(spec$priors$`(Intercept)`$family, "gamma")
  testthat::expect_equal(spec$priors$temp$family, "beta")
})


testthat::test_that("linkage_spec accepts programmatic prior_* construction", {
  pl <- list(temp = Rceattle::prior_normal(0, 1))
  spec <- Rceattle:::linkage_spec(
    formula = ~ temp,
    param   = "log_M",
    priors  = pl
  )
  testthat::expect_equal(spec$priors$temp$family, "normal")
})


testthat::test_that("linkage_spec rejects non-prior values in priors", {
  testthat::expect_error(
    Rceattle:::linkage_spec(
      formula = ~ temp,
      param   = "log_M",
      priors  = list(temp = c(0, 1))
    ),
    "must be an Rceattle_prior"
  )
  testthat::expect_error(
    Rceattle:::linkage_spec(
      formula = ~ temp,
      param   = "log_M",
      priors  = list(normal(0, 1))   # unnamed at top level
    ),
    "named list"
  )
})


testthat::test_that("species-specific priors materialize per row", {
  env <- data.frame(Year = 2000:2004, temp = stats::rnorm(5))
  spec <- Rceattle:::linkage_spec(
    formula = ~ temp,
    param   = "log_alpha",
    by      = ~ species,
    priors  = list(
      temp = list(
        `1` = normal(0, 0.1),
        `2` = normal(0, 0.5)
        # species 3 deliberately absent -> no prior
      )
    )
  )
  rows <- Rceattle:::materialize_linkage(
    spec, process = "recruitment",
    env_data = env,
    strata   = list(species = 1:3)
  )
  temp_rows <- rows[rows$X_col == 2L, ]
  testthat::expect_equal(nrow(temp_rows), 3L)

  sp1 <- temp_rows[temp_rows$species == 1L, ]
  sp2 <- temp_rows[temp_rows$species == 2L, ]
  sp3 <- temp_rows[temp_rows$species == 3L, ]
  testthat::expect_equal(sp1$prior_family, "normal")
  testthat::expect_equal(sp1$prior_p2, 0.1)
  testthat::expect_equal(sp2$prior_family, "normal")
  testthat::expect_equal(sp2$prior_p2, 0.5)
  testthat::expect_equal(sp3$prior_family, "none")
  testthat::expect_true(is.na(sp3$prior_p1))
})


testthat::test_that("species-specific priors validate input shape", {
  testthat::expect_error(
    Rceattle:::linkage_spec(
      formula = ~ temp,
      param   = "log_M",
      priors  = list(temp = list(normal(0, 1), normal(0, 0.5)))   # unnamed
    ),
    "named list keyed by species id"
  )
  testthat::expect_error(
    Rceattle:::linkage_spec(
      formula = ~ temp,
      param   = "log_M",
      priors  = list(temp = list(`1` = c(0, 1)))   # not an Rceattle_prior
    ),
    "is not an Rceattle_prior"
  )
})


testthat::test_that("scalar prior still applies across all species", {
  env <- data.frame(Year = 2000:2004, temp = stats::rnorm(5))
  spec <- Rceattle:::linkage_spec(
    formula = ~ temp,
    param   = "log_M",
    by      = ~ species,
    priors  = list(temp = normal(0, 0.3))
  )
  rows <- Rceattle:::materialize_linkage(
    spec, process = "M",
    env_data = env,
    strata   = list(species = 1:3)
  )
  temp_rows <- rows[rows$X_col == 2L, ]
  testthat::expect_true(all(temp_rows$prior_family == "normal"))
  testthat::expect_true(all(temp_rows$prior_p2 == 0.3))
})


testthat::test_that("materialize_linkage propagates prior cols to rows", {
  env <- data.frame(Year = 2000:2009, temp = stats::rnorm(10))
  spec <- Rceattle:::linkage_spec(
    formula = ~ temp,                  # 2 design cols (intercept + temp)
    param   = "log_alpha",
    by      = ~species,
    link    = "log",
    priors  = list(temp = normal(0, 0.3))
  )
  rows <- Rceattle:::materialize_linkage(
    spec, process = "recruitment",
    env_data = env,
    strata   = list(species = 1:2)
  )
  testthat::expect_equal(nrow(rows), 4L)        # 2 cols x 2 species

  intercept_rows <- rows[rows$X_col == 1L, ]
  temp_rows      <- rows[rows$X_col == 2L, ]
  testthat::expect_true(all(intercept_rows$prior_family == "none"))
  testthat::expect_true(all(is.na(intercept_rows$prior_p1)))
  testthat::expect_true(all(temp_rows$prior_family == "normal"))
  testthat::expect_true(all(temp_rows$prior_p1 == 0))
  testthat::expect_true(all(temp_rows$prior_p2 == 0.3))
})


testthat::test_that("validate_linkage_table catches inconsistent prior rows", {
  good <- Rceattle:::linkage_row(
    process = "M", param = "log_M", X_col = 1L,
    prior_family = "normal", prior_p1 = 0, prior_p2 = 1
  )
  testthat::expect_silent(Rceattle:::validate_linkage_table(good))

  # prior_family != "none" but params are NA
  # Construct directly to bypass linkage_row()'s own validation.
  bad <- good
  bad$prior_p1 <- NA_real_
  bad$prior_p2 <- NA_real_
  testthat::expect_error(
    Rceattle:::validate_linkage_table(bad),
    "non-NA"
  )

  # unknown family
  bad2 <- good
  bad2$prior_family <- "cauchy"
  testthat::expect_error(
    Rceattle:::validate_linkage_table(bad2),
    "unknown prior family"
  )
})
