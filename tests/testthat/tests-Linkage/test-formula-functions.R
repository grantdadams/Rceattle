testthat::test_that("I() squared term materializes a separate design col", {
  env <- data.frame(Year = 2000:2009, temp = stats::rnorm(10))
  spec <- Rceattle::linkage_spec(
    formula = ~ temp + I(temp^2),
    param   = "K",
    by      = ~ species
  )
  rows <- Rceattle:::materialize_linkage(
    spec, process = "growth",
    env_data = env, strata = list(species = 1:2)
  )
  testthat::expect_setequal(attr(rows, "design_colnames"),
                            c("(Intercept)", "temp", "I(temp^2)"))
  # 3 design cols x 2 species = 6 rows
  testthat::expect_equal(nrow(rows), 6L)
})


testthat::test_that("poly(temp, 4) materializes the 4 orthogonal columns", {
  env <- data.frame(Year = 2000:2019, temp = stats::rnorm(20))
  spec <- Rceattle::linkage_spec(
    formula = ~ poly(temp, 4),
    param   = "K",
    by      = ~ species
  )
  rows <- Rceattle:::materialize_linkage(
    spec, process = "growth",
    env_data = env, strata = list(species = 1L)
  )
  cols <- attr(rows, "design_colnames")
  testthat::expect_equal(length(cols), 5L)
  testthat::expect_true("(Intercept)" %in% cols)
  testthat::expect_true(all(c("poly(temp, 4)1", "poly(temp, 4)2",
                              "poly(temp, 4)3", "poly(temp, 4)4") %in% cols))

  x_mat <- attr(rows, "design_matrix")
  testthat::expect_equal(nrow(x_mat), nrow(env))
  # Orthogonal polynomial columns are zero-mean over the env_data rows.
  poly_cols <- grep("^poly", colnames(x_mat), value = TRUE)
  testthat::expect_true(
    all(abs(colMeans(x_mat[, poly_cols, drop = FALSE])) < 1e-10)
  )
})


testthat::test_that("formula functions feed cleanly through pool_linkages", {
  env <- data.frame(Year = 2000:2019, temp = stats::rnorm(20))
  s1 <- Rceattle::linkage_spec(
    formula = ~ temp + I(temp^2),
    by      = ~ species,
    species = 1L
  )
  s2 <- Rceattle::linkage_spec(
    formula = ~ poly(temp, 3),
    by      = ~ species,
    species = 2L
  )
  out <- Rceattle:::pool_linkages(
    spec_groups = list(growth = list(K = list(s1, s2))),
    env_data    = env,
    strata      = list(species = 1:2)
  )

  # All design column names from both specs are unioned in the global X.
  testthat::expect_true("I(temp^2)" %in% out$design_names)
  testthat::expect_true("poly(temp, 3)1" %in% out$design_names)
  # X has the right shape for the pooled design.
  testthat::expect_equal(nrow(out$X), nrow(env))
  testthat::expect_equal(ncol(out$X), length(out$design_names))
})


testthat::test_that("priors and inits key off formula-function colnames", {
  env <- data.frame(Year = 2000:2019, temp = stats::rnorm(20))
  spec <- Rceattle::linkage_spec(
    formula = ~ temp + I(temp^2),
    param   = "K",
    by      = ~ species,
    init    = list(`I(temp^2)` = -0.05),
    priors  = list(`I(temp^2)` = normal(0, 0.1))
  )
  rows <- Rceattle:::materialize_linkage(
    spec, process = "growth",
    env_data = env, strata = list(species = 1:2)
  )
  sq_idx <- match("I(temp^2)", attr(rows, "design_colnames"))
  sq_rows <- rows[rows$X_col == sq_idx, ]
  testthat::expect_true(all(sq_rows$init == -0.05))
  testthat::expect_true(all(sq_rows$prior_family == "normal"))
  testthat::expect_true(all(sq_rows$prior_p2 == 0.1))
})
