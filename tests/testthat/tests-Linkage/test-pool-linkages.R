testthat::test_that("pool_linkages() returns empty pool with no specs", {
  env <- data.frame(Year = 2000:2009, temp = stats::rnorm(10))
  out <- Rceattle:::pool_linkages(spec_groups = NULL, env_data = env)
  testthat::expect_s3_class(out$table, "Rceattle_linkage_table")
  testthat::expect_equal(nrow(out$table), 0L)
  testthat::expect_equal(ncol(out$X), 0L)
  testthat::expect_equal(nrow(out$X), 10L)

  out2 <- Rceattle:::pool_linkages(
    spec_groups = list(growth = list()),
    env_data = env
  )
  testthat::expect_equal(nrow(out2$table), 0L)
})


testthat::test_that("pool_linkages() errors when specs require env_data", {
  spec <- Rceattle::linkage_spec(formula = ~ temp, param = "K")
  testthat::expect_error(
    Rceattle:::pool_linkages(
      spec_groups = list(growth = list(K = spec)),
      env_data    = NULL,
      strata      = list(species = 1L)
    ),
    "env_data"
  )
})


testthat::test_that("pool_linkages() unifies design cols across specs", {
  env <- data.frame(
    Year = 2000:2004,
    temp = c(0.1, 0.2, 0.0, -0.1, 0.3),
    PDO  = c(-1, 0, 1, 0, -1)
  )
  # Two specs share `(Intercept)` and `temp`; second adds `PDO`.
  s_K <- Rceattle::linkage_spec(
    formula = ~ temp,
    param   = "K",
    by      = ~ species
  )
  s_Linf <- Rceattle::linkage_spec(
    formula = ~ temp + PDO,
    param   = "Linf",
    by      = ~ species
  )
  out <- Rceattle:::pool_linkages(
    spec_groups = list(growth = list(K = s_K, Linf = s_Linf)),
    env_data    = env,
    strata      = list(species = 1:2)
  )

  testthat::expect_s3_class(out$table, "Rceattle_linkage_table")
  # K: 2 cols (intercept + temp) x 2 species = 4
  # Linf: 3 cols (intercept + temp + PDO) x 2 species = 6
  testthat::expect_equal(nrow(out$table), 10L)
  testthat::expect_setequal(out$design_names, c("(Intercept)", "temp", "PDO"))
  testthat::expect_equal(ncol(out$X), 3L)
  testthat::expect_equal(nrow(out$X), nrow(env))
  testthat::expect_equal(unname(out$X[, "temp"]), env$temp)
  testthat::expect_equal(unname(out$X[, "PDO"]),  env$PDO)
  testthat::expect_equal(unname(out$X[, "(Intercept)"]),
                         rep(1, nrow(env)))

  # X_col must reference valid global columns
  testthat::expect_true(all(out$table$X_col %in% seq_along(out$design_names)))
  # Each row's X_col, when looked up in design_names, matches the
  # column name expected for that param. K only ever uses cols
  # (Intercept, temp); Linf can additionally use PDO.
  k_cols    <- out$design_names[out$table$X_col[out$table$param == "K"]]
  linf_cols <- out$design_names[out$table$X_col[out$table$param == "Linf"]]
  testthat::expect_setequal(k_cols, c("(Intercept)", "temp"))
  testthat::expect_setequal(linf_cols, c("(Intercept)", "temp", "PDO"))
})


testthat::test_that("pool_linkages() preserves prior info on rows", {
  env <- data.frame(Year = 2000:2004, temp = stats::rnorm(5))
  s_K <- Rceattle::linkage_spec(
    formula = ~ temp,
    param   = "K",
    by      = ~ species,
    priors  = list(temp = normal(0, 0.5))
  )
  out <- Rceattle:::pool_linkages(
    spec_groups = list(growth = list(K = s_K)),
    env_data    = env,
    strata      = list(species = 1:3)
  )
  temp_global <- match("temp", out$design_names)
  temp_rows <- out$table[out$table$X_col == temp_global, ]
  testthat::expect_equal(nrow(temp_rows), 3L)        # 1 col x 3 species
  testthat::expect_true(all(temp_rows$prior_family == "normal"))
  testthat::expect_true(all(temp_rows$prior_p1 == 0))
  testthat::expect_true(all(temp_rows$prior_p2 == 0.5))
})


testthat::test_that("pool_linkages drops the per-spec design attrs on output", {
  env <- data.frame(Year = 2000:2004, temp = stats::rnorm(5))
  s <- Rceattle::linkage_spec(formula = ~ temp, param = "K", by = ~ species)
  out <- Rceattle:::pool_linkages(
    spec_groups = list(growth = list(K = s)),
    env_data    = env,
    strata      = list(species = 1:2)
  )
  testthat::expect_null(attr(out$table, "design_colnames"))
  testthat::expect_null(attr(out$table, "design_matrix"))
})
