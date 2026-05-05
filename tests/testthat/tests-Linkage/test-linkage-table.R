testthat::test_that("new_linkage_table() yields an empty schema table", {
  tbl <- Rceattle:::new_linkage_table()
  testthat::expect_s3_class(tbl, "Rceattle_linkage_table")
  testthat::expect_equal(nrow(tbl), 0L)
  testthat::expect_setequal(names(tbl), names(Rceattle:::LINKAGE_COLS))
})


testthat::test_that("validate_linkage_table catches bad cols/enums", {
  good <- Rceattle:::linkage_row(process = "M", param = "log_M", X_col = 1L)
  testthat::expect_silent(Rceattle:::validate_linkage_table(good))

  bad1 <- good
  bad1$process <- NULL
  testthat::expect_error(
    Rceattle:::validate_linkage_table(bad1),
    "missing columns"
  )

  bad2 <- good
  bad2$process <- "predation"
  testthat::expect_error(
    Rceattle:::validate_linkage_table(bad2),
    "unknown process"
  )

  bad3 <- good
  bad3$lower <- 1
  bad3$upper <- 0
  testthat::expect_error(
    Rceattle:::validate_linkage_table(bad3),
    "lower > upper"
  )
})


testthat::test_that("linkage_row() builds a one-row table with defaults", {
  r <- Rceattle:::linkage_row(process = "recruitment",
                              param = "log_alpha",
                              X_col = 2L,
                              species = 1L,
                              link = "log",
                              init = -1)
  testthat::expect_equal(nrow(r), 1L)
  testthat::expect_equal(r$process, "recruitment")
  testthat::expect_equal(r$param, "log_alpha")
  testthat::expect_equal(r$X_col, 2L)
  testthat::expect_true(is.na(r$sex))
  testthat::expect_true(is.na(r$age_bin))
  testthat::expect_equal(r$link, "log")
  testthat::expect_equal(r$init, -1)
  testthat::expect_equal(r$est_phase, 1L)
})


testthat::test_that("bind_linkage() preserves schema and class across rbind", {
  a <- Rceattle:::linkage_row("M", "log_M", X_col = 1L, species = 1L)
  b <- Rceattle:::linkage_row("M", "log_M", X_col = 1L, species = 2L)
  c <- Rceattle:::linkage_row("growth", "log_K", X_col = 1L, species = 1L)
  out <- Rceattle:::bind_linkage(a, b, c)
  testthat::expect_s3_class(out, "Rceattle_linkage_table")
  testthat::expect_equal(nrow(out), 3L)
  testthat::expect_setequal(names(out), names(Rceattle:::LINKAGE_COLS))

  # List form
  out2 <- Rceattle:::bind_linkage(list(a, b, c))
  testthat::expect_equal(out, out2)

  # Empty input -> empty schema
  empty <- Rceattle:::bind_linkage()
  testthat::expect_equal(nrow(empty), 0L)
  testthat::expect_s3_class(empty, "Rceattle_linkage_table")
})


testthat::test_that("linkage_spec() validates formula shape", {
  testthat::expect_error(
    Rceattle:::linkage_spec(formula = "temp", param = "log_M"),
    "must be an R formula"
  )
  testthat::expect_error(
    Rceattle:::linkage_spec(formula = y ~ temp, param = "log_M"),
    "one-sided"
  )
  spec <- Rceattle:::linkage_spec(formula = ~ temp + PDO,
                                  param = "log_alpha",
                                  by = ~species,
                                  link = "log")
  testthat::expect_s3_class(spec, "Rceattle_linkage_spec")
  testthat::expect_equal(spec$param, "log_alpha")
})


testthat::test_that("materialize_linkage row count = ncol(X) x by", {
  env <- data.frame(
    Year = 2000:2009,
    temp = stats::rnorm(10),
    PDO  = stats::rnorm(10)
  )
  spec <- Rceattle:::linkage_spec(
    formula = ~ temp + PDO,
    param   = "log_alpha",
    by      = ~species,
    link    = "log",
    init    = list(`(Intercept)` = -1, temp = 0.1)
  )
  rows <- Rceattle:::materialize_linkage(
    spec, process = "recruitment",
    env_data = env,
    strata   = list(species = 1:3)
  )
  testthat::expect_s3_class(rows, "Rceattle_linkage_table")
  # 3 design columns (intercept + temp + PDO) x 3 species = 9 rows
  testthat::expect_equal(nrow(rows), 9L)
  testthat::expect_setequal(rows$species, 1:3)
  testthat::expect_setequal(rows$X_col, 1:3)
  testthat::expect_true(all(is.na(rows$sex)))
  testthat::expect_true(all(rows$prior_family == "none"))

  # Initial values keyed on design column name propagate
  intercept_rows <- rows[rows$X_col == 1L, ]
  temp_rows      <- rows[rows$X_col == 2L, ]
  pdo_rows       <- rows[rows$X_col == 3L, ]
  testthat::expect_true(all(intercept_rows$init == -1))
  testthat::expect_true(all(temp_rows$init      == 0.1))
  testthat::expect_true(all(pdo_rows$init       == 0))   # default

  # Design colnames stashed for later remapping
  testthat::expect_equal(
    attr(rows, "design_colnames"),
    c("(Intercept)", "temp", "PDO")
  )
})


testthat::test_that("materialize_linkage() handles species + sex grouping", {
  env <- data.frame(Year = 2000:2004, temp = stats::rnorm(5))
  spec <- Rceattle:::linkage_spec(
    formula = ~ temp,        # intercept + temp = 2 cols
    param   = "log_M",
    by      = ~species + sex,
    link    = "log"
  )
  rows <- Rceattle:::materialize_linkage(
    spec, process = "M",
    env_data = env,
    strata   = list(species = 1:2, sex = 1:2)
  )
  # 2 cols x 2 species x 2 sexes = 8 rows
  testthat::expect_equal(nrow(rows), 8L)
  testthat::expect_setequal(rows$species, 1:2)
  testthat::expect_setequal(rows$sex, 1:2)
})


testthat::test_that("materialize_linkage rejects unknown grouping vars", {
  env <- data.frame(Year = 2000:2002, temp = 1:3)
  spec <- Rceattle:::linkage_spec(
    formula = ~ temp,
    param   = "log_M",
    by      = ~stock           # not allowed
  )
  testthat::expect_error(
    Rceattle:::materialize_linkage(
      spec, process = "M",
      env_data = env,
      strata   = list(stock = 1:2)
    ),
    "unknown grouping variable"
  )
})
