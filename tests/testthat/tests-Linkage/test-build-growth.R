testthat::test_that("build_growth() defaults to empirical", {
  g <- Rceattle::build_growth()
  testthat::expect_type(g, "list")
  testthat::expect_equal(g$fun, "empirical")
  testthat::expect_equal(g$growth_model, 0L)
  testthat::expect_null(g$linkages)

  g2 <- Rceattle::build_growth(fun = "vonBertalanffy")
  testthat::expect_equal(g2$fun, "vonBertalanffy")
  testthat::expect_equal(g2$growth_model, 1L)
  testthat::expect_null(g2$linkages)

  g3 <- Rceattle::build_growth(fun = "Richards")
  testthat::expect_equal(g3$fun, "Richards")
  testthat::expect_equal(g3$growth_model, 2L)
})


testthat::test_that("build_growth() rejects unknown growth functions", {
  testthat::expect_error(
    Rceattle::build_growth(fun = "Schnute"),
    "unknown growth fun"
  )
  testthat::expect_error(
    Rceattle::build_growth(fun = 7),
    "integer `fun` must be one of"
  )
  testthat::expect_error(
    Rceattle::build_growth(fun = TRUE),
    "must be a string or integer"
  )
})


testthat::test_that("build_growth() accepts integer fun (parity with strings)", {
  testthat::expect_equal(Rceattle::build_growth(fun = 0)$fun, "empirical")
  testthat::expect_equal(Rceattle::build_growth(fun = 1)$fun, "vonBertalanffy")
  testthat::expect_equal(Rceattle::build_growth(fun = 2)$fun, "Richards")

  # Integer code on the result matches the canonical mapping
  testthat::expect_equal(Rceattle::build_growth(fun = 1)$growth_model, 1L)
  testthat::expect_equal(Rceattle::build_growth(fun = "Richards")$growth_model,
                         2L)
})


testthat::test_that("build_growth() accepts vector fun for per-species growth", {
  g <- Rceattle::build_growth(fun = c("vonBertalanffy", "Richards"))
  testthat::expect_equal(g$fun, c("vonBertalanffy", "Richards"))
  testthat::expect_equal(g$growth_model, c(1L, 2L))

  # Integer-vector form is equivalent
  g2 <- Rceattle::build_growth(fun = c(1, 2))
  testthat::expect_equal(g2$fun, c("vonBertalanffy", "Richards"))
  testthat::expect_equal(g2$growth_model, c(1L, 2L))

  # Mixed bad entries still error
  testthat::expect_error(
    Rceattle::build_growth(fun = c("vonBertalanffy", "Schnute")),
    "Schnute"
  )
  testthat::expect_error(
    Rceattle::build_growth(fun = c(1, 99)),
    "integer `fun` must be one of"
  )
})


testthat::test_that("build_growth(linkages = ...) attaches and infers param", {
  g <- Rceattle::build_growth(
    fun = "vonBertalanffy",
    linkages = list(
      log_K    = Rceattle::linkage_spec(
        formula = ~ temp,
        by      = ~ species + sex,
        link    = "identity",
        priors  = list(temp = normal(0, 1))
      ),
      log_Linf = Rceattle::linkage_spec(
        formula = ~ 1,
        by      = ~ species + sex
      )
    )
  )
  testthat::expect_named(g$linkages, c("log_K", "log_Linf"))
  testthat::expect_equal(g$linkages$log_K$param, "log_K")
  testthat::expect_equal(g$linkages$log_Linf$param, "log_Linf")
  testthat::expect_equal(g$linkages$log_K$priors$temp$family, "normal")
})


testthat::test_that("build_growth() rejects unknown parameter names", {
  testthat::expect_error(
    Rceattle::build_growth(
      fun = "vonBertalanffy",
      linkages = list(
        log_K    = Rceattle::linkage_spec(formula = ~ 1),
        log_BOGUS = Rceattle::linkage_spec(formula = ~ 1)
      )
    ),
    "unknown growth linkage parameter"
  )
})


testthat::test_that("build_growth() requires a named list", {
  testthat::expect_error(
    Rceattle::build_growth(
      fun = "vonBertalanffy",
      linkages = list(Rceattle::linkage_spec(formula = ~ 1))
    ),
    "named list"
  )
  testthat::expect_error(
    Rceattle::build_growth(
      fun = "vonBertalanffy",
      linkages = "log_K"
    ),
    "named list"
  )
})


testthat::test_that("build_growth() rejects non-spec list entries", {
  testthat::expect_error(
    Rceattle::build_growth(
      fun = "vonBertalanffy",
      linkages = list(log_K = list(formula = ~ temp))
    ),
    "linkage_spec"
  )
})


testthat::test_that("build_growth() rejects log_m under von Bertalanffy", {
  testthat::expect_error(
    Rceattle::build_growth(
      fun = "vonBertalanffy",
      linkages = list(
        log_m = Rceattle::linkage_spec(formula = ~ temp)
      )
    ),
    "Richards"
  )
})


testthat::test_that("build_growth() warns when linkages used with empirical", {
  testthat::expect_warning(
    Rceattle::build_growth(
      fun = "empirical",
      linkages = list(log_K = Rceattle::linkage_spec(formula = ~ 1))
    ),
    "empirical"
  )
})


testthat::test_that("linkage_spec(param=NULL) routes through param setter", {
  s <- Rceattle::linkage_spec(formula = ~ temp)
  testthat::expect_true(is.na(s$param))

  s2 <- Rceattle:::.set_linkage_param(s, "log_K")
  testthat::expect_equal(s2$param, "log_K")

  # Re-setting to a conflicting param errors
  s3 <- Rceattle::linkage_spec(formula = ~ temp, param = "log_K")
  testthat::expect_error(
    Rceattle:::.set_linkage_param(s3, "log_Linf"),
    "conflicts"
  )
})


testthat::test_that("materialize_linkage requires param to be set", {
  bare <- Rceattle::linkage_spec(formula = ~ temp)  # no param
  env  <- data.frame(Year = 2000:2002, temp = 1:3)
  testthat::expect_error(
    Rceattle:::materialize_linkage(
      bare, process = "growth",
      env_data = env, strata = list()
    ),
    "missing a `param`"
  )
})


# --- SD endpoint linkages (log_sd_L1, log_sd_Linf) -------------------------
# These plug into `growth_log_sd[sp, sex, k]` rather than `log_growth_pars`
# and only honor intercept-bearing formulas (no year-dim offset path).

testthat::test_that("build_growth() accepts log_sd_L1 / log_sd_Linf as intercept-only", {
  g <- Rceattle::build_growth(
    fun = "vonBertalanffy",
    linkages = list(
      log_sd_L1   = Rceattle::linkage_spec(
        formula = ~ 1,
        init    = list(`(Intercept)` = log(5)),
        priors  = list(`(Intercept)` = normal(log(5), 0.3))
      ),
      log_sd_Linf = Rceattle::linkage_spec(formula = ~ 1)
    )
  )
  testthat::expect_named(g$linkages, c("log_sd_L1", "log_sd_Linf"))
  testthat::expect_equal(g$linkages$log_sd_L1$param,   "log_sd_L1")
  testthat::expect_equal(g$linkages$log_sd_Linf$param, "log_sd_Linf")
  testthat::expect_equal(g$linkages$log_sd_L1$priors$`(Intercept)`$family,
                         "normal")
})


testthat::test_that("build_growth() errors on slope-only formulas for SD endpoints", {
  testthat::expect_error(
    Rceattle::build_growth(
      fun = "vonBertalanffy",
      linkages = list(
        log_sd_L1 = Rceattle::linkage_spec(formula = ~ 0 + temp)
      )
    ),
    "intercept-bearing"
  )
})


testthat::test_that("build_growth() warns on intercept + slope formulas for SD endpoints", {
  testthat::expect_warning(
    Rceattle::build_growth(
      fun = "vonBertalanffy",
      linkages = list(
        log_sd_Linf = Rceattle::linkage_spec(formula = ~ 1 + temp)
      )
    ),
    "slope terms"
  )
})


testthat::test_that("LINKAGE_PARAM_CODES$growth pins SD endpoints to 4 / 5", {
  # These codes are consumed by linkage.hpp; the env-offset loop guards
  # `param >= RCEATTLE_N_GROWTH_PARAMS (=4)` to drop SD rows, and the
  # prior loop re-targets >=4 to `growth_log_sd(sp, sex, param - 4)`.
  testthat::expect_equal(Rceattle:::LINKAGE_PARAM_CODES$growth[["log_sd_L1"]],   4L)
  testthat::expect_equal(Rceattle:::LINKAGE_PARAM_CODES$growth[["log_sd_Linf"]], 5L)
})


testthat::test_that(".GROWTH_SD_PARAM_TO_INDEX maps names to growth_log_sd slices", {
  testthat::expect_equal(Rceattle:::.GROWTH_SD_PARAM_TO_INDEX,
                         c(log_sd_L1 = 1L, log_sd_Linf = 2L))
})
