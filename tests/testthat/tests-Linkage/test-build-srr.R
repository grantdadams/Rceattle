testthat::test_that("build_srr() default returns no linkages", {
  s <- Rceattle::build_srr()
  testthat::expect_type(s, "list")
  testthat::expect_equal(s$srr_fun, 0)
  testthat::expect_null(s$linkages)
})


testthat::test_that("build_srr(linkages = ...) attaches and infers param", {
  s <- Rceattle::build_srr(
    srr_fun = 2,                   # Beverton-Holt
    linkages = list(
      log_alpha = Rceattle::linkage_spec(
        formula = ~ temp,
        by      = ~ species,
        priors  = list(temp = normal(0, 0.5))
      )
    )
  )
  testthat::expect_equal(s$linkages$log_alpha$param, "log_alpha")
  testthat::expect_equal(s$linkages$log_alpha$priors$temp$family, "normal")
})


testthat::test_that("build_srr() accepts log_R0 and log_beta keys", {
  s <- Rceattle::build_srr(
    srr_fun = 2,
    linkages = list(
      log_R0   = Rceattle::linkage_spec(formula = ~ temp, by = ~ species),
      log_beta = Rceattle::linkage_spec(formula = ~ temp, by = ~ species)
    )
  )
  testthat::expect_equal(s$linkages$log_R0$param, "log_R0")
  testthat::expect_equal(s$linkages$log_beta$param, "log_beta")
})


testthat::test_that("build_srr() rejects unknown linkage parameter names", {
  testthat::expect_error(
    Rceattle::build_srr(
      linkages = list(log_K = Rceattle::linkage_spec(formula = ~ 1))
    ),
    "unknown recruitment linkage parameter"
  )
})


testthat::test_that("build_srr() requires a named list", {
  testthat::expect_error(
    Rceattle::build_srr(
      linkages = list(Rceattle::linkage_spec(formula = ~ 1))
    ),
    "named list"
  )
})


testthat::test_that("build_srr() rejects non-spec list entries", {
  testthat::expect_error(
    Rceattle::build_srr(
      linkages = list(log_R0 = list(formula = ~ temp))   # not a spec
    ),
    "must be a linkage_spec"
  )
})


testthat::test_that("build_srr() warns on alpha/beta linkages with srr_fun = 0", {
  testthat::expect_warning(
    Rceattle::build_srr(
      srr_fun = 0,
      linkages = list(
        log_alpha = Rceattle::linkage_spec(formula = ~ temp, by = ~ species)
      )
    ),
    "does not use alpha"
  )
  # log_R0 alone with srr_fun = 0 is fine -- mean recruitment uses R0.
  testthat::expect_silent(
    Rceattle::build_srr(
      srr_fun = 0,
      linkages = list(
        log_R0 = Rceattle::linkage_spec(formula = ~ temp, by = ~ species)
      )
    )
  )
})


testthat::test_that("build_srr() accepts list of specs per param", {
  s <- Rceattle::build_srr(
    srr_fun = 2,
    linkages = list(
      log_alpha = list(
        Rceattle::linkage_spec(formula = ~ temp,
                               by      = ~ species,
                               species = 1L),
        Rceattle::linkage_spec(formula = ~ temp + PDO,
                               by      = ~ species,
                               species = 2L)
      )
    )
  )
  testthat::expect_length(s$linkages$log_alpha, 2L)
  testthat::expect_equal(s$linkages$log_alpha[[1]]$param, "log_alpha")
  testthat::expect_equal(s$linkages$log_alpha[[2]]$param, "log_alpha")
  testthat::expect_equal(s$linkages$log_alpha[[1]]$species, 1L)
  testthat::expect_equal(s$linkages$log_alpha[[2]]$species, 2L)
})
