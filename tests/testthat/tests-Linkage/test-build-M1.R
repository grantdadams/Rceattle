testthat::test_that("build_M1() default returns canonical integer codes", {
  m <- Rceattle::build_M1()
  testthat::expect_equal(m$M1_model, 0L)
  testthat::expect_equal(m$M1_re,    0L)
  testthat::expect_null(m$linkages)
})


testthat::test_that("build_M1() accepts string M1_model (parity with int)", {
  testthat::expect_equal(
    Rceattle::build_M1(M1_model = "fixed")$M1_model, 0L
  )
  testthat::expect_equal(
    Rceattle::build_M1(M1_model = "sex_age_invariant")$M1_model, 1L
  )
  testthat::expect_equal(
    Rceattle::build_M1(M1_model = "sex_specific")$M1_model, 2L
  )
  testthat::expect_equal(
    Rceattle::build_M1(M1_model = "sex_age_specific")$M1_model, 3L
  )
  # No string alias for the soft-deprecated env-driven integer codes
  # (4, 5); their use is discouraged in new code.
  testthat::expect_error(
    Rceattle::build_M1(M1_model = "env_sex_invariant"),
    "unknown `M1_model`"
  )
  testthat::expect_error(
    Rceattle::build_M1(M1_model = "env_sex_specific"),
    "unknown `M1_model`"
  )
})


testthat::test_that("build_M1(M1_model = 4 or 5) still works but warns", {
  testthat::expect_warning(
    res4 <- Rceattle::build_M1(M1_model = 4),
    "soft-deprecated"
  )
  testthat::expect_equal(res4$M1_model, 4L)
  testthat::expect_warning(
    res5 <- Rceattle::build_M1(M1_model = 5),
    "soft-deprecated"
  )
  testthat::expect_equal(res5$M1_model, 5L)

  # Mixed vector: only the deprecated entries trigger the warning.
  testthat::expect_warning(
    Rceattle::build_M1(M1_model = c(1, 2, 4)),
    "soft-deprecated"
  )
})


testthat::test_that("build_M1(M1_indices = ...) emits deprecation warning", {
  # NA / NA_integer_ / NA_real_ are "not supplied" -> no warning.
  testthat::expect_silent(Rceattle::build_M1(M1_indices = NA))
  testthat::expect_silent(Rceattle::build_M1(M1_indices = NA_integer_))

  # Any other value triggers the soft-deprecation warning.
  testthat::expect_warning(
    Rceattle::build_M1(M1_indices = 5),
    "soft-deprecated"
  )
  testthat::expect_warning(
    Rceattle::build_M1(M1_indices = c(2, 3)),
    "soft-deprecated"
  )
})


testthat::test_that("build_M1() accepts string M1_re (parity with int)", {
  testthat::expect_equal(
    Rceattle::build_M1(M1_re = "none")$M1_re, 0L
  )
  testthat::expect_equal(
    Rceattle::build_M1(M1_re = "iid_age")$M1_re, 1L
  )
  testthat::expect_equal(
    Rceattle::build_M1(M1_re = "iid_year")$M1_re, 2L
  )
  testthat::expect_equal(
    Rceattle::build_M1(M1_re = "iid_age_year")$M1_re, 3L
  )
  testthat::expect_equal(
    Rceattle::build_M1(M1_re = "ar1_age")$M1_re, 4L
  )
  testthat::expect_equal(
    Rceattle::build_M1(M1_re = "ar1_year")$M1_re, 5L
  )
  testthat::expect_equal(
    Rceattle::build_M1(M1_re = "ar1_age_year")$M1_re, 6L
  )
})


testthat::test_that("build_M1() accepts vector forms (per-species)", {
  m <- Rceattle::build_M1(
    M1_model = c("sex_age_invariant", "sex_specific"),
    M1_re    = c("none", "iid_age")
  )
  testthat::expect_equal(m$M1_model, c(1L, 2L))
  testthat::expect_equal(m$M1_re,    c(0L, 1L))

  # Integer-vector form is equivalent.
  m2 <- Rceattle::build_M1(M1_model = c(1, 2), M1_re = c(0, 1))
  testthat::expect_equal(m2$M1_model, c(1L, 2L))
  testthat::expect_equal(m2$M1_re,    c(0L, 1L))
})


testthat::test_that("build_M1() rejects unknown string values", {
  testthat::expect_error(
    Rceattle::build_M1(M1_model = "Lorenzen"),
    "unknown `M1_model`"
  )
  testthat::expect_error(
    Rceattle::build_M1(M1_re = "rw_age"),
    "unknown `M1_re`"
  )
  testthat::expect_error(
    Rceattle::build_M1(M1_model = 99),
    "integer `M1_model` must be one of"
  )
  testthat::expect_error(
    Rceattle::build_M1(M1_model = TRUE),
    "must be a string or integer"
  )
})


testthat::test_that("build_M1(linkages = ...) attaches and infers param", {
  m <- Rceattle::build_M1(
    M1_model = "sex_age_invariant",
    linkages = list(
      log_M = Rceattle::linkage_spec(
        formula = ~ temp,
        by      = ~ species,
        priors  = list(temp = normal(0, 0.5))
      )
    )
  )
  testthat::expect_equal(m$linkages$log_M$param, "log_M")
  testthat::expect_equal(m$linkages$log_M$priors$temp$family, "normal")
})


testthat::test_that("build_M1() rejects unknown linkage parameter names", {
  testthat::expect_error(
    Rceattle::build_M1(
      linkages = list(log_K = Rceattle::linkage_spec(formula = ~ 1))
    ),
    "unknown M linkage parameter"
  )
})


testthat::test_that("build_M1() accepts list of specs per param", {
  m <- Rceattle::build_M1(
    linkages = list(
      log_M = list(
        Rceattle::linkage_spec(formula = ~ temp,
                               by      = ~ species,
                               species = 1L),
        Rceattle::linkage_spec(formula = ~ temp + PDO,
                               by      = ~ species,
                               species = 2L)
      )
    )
  )
  testthat::expect_length(m$linkages$log_M, 2L)
  testthat::expect_equal(m$linkages$log_M[[1]]$param, "log_M")
  testthat::expect_equal(m$linkages$log_M[[2]]$param, "log_M")
  testthat::expect_equal(m$linkages$log_M[[1]]$species, 1L)
  testthat::expect_equal(m$linkages$log_M[[2]]$species, 2L)
})
