testthat::test_that("build_srr() default returns no linkages", {
  s <- Rceattle::build_srr()
  testthat::expect_type(s, "list")
  testthat::expect_equal(s$srr_fun, 0L)
  testthat::expect_null(s$linkages)
})


testthat::test_that("build_srr() accepts string srr_fun (parity with int)", {
  testthat::expect_equal(Rceattle::build_srr(srr_fun = "mean")$srr_fun, 0L)
  testthat::expect_equal(
    Rceattle::build_srr(srr_fun = "BevertonHolt")$srr_fun, 2L
  )
  testthat::expect_equal(Rceattle::build_srr(srr_fun = "Ricker")$srr_fun, 4L)
  # No string aliases for the soft-deprecated env-driven codes.
  testthat::expect_error(
    Rceattle::build_srr(srr_fun = "mean_env"),
    "unknown `srr_fun`"
  )
})


testthat::test_that("build_srr(srr_fun = 1|3|5) still works but warns", {
  testthat::expect_warning(
    res1 <- Rceattle::build_srr(srr_fun = 1),
    "soft-deprecated"
  )
  testthat::expect_equal(res1$srr_fun, 1L)
  testthat::expect_warning(
    res3 <- Rceattle::build_srr(srr_fun = 3),
    "soft-deprecated"
  )
  testthat::expect_warning(
    res5 <- Rceattle::build_srr(srr_fun = 5),
    "soft-deprecated"
  )
})


testthat::test_that("build_srr(srr_indices = ...) emits deprecation warning", {
  # NA / not supplied -> no warning.
  testthat::expect_silent(Rceattle::build_srr(srr_indices = NA))

  # Any other value triggers the soft-deprecation warning.
  testthat::expect_warning(
    Rceattle::build_srr(srr_indices = 1),
    "deprecated"
  )
  testthat::expect_warning(
    Rceattle::build_srr(srr_indices = c(1, 2, 3)),
    "deprecated"
  )
})


testthat::test_that("build_srr(linkages = ...) attaches and infers param", {
  s <- Rceattle::build_srr(
    srr_fun = 2,                   # Beverton-Holt
    linkages = list(
      alpha = Rceattle::linkage_spec(
        formula = ~ temp,
        by      = ~ species,
        priors  = list(temp = normal(0, 0.5))
      )
    )
  )
  testthat::expect_equal(s$linkages$alpha$param, "alpha")
  fam <- s$linkages$alpha$priors$temp$family
  testthat::expect_equal(fam, "normal")
})


testthat::test_that("build_srr() accepts R0 and beta keys", {
  s <- Rceattle::build_srr(
    srr_fun = 2,
    linkages = list(
      R0   = Rceattle::linkage_spec(formula = ~ temp, by = ~ species),
      beta = Rceattle::linkage_spec(formula = ~ temp, by = ~ species)
    )
  )
  testthat::expect_equal(s$linkages$R0$param, "R0")
  testthat::expect_equal(s$linkages$beta$param, "beta")
})


testthat::test_that("build_srr() rejects unknown linkage parameter names", {
  testthat::expect_error(
    Rceattle::build_srr(
      linkages = list(K = Rceattle::linkage_spec(formula = ~ 1))
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
      linkages = list(R0 = list(formula = ~ temp))   # not a spec
    ),
    "must be a linkage_spec"
  )
})


testthat::test_that("build_srr() warns on alpha/beta linkages with srr_fun = 0", {
  testthat::expect_warning(
    Rceattle::build_srr(
      srr_fun = 0,
      linkages = list(
        alpha = Rceattle::linkage_spec(formula = ~ temp, by = ~ species)
      )
    ),
    "does not use alpha"
  )
  # R0 alone with srr_fun = 0 is fine -- mean recruitment uses R0.
  testthat::expect_silent(
    Rceattle::build_srr(
      srr_fun = 0,
      linkages = list(
        R0 = Rceattle::linkage_spec(formula = ~ temp, by = ~ species)
      )
    )
  )
})


testthat::test_that("build_srr() accepts list of specs per param", {
  s <- Rceattle::build_srr(
    srr_fun = 2,
    linkages = list(
      alpha = list(
        Rceattle::linkage_spec(formula = ~ temp,
                               by      = ~ species,
                               species = 1L),
        Rceattle::linkage_spec(formula = ~ temp + PDO,
                               by      = ~ species,
                               species = 2L)
      )
    )
  )
  testthat::expect_length(s$linkages$alpha, 2L)
  testthat::expect_equal(s$linkages$alpha[[1]]$param, "alpha")
  testthat::expect_equal(s$linkages$alpha[[2]]$param, "alpha")
  testthat::expect_equal(s$linkages$alpha[[1]]$species, 1L)
  testthat::expect_equal(s$linkages$alpha[[2]]$species, 2L)
})
