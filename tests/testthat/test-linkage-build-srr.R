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
  # No string aliases for the retired env-driven codes.
  testthat::expect_error(
    Rceattle::build_srr(srr_fun = "mean_env"),
    "unknown `srr_fun`"
  )
})


# From 4.4.0 the template no longer applied the environmental recruitment term,
# but srr_fun 1/3/5 and srr_indices kept running with only a soft-deprecation
# warning, so a model was fitted without its covariate and reported no error.
# They now stop and point to the linkage that replaces them.
testthat::test_that("build_srr(srr_fun / srr_pred_fun = 1|3|5) is an error", {
  for (code in c(1, 3, 5)) {
    testthat::expect_error(Rceattle::build_srr(srr_fun = code), "linkages",
                           info = code)
    testthat::expect_error(Rceattle::build_srr(srr_pred_fun = code), "linkages",
                           info = code)
  }
})


testthat::test_that("build_srr(srr_indices = ...) is an error; NA or NULL is not supplied", {
  # .refit_like() passes data_list$srr_indices, which a current fit leaves NULL.
  testthat::expect_silent(Rceattle::build_srr(srr_indices = NA))
  testthat::expect_silent(Rceattle::build_srr(srr_indices = NULL))
  testthat::expect_error(Rceattle::build_srr(srr_indices = 1), "linkages")
  testthat::expect_error(Rceattle::build_srr(srr_indices = c(1, 2, 3)), "linkages")
})


testthat::test_that("a refit maps a stored code 1/3/5 to the form it fitted, with a warning", {
  # A fit made before 5.32.0 can carry code 1, 3 or 5; from 4.4.0 those fitted
  # 0, 2 or 4, so the refitting diagnostics must still rebuild it.
  for (code in c(1L, 3L, 5L)) {
    testthat::expect_warning(
      mapped <- Rceattle:::.srr_fun_structural(code), "no effect since 4.4.0"
    )
    testthat::expect_identical(mapped, code - 1L)
  }
  for (code in c(0L, 2L, 4L)) {
    testthat::expect_silent(mapped <- Rceattle:::.srr_fun_structural(code))
    testthat::expect_identical(mapped, code)
  }
  testthat::expect_null(Rceattle:::.srr_fun_structural(NULL))
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
