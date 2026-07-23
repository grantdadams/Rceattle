# data_requirements() — the exported introspection reader over the requirement
# table. Pure classification (no fit), so unguarded/fast.

testthat::test_that("single-species: predation Ignored, comp/caal Optional", {
  r <- data_requirements(msmMode = 0)
  st <- setNames(r$status, r$element)
  testthat::expect_identical(unname(st["diet_data"]), "Ignored")
  testthat::expect_identical(unname(st["ration_data"]), "Ignored")
  testthat::expect_identical(unname(st["bioenergetics"]), "Ignored")
  testthat::expect_identical(unname(st["comp_data"]), "Optional")
  testthat::expect_identical(unname(st["caal_data"]), "Optional")
  # Core backbone always Required.
  for (nm in c("nspp", "styr", "endyr", "fleet_control", "weight",
               "maturity", "sex_ratio", "M1_base"))
    testthat::expect_identical(unname(st[nm]), "Required", info = nm)
})

testthat::test_that("multispecies flips the predation trio to Required", {
  r0 <- data_requirements(msmMode = 0)
  r1 <- data_requirements(msmMode = 1)
  s0 <- setNames(r0$status, r0$element); s1 <- setNames(r1$status, r1$element)
  for (nm in c("diet_data", "ration_data", "bioenergetics")) {
    testthat::expect_identical(unname(s0[nm]), "Ignored", info = nm)
    testthat::expect_identical(unname(s1[nm]), "Required", info = nm)
  }
})

testthat::test_that("conditional requirements respond to the switches", {
  # growth_model > 0 -> caal_data Required
  r <- data_requirements(growth_model = 1)
  testthat::expect_identical(
    r$status[r$element == "caal_data"], "Required")

  # estDynamics > 0 -> NByageFixed Required
  r <- data_requirements(estDynamics = 1)
  testthat::expect_identical(
    r$status[r$element == "NByageFixed"], "Required")

  # Fixed selectivity -> emp_sel Required; else Ignored
  r <- data_requirements(selectivity = c("Logistic", "Fixed"))
  testthat::expect_identical(r$status[r$element == "emp_sel"], "Required")
  r <- data_requirements(selectivity = c("Logistic", "Logistic"))
  testthat::expect_identical(r$status[r$element == "emp_sel"], "Ignored")

  # MVN index -> index_cov Required; else Ignored
  r <- data_requirements(index_loglike = c("Lognormal", "MVN"))
  testthat::expect_identical(r$status[r$element == "index_cov"], "Required")
  r <- data_requirements(index_loglike = c("Lognormal", "Lognormal"))
  testthat::expect_identical(r$status[r$element == "index_cov"], "Ignored")

  # Ceq > 1 -> env_data Required; else Optional
  r <- data_requirements(Ceq = 2)
  testthat::expect_identical(r$status[r$element == "env_data"], "Required")
  r <- data_requirements(Ceq = 1)
  testthat::expect_identical(r$status[r$element == "env_data"], "Optional")
})

testthat::test_that("return shape is a tidy ordered data.frame", {
  r <- data_requirements(msmMode = 0)
  testthat::expect_s3_class(r, "data.frame")
  testthat::expect_setequal(
    colnames(r), c("element", "category", "status", "condition", "default"))
  testthat::expect_true(all(r$status %in% c("Required", "Optional", "Ignored")))
  # Required rows come first.
  rank <- match(r$status, c("Required", "Optional", "Ignored"))
  testthat::expect_false(is.unsorted(rank))
  # Optional rows carry a default description; others do not.
  testthat::expect_true(all(nzchar(r$default[r$status == "Optional"])))
  testthat::expect_true(all(!nzchar(r$default[r$status != "Optional"])))
})

testthat::test_that("accepts a real (bundled) data list", {
  testthat::skip_if_not_installed("Rceattle")
  # A bundled data list may STORE a different msmMode than the canonical fit
  # uses (BS2017SS stores msmMode = 1). By default data_requirements() reflects
  # the stored config...
  r <- suppressWarnings(suppressMessages(data_requirements(BS2017SS)))
  st <- setNames(r$status, r$element)
  testthat::expect_identical(unname(st["fleet_control"]), "Required")
  testthat::expect_identical(unname(st["diet_data"]),
                             if (isTRUE(BS2017SS$msmMode > 0)) "Required" else "Ignored")

  # ...and an explicit switch arg overrides it (fit_mod precedence), so a
  # single-species preview Ignores the predation inputs.
  r0 <- suppressWarnings(suppressMessages(data_requirements(BS2017SS, msmMode = 0)))
  st0 <- setNames(r0$status, r0$element)
  testthat::expect_identical(unname(st0["diet_data"]), "Ignored")
  testthat::expect_identical(unname(st0["bioenergetics"]), "Ignored")
})
