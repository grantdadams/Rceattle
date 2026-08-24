# Provenance: adversarial review of PR #117 (Tier 0 cleanup), defect 2.
#
# comp_data$Sex is an encoding -- 0 combined, 1 female, 2 male, 3 joint female
# and male -- and until 5.13.0 nothing checked it against nsex, though M1_base,
# weight and ration_data are all checked.
#
# The joint case is the damaging one, and it is a disagreement between two
# registries rather than a single wrong line. check_composition_data()
# (R/5-rearrange_data.R) computes its joint_adjust as
# `nsex[Species] == 2 & Sex == 3`, so for a one-sex species it requires only
# nages composition columns. The template keys on `flt_sex == 3` alone, so it
# writes that row's predicted composition out to 2 * nages. The surplus lands in
# the next observation's row of comp_hat -- a silently wrong likelihood for an
# observation that looked fine, not a crash. Sex = 2 on a one-sex species is the
# same class: the template reads F_flt_age at sex index 1, which does not exist.
#
# Refused at the boundary rather than reconciled in the template: a joint
# composition for a species the model runs with one sex is not a configuration
# with a right answer.

test_that("data_check() refuses joint or male-only comps on a one-sex species", {
  d <- make_test_data(nyrs = 8, nages = 5, seed = 1)
  expect_true(all(d$nsex == 1))

  for (bad in c(2, 3)) {
    dd <- d
    dd$comp_data$Sex[1] <- bad
    expect_error(suppressMessages(suppressWarnings(data_check(dd))),
                 "one-sex species", info = paste("Sex =", bad))
  }
})


test_that("data_check() allows the sex codes a one-sex species can carry", {
  d <- make_test_data(nyrs = 8, nages = 5, seed = 1)

  # 0 (combined) is the normal case; 1 (female only) resolves to sex index 0,
  # which exists, so it is left alone rather than second-guessed.
  for (ok in c(0, 1)) {
    dd <- d
    dd$comp_data$Sex[1] <- ok
    expect_no_error(suppressMessages(suppressWarnings(data_check(dd))))
  }
})


test_that("a two-sex species may carry a joint composition", {
  # The check is on the pairing, not on the code, so it must not fire when the
  # second sex exists. GOA2018SS is the bundled joint-sex case: every Sex = 3
  # row belongs to a species the model runs with two sexes.
  expect_true(any(GOA2018SS$comp_data$Sex == 3))
  joint_spp <- unique(GOA2018SS$comp_data$Species[GOA2018SS$comp_data$Sex == 3])
  expect_true(all(GOA2018SS$nsex[joint_spp] == 2))

  # data_check() on the bundled object raises unrelated errors -- it is written
  # for a data_list that has been through the fit_mod() pipeline -- so pin that
  # THIS check is not among them rather than that nothing is raised.
  err <- tryCatch({ suppressMessages(suppressWarnings(data_check(GOA2018SS))); NULL },
                  error = conditionMessage)
  expect_false(isTRUE(grepl("one-sex species", err)))

  # And the same species with nsex forced to 1 does trip it, so the pairing is
  # what the check reads.
  bad <- GOA2018SS
  bad$nsex[joint_spp] <- 1
  err_bad <- tryCatch({ suppressMessages(suppressWarnings(data_check(bad))); NULL },
                      error = conditionMessage)
  expect_true(isTRUE(grepl("one-sex species", err_bad)))
})
