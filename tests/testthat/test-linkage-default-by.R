# An omitted `by` in linkage_spec() resolves to the base stratum of whichever
# process the spec is attached to: ~fleet for catchability / selectivity / the
# fleet composition weights (theta_comp, theta_caal), ~species for recruitment,
# M, growth, and the diet weight (theta_diet). An explicit `by` -- a formula or
# NULL (shared) -- is always kept as the user gave it. (Fast: no model fit.)

testthat::test_that("omitted by resolves to ~fleet for the fleet-keyed processes", {
  q  <- Rceattle::build_catchability(linkages = list(
          q = Rceattle::linkage_spec(~ temp)))$linkages$q
  se <- Rceattle::build_selectivity(linkages = list(
          inf_asc = Rceattle::linkage_spec(~ temp)))$linkages$inf_asc
  testthat::expect_identical(all.vars(q$by),  "fleet")
  testthat::expect_identical(all.vars(se$by), "fleet")
})

testthat::test_that("omitted by resolves to ~species for the species-keyed processes", {
  r <- Rceattle::build_srr(srr_fun = "BevertonHolt", linkages = list(
          R0 = Rceattle::linkage_spec(~ temp)))$linkages$R0
  m <- Rceattle::build_M1(linkages = list(
          M1 = Rceattle::linkage_spec(~ temp)))$linkages$M1
  g <- Rceattle::build_growth(linkages = list(
          Linf = Rceattle::linkage_spec(~ temp)))$linkages$Linf
  testthat::expect_identical(all.vars(r$by), "species")
  testthat::expect_identical(all.vars(m$by), "species")
  testthat::expect_identical(all.vars(g$by), "species")
})

testthat::test_that("composition's default by splits by parameter (fleet weights vs diet)", {
  tc <- Rceattle::build_composition(linkages = list(
          theta_comp = Rceattle::linkage_spec(~ 1,
            priors = list(`(Intercept)` = gamma(2, 0.5)))))$linkages$theta_comp
  td <- Rceattle::build_composition(linkages = list(
          theta_diet = Rceattle::linkage_spec(~ 1,
            priors = list(`(Intercept)` = gamma(2, 0.5)))))$linkages$theta_diet
  testthat::expect_identical(all.vars(tc$by), "fleet")     # fleet DM weight
  testthat::expect_identical(all.vars(td$by), "species")   # per-predator diet weight
})

testthat::test_that("an explicit by (including NULL) is always kept, not auto-filled", {
  q_sp   <- Rceattle::build_catchability(linkages = list(
              q = Rceattle::linkage_spec(~ temp, by = ~ species)))$linkages$q
  q_null <- Rceattle::build_catchability(linkages = list(
              q = Rceattle::linkage_spec(~ temp, by = NULL)))$linkages$q
  testthat::expect_identical(all.vars(q_sp$by), "species")  # not overridden to ~fleet
  testthat::expect_null(q_null$by)                          # shared coefficient kept
})

testthat::test_that("auto-by matches an explicit by = ~fleet end to end", {
  auto <- Rceattle::build_catchability(linkages = list(
            q = Rceattle::linkage_spec(~ temp, fleet = 1L)))$linkages$q
  expl <- Rceattle::build_catchability(linkages = list(
            q = Rceattle::linkage_spec(~ temp, by = ~ fleet, fleet = 1L)))$linkages$q
  auto$by_auto <- expl$by_auto <- NULL   # the only intended difference
  testthat::expect_equal(all.vars(auto$by), all.vars(expl$by))
  testthat::expect_equal(auto[setdiff(names(auto), "by")],
                         expl[setdiff(names(expl), "by")])
})
