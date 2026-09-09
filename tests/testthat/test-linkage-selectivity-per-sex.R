# A per-sex selectivity linkage masks the base parameter for that sex only.
#
# Provenance: the `sel` branch of map_linkage_adjuster() masked every sex, so a
# linkage stratified on one sex also fixed the other's inflection or slope --
# the reference sex a Stock Synthesis-style offset is measured from, which has
# to stay estimated for the offset to mean anything.
#
# GOAatf is the fixture: the only bundled two-sex fit that converges
# (inst/dev/TRAPS.md). estimateMode = 3 builds the map without optimizing.

sel_map_for <- function(spec, flt = 3L) {
  # GOAatf ships this fishery as NonParametric, which carries no linkage; the
  # wired forms are Logistic / DoubleLogistic / DescendingLogistic /
  # DoubleNormal / LogisticPM.
  d <- Rceattle::GOAatf
  d$fleet_control$Selectivity[flt] <- "Logistic"
  # GOAatf ships Selectivity_index 2 on Fleet_code 3, which reads as a mirror;
  # a prior is refused on a mirrored fleet because the shared block would be
  # penalized once per sharing fleet.
  d$fleet_control$Selectivity_index[flt] <- flt
  fit <- suppressMessages(suppressWarnings(Rceattle::fit_mod(
    data_list = d, msmMode = 0, estimateMode = 3,
    selFun = spec,
    fit_control = Rceattle::fit_control(getsd = FALSE, verbose = 0))))
  # mapList keeps the [slot, fleet, sex] shape; mapFactor is flattened.
  fit$map$mapList$sel_inf
}

testthat::test_that("a per-sex selectivity linkage leaves the other sex estimated", {
  testthat::skip_on_cran()

  flt <- 3L   # GOA_atf_fishery, the two-sex fishery
  # est_phase = 0 is the fix-at-init contract: the intercept holds the level
  # and the base parameter is mapped out. That is the branch under test.
  spec <- Rceattle::build_selectivity(linkages = list(
    inf_asc = Rceattle::linkage_spec(~ 1, by = ~ fleet + sex,
                                     fleet = flt, sex = 2L,
                                     est_phase = 0L)))
  m <- sel_map_for(spec)

  # Slot 1 is the ascending inflection. The linked sex is fixed; the reference
  # sex is not -- that is the whole point of an offset parameterization.
  testthat::expect_true(is.na(m[1, flt, 2]))
  testthat::expect_false(is.na(m[1, flt, 1]))
})


testthat::test_that("a linkage with no sex stratum still masks both sexes", {
  testthat::skip_on_cran()

  # The unstratified case must not change: `by = ~ fleet` means the offset
  # applies to the whole fleet, so both sexes' base parameters are replaced.
  flt <- 3L
  spec <- Rceattle::build_selectivity(linkages = list(
    inf_asc = Rceattle::linkage_spec(~ 1, by = ~ fleet, fleet = flt,
                                     est_phase = 0L)))
  m <- sel_map_for(spec)

  testthat::expect_true(all(is.na(m[1, flt, ])))
})
