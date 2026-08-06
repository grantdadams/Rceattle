# Round-trip fidelity of the run-config <-> plain-list converters that back
# save_config()/load_config(). Everything serializes as plain scalars/lists
# except the linkage `formula` objects and `Rceattle_prior` S3 objects, so the
# tests exercise a config carrying all of those.

# A run config exercising every serialization-hostile case: a selectivity
# linkage with a fixed + random-effect formula, a non-default `by`, a prior,
# init, bounds, and an estimation phase; a non-default HCR; non-default modes;
# non-default estimation controls and fit_control.
.rich_run_config <- function() {
  sel <- Rceattle::build_selectivity(linkages = list(
    inf_asc = Rceattle::linkage_spec(
      ~ cold_pool + (1 | Year), by = ~ fleet,
      priors = list("(Intercept)" = Rceattle::prior_normal(0, 3)),
      init   = list("(Intercept)" = 0.5),
      bounds = list("(Intercept)" = c(-5, 5)),
      est_phase = 2L)))
  mc <- Rceattle::model_config(
    msmMode = 1, initMode = "OffsetEquilibrium",
    HCR = Rceattle::build_hcr(HCR = "NPFMC", Ptarget = 0.4),
    selFun = sel)
  Rceattle:::.rce_run_config(
    mc = mc, estimateMode = 1, random_rec = TRUE,
    fc = Rceattle::fit_control(getsd = FALSE, phase = TRUE, loopnum = 1))
}

testthat::test_that("run-config to-list / from-list is a round-trip fixed point", {
  testthat::skip_if_not_installed("Rceattle")
  rc <- .rich_run_config()
  l1 <- Rceattle:::.rce_run_config_to_list(rc)
  # Rebuilding from the list and re-serializing must reproduce the list exactly.
  l2 <- Rceattle:::.rce_run_config_to_list(
    Rceattle:::.rce_run_config_from_list(l1))
  testthat::expect_identical(l1, l2)
})

testthat::test_that("a real YAML file round-trips the config losslessly", {
  testthat::skip_if_not_installed("yaml")
  testthat::skip_if_not_installed("Rceattle")
  rc <- .rich_run_config()
  f <- withr::local_tempfile(fileext = ".yaml")
  yaml::write_yaml(Rceattle:::.rce_run_config_to_list(rc), f)
  rc2 <- Rceattle:::.rce_run_config_from_list(yaml::read_yaml(f))

  # Scalars.
  testthat::expect_equal(rc2$model_config$msmMode, 1)
  testthat::expect_equal(rc2$model_config$initMode, "OffsetEquilibrium")
  testthat::expect_equal(rc2$estimateMode, 1)
  testthat::expect_true(rc2$random_rec)
  testthat::expect_false(rc2$fit_control$getsd)
  testthat::expect_true(rc2$fit_control$phase)
  testthat::expect_equal(rc2$fit_control$loopnum, 1)
  testthat::expect_equal(rc2$model_config$HCR$HCR, "NPFMC")

  # The linkage spec: formula/by/prior/init/bounds/phase survive.
  sp <- rc2$model_config$selFun$linkages$inf_asc
  testthat::expect_equal(deparse(sp$formula), "~cold_pool + (1 | Year)")
  testthat::expect_equal(deparse(sp$by), "~fleet")
  testthat::expect_equal(sp$param, "inf_asc")            # re-stamped from the key
  testthat::expect_equal(sp$priors[["(Intercept)"]]$family, "normal")
  testthat::expect_equal(sp$priors[["(Intercept)"]]$p1, 0)
  testthat::expect_equal(sp$priors[["(Intercept)"]]$p2, 3)
  testthat::expect_equal(sp$init[["(Intercept)"]], 0.5)
  testthat::expect_equal(sp$bounds[["(Intercept)"]], c(-5, 5))
  testthat::expect_equal(sp$est_phase, 2L)
})

testthat::test_that("a default run config serializes to an (almost) empty list", {
  testthat::skip_if_not_installed("Rceattle")
  rc <- Rceattle:::.rce_run_config()   # all defaults
  l <- Rceattle:::.rce_run_config_to_list(rc)
  # Default-omission: nothing to record when everything is at its default.
  testthat::expect_length(l, 0)
  # And it still rebuilds to a valid default run config.
  rc2 <- Rceattle:::.rce_run_config_from_list(l)
  testthat::expect_s3_class(rc2, "Rceattle_run_config")
  testthat::expect_equal(rc2$model_config$msmMode, 0)
})

testthat::test_that("save_config writes a documented YAML that load_config rebuilds", {
  testthat::skip_if_not_installed("yaml")
  testthat::skip_if_not_installed("Rceattle")
  rc <- .rich_run_config()
  f <- withr::local_tempfile(fileext = ".yaml")
  ret <- Rceattle::save_config(rc, f)
  testthat::expect_s3_class(ret, "Rceattle_run_config")
  testthat::expect_true(file.exists(f))

  lines <- readLines(f)
  # Provenance header + at least one doc comment.
  testthat::expect_true(any(grepl("^# Rceattle run configuration", lines)))
  testthat::expect_true(any(grepl("Predation-mortality mode", lines)))

  rc2 <- Rceattle::load_config(f)
  testthat::expect_identical(Rceattle:::.rce_run_config_to_list(rc),
                             Rceattle:::.rce_run_config_to_list(rc2))
})

testthat::test_that("run_config() normalizes a model_config, run_config, and data list", {
  testthat::skip_if_not_installed("Rceattle")
  mc <- Rceattle::model_config(msmMode = 1)
  testthat::expect_equal(Rceattle::run_config(mc)$model_config$msmMode, 1)

  rc <- Rceattle:::.rce_run_config(mc = mc, estimateMode = 2)
  testthat::expect_identical(Rceattle::run_config(rc), rc)

  dl <- list(model_config = mc, estimateMode = 3)
  cfg <- Rceattle::run_config(dl)
  testthat::expect_equal(cfg$model_config$msmMode, 1)
  testthat::expect_equal(cfg$estimateMode, 3)

  # ... overrides win.
  testthat::expect_equal(Rceattle::run_config(mc, estimateMode = 9)$estimateMode, 9)
})

testthat::test_that("load_config errors clearly on a missing file", {
  testthat::expect_error(Rceattle::load_config(tempfile()), "not found")
})

testthat::test_that("a shared-coefficient linkage (by = NULL) round-trips as NULL, not ~species", {
  testthat::skip_if_not_installed("yaml")
  testthat::skip_if_not_installed("Rceattle")
  # by = NULL means one coefficient set shared across all species -- a distinct
  # model from the default by = ~species (one per species).
  q <- Rceattle::build_catchability(linkages = list(
    q = Rceattle::linkage_spec(~ cold_pool, by = NULL)))
  rc <- Rceattle:::.rce_run_config(mc = Rceattle::model_config(qFun = q))
  f <- withr::local_tempfile(fileext = ".yaml")
  Rceattle::save_config(rc, f)
  sp <- Rceattle::load_config(f)$model_config$qFun$linkages$q
  testthat::expect_null(sp$by)                 # shared, NOT ~species

  # And a spec with no priors does not emit an empty `priors:` field.
  l <- Rceattle:::.rce_run_config_to_list(rc)
  testthat::expect_false("priors" %in% names(l$model$qFun$linkages$q))
})

testthat::test_that("an explicit by = ~ species on a fleet process survives the round-trip", {
  testthat::skip_if_not_installed("yaml")
  testthat::skip_if_not_installed("Rceattle")
  # The serializer keys default-omission on `by_auto`, not on equality with
  # ~species: an *explicit* by = ~ species on catchability (whose auto default is
  # ~ fleet) must not be dropped and silently re-resolved to ~ fleet on reload.
  q <- Rceattle::build_catchability(linkages = list(
    q = Rceattle::linkage_spec(~ cold_pool, by = ~ species)))
  rc <- Rceattle:::.rce_run_config(mc = Rceattle::model_config(qFun = q))
  f <- withr::local_tempfile(fileext = ".yaml")
  Rceattle::save_config(rc, f)
  sp <- Rceattle::load_config(f)$model_config$qFun$linkages$q
  testthat::expect_identical(all.vars(sp$by), "species")   # kept, NOT re-resolved to ~fleet

  # while an *omitted* by on the same fleet process reloads as the ~fleet default
  q2 <- Rceattle::build_catchability(linkages = list(
    q = Rceattle::linkage_spec(~ cold_pool)))
  rc2 <- Rceattle:::.rce_run_config(mc = Rceattle::model_config(qFun = q2))
  f2 <- withr::local_tempfile(fileext = ".yaml")
  Rceattle::save_config(rc2, f2)
  sp2 <- Rceattle::load_config(f2)$model_config$qFun$linkages$q
  testthat::expect_identical(all.vars(sp2$by), "fleet")
})
