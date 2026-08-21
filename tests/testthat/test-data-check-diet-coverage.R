# Under empirical suitability (suitMode = 0) the suitability array is read
# straight out of diet_data, so an age with no diet row is switched off rather
# than estimated: suit_other goes to 1 for a predator age with no rows, and a
# prey age with no rows is never eaten. Neither raises an error or moves the
# likelihood, so these tests pin that data_check() says so.
#
# The subtlety the check must respect is that Pred_age / Prey_age below the
# species' minage are not ages -- they select the aggregated diet formats, which
# organize_diet_obs() skips when filling the suitability array. Counting them as
# coverage would hide the gap this warns about.

prep_ms <- function(x) {
  x$HCR <- 0                    # build_hcr() would set this
  x$msmMode <- 1                # fit_mod() copies its arguments onto the
  x$suitMode <- rep(0, x$nspp)  # data_list before calling data_check()
  suppressMessages(suppressWarnings(
    Rceattle::switch_check(Rceattle::clean_data(x))))
}

quiet_check <- function(x) suppressMessages(data_check(x))


testthat::test_that("full age coverage passes without a diet-coverage warning", {
  testthat::skip_on_cran()

  data("BS2017MS", package = "Rceattle", envir = environment())
  cleaned <- prep_ms(BS2017MS)

  # Every age of every species appears as both predator and prey, so there is
  # nothing to warn about. Guard the premise so a change to the bundled data
  # cannot quietly turn this into a vacuous test.
  dd <- cleaned$diet_data
  testthat::expect_true(all(seq_len(cleaned$nages[3]) %in%
                            dd$Pred_age[dd$Pred == 3]))

  testthat::expect_no_warning(quiet_check(cleaned))
})


testthat::test_that("data_check() warns when a predator's diet data is truncated", {
  testthat::skip_on_cran()

  data("BS2017MS", package = "Rceattle", envir = environment())
  cleaned <- prep_ms(BS2017MS)

  # Arrowtooth (species 3) has 21 ages against 12 for the other two. Cut its
  # diet rows above age 10 and leave pollock and cod whole, which isolates the
  # gap the GOA data file had.
  truncated <- cleaned
  dd <- truncated$diet_data
  truncated$diet_data <- dd[!((dd$Pred == 3 & dd$Pred_age > 10) |
                              (dd$Prey == 3 & dd$Prey_age > 10)), ]
  truncated <- suppressMessages(suppressWarnings(Rceattle::clean_data(truncated)))

  testthat::expect_warning(quiet_check(truncated), "Arrowtooth flounder")
  w <- testthat::capture_warnings(quiet_check(truncated))
  cov_warn <- w[grepl("empirical suitability", w)]
  testthat::expect_length(cov_warn, 1)
  # Both roles are reported, and the missing ages are collapsed into a run.
  testthat::expect_match(cov_warn, "as PREDATOR: no diet data at age 11-21")
  testthat::expect_match(cov_warn, "as PREY: no diet data at age 11-21")
  # Species that are still fully covered are not named as gaps.
  testthat::expect_false(grepl("Pollock", cov_warn))
  testthat::expect_false(grepl("Cod", cov_warn))
})


testthat::test_that("aggregated diet rows do not count as age coverage", {
  testthat::skip_on_cran()

  data("BS2017MS", package = "Rceattle", envir = environment())
  cleaned <- prep_ms(BS2017MS)
  dd <- cleaned$diet_data

  # Replace arrowtooth's age 11+ rows with the aggregated formats: prey summed
  # across ages (Prey_age < minage), predator averaged across ages (both < minage),
  # and the abundance-weighted variant (Pred_age < -500). organize_diet_obs()
  # skips all three when building suitability, so they must not silence the
  # warning -- if they did, a model would look covered while carrying no
  # predation from those ages.
  keep <- dd[dd$Pred_age <= 10 & dd$Prey_age <= 10, ]
  agg <- keep[rep(1, 3), ]
  agg$Pred <- 3; agg$Prey <- 1
  agg$Pred_age <- c(11L, -1L, -999L)
  agg$Prey_age <- c(-1L, -1L, -1L)
  agg$Stomach_proportion_by_weight <- 0.1

  aggregated <- cleaned
  aggregated$diet_data <- rbind(keep, agg)
  aggregated <- suppressMessages(suppressWarnings(Rceattle::clean_data(aggregated)))

  w <- testthat::capture_warnings(quiet_check(aggregated))
  cov_warn <- w[grepl("empirical suitability", w)]
  testthat::expect_length(cov_warn, 1)
  # Age 11 has an aggregated predator row but still no prey-at-age-in-
  # predator-at-age row, so it is still an uncovered age.
  testthat::expect_match(cov_warn, "as PREDATOR: no diet data at age 11-21")

  # And the aggregated rows alone raise nothing when real coverage is intact.
  covered <- cleaned
  covered$diet_data <- rbind(cleaned$diet_data, agg)
  covered <- suppressMessages(suppressWarnings(Rceattle::clean_data(covered)))
  testthat::expect_no_warning(quiet_check(covered))
})


testthat::test_that("a predator left out of the diet table entirely is named as such", {
  testthat::skip_on_cran()

  data("BS2017MS", package = "Rceattle", envir = environment())
  cleaned <- prep_ms(BS2017MS)

  # Drop every row where cod (species 2) is the predator. Forgetting one
  # predator's diet reads differently from truncating an age range, so it gets
  # its own wording rather than "no diet data at age 1-12".
  dropped <- cleaned
  dropped$diet_data <- cleaned$diet_data[cleaned$diet_data$Pred != 2, ]
  dropped <- suppressMessages(suppressWarnings(Rceattle::clean_data(dropped)))

  w <- testthat::capture_warnings(quiet_check(dropped))
  cov_warn <- w[grepl("empirical suitability", w)]
  testthat::expect_length(cov_warn, 1)
  testthat::expect_match(cov_warn,
                         "Cod as PREDATOR: no diet data at any age", fixed = TRUE)
  # Cod is still eaten by the other two, so only the predator role is reported.
  testthat::expect_false(grepl("Cod as PREY", cov_warn))
})


testthat::test_that("prey coverage is checked for every species, not just empirical predators", {
  testthat::skip_on_cran()

  data("BS2017MS", package = "Rceattle", envir = environment())
  cleaned <- prep_ms(BS2017MS)

  # Cod (species 2) derives its OWN suitability parametrically, so it owes no
  # rows as a predator -- but pollock and arrowtooth still eat it empirically,
  # so it owes rows as prey. suitMode says how a species feeds, not how it is
  # fed on, so scoping the prey check to empirical species would miss this.
  #
  # Truncate cod's prey rows rather than removing them: a species with no prey
  # rows at all is one nothing eats, which is a modelling choice and is not
  # reported. A partial gap is the truncation this warns about.
  top <- max(cleaned$diet_data$Prey_age[cleaned$diet_data$Prey == 2])
  mixed <- cleaned
  mixed$suitMode <- c(0, 2, 0)
  mixed$diet_data <-
    cleaned$diet_data[!(cleaned$diet_data$Prey == 2 &
                          cleaned$diet_data$Prey_age == top), ]
  mixed <- suppressMessages(suppressWarnings(Rceattle::clean_data(mixed)))

  w <- testthat::capture_warnings(quiet_check(mixed))
  cov_warn <- w[grepl("empirical suitability", w)]
  testthat::expect_length(cov_warn, 1)
  testthat::expect_match(cov_warn, paste0("Cod as PREY: no diet data at age ", top),
                         fixed = TRUE)
  # It owes nothing as a predator, so that role is not reported for it.
  testthat::expect_false(grepl("Cod as PREDATOR", cov_warn))
})


testthat::test_that("a species nothing eats is not reported as a prey gap", {
  testthat::skip_on_cran()

  # An apex predator in a multispecies run has no prey rows at any age. That is
  # the model the author specified, not a truncated diet table, and reporting it
  # would warn on every fit of a correctly built model.
  data("BS2017MS", package = "Rceattle", envir = environment())
  cleaned <- prep_ms(BS2017MS)

  uneaten <- cleaned
  uneaten$diet_data <- cleaned$diet_data[cleaned$diet_data$Prey != 3, ]
  uneaten <- suppressMessages(suppressWarnings(Rceattle::clean_data(uneaten)))

  w <- testthat::capture_warnings(quiet_check(uneaten))
  testthat::expect_length(grep("as PREY", w), 0)

  # A predator that supplies no diet data at all is still reported: it asked for
  # empirical suitability and gets none, which it did not choose.
  nopred <- cleaned
  nopred$diet_data <- cleaned$diet_data[cleaned$diet_data$Pred != 3, ]
  nopred <- suppressMessages(suppressWarnings(Rceattle::clean_data(nopred)))

  w2 <- testthat::capture_warnings(quiet_check(nopred))
  testthat::expect_match(paste(w2, collapse = "\n"),
                         "as PREDATOR: no diet data at any age")
})


testthat::test_that("a one-sex species' diet rows count whatever the sex column says", {
  testthat::skip_on_cran()

  data("BS2017MS", package = "Rceattle", envir = environment())
  cleaned <- prep_ms(BS2017MS)
  testthat::expect_true(all(cleaned$nsex == 1))       # premise

  # organize_diet_obs() reads the sex column only when nsex == 2, so for these
  # species every row is used regardless of what Pred_sex / Prey_sex holds. The
  # check has to agree, or it reports a gap in data the model does read.
  sexed <- cleaned
  sexed$diet_data$Pred_sex <- 1L
  sexed$diet_data$Prey_sex <- 2L
  sexed <- suppressMessages(suppressWarnings(Rceattle::clean_data(sexed)))

  testthat::expect_no_warning(quiet_check(sexed))
})


testthat::test_that("the coverage warning is scoped to empirical suitability and predation", {
  testthat::skip_on_cran()

  data("BS2017MS", package = "Rceattle", envir = environment())
  cleaned <- prep_ms(BS2017MS)
  dd <- cleaned$diet_data
  truncate <- function(x) {
    x$diet_data <- dd[dd$Pred_age <= 10 & dd$Prey_age <= 10, ]
    suppressMessages(suppressWarnings(Rceattle::clean_data(x)))
  }

  # A parametric-suitability predator does not read the diet data for its
  # suitability, so a gap in it costs that predator nothing.
  parametric <- cleaned
  parametric$suitMode <- rep(2, parametric$nspp)
  parametric <- truncate(parametric)
  w <- testthat::capture_warnings(quiet_check(parametric))
  testthat::expect_length(w[grepl("empirical suitability", w)], 0)

  # Neither does a single-species model, which never forms suitability at all.
  single <- cleaned
  single$msmMode <- 0
  single <- truncate(single)
  w <- testthat::capture_warnings(quiet_check(single))
  testthat::expect_length(w[grepl("empirical suitability", w)], 0)
})


testthat::test_that("age coverage follows minage, not the age count", {
  testthat::skip_on_cran()

  # `nages` counts age bins; the ages run minage .. minage+nages-1. Writing the
  # expected set as `minage:nages` coincides with that only at minage = 1, which
  # is what both bundled multispecies datasets use -- so every other test here
  # passes either way. minage = 0 is supported (see test-growth-minage0.R).
  data("BS2017MS", package = "Rceattle", envir = environment())
  cleaned <- prep_ms(BS2017MS)
  testthat::expect_true(all(cleaned$minage == 1))     # premise
  testthat::expect_no_warning(quiet_check(cleaned))

  # Shift species 1 to minage = 0 and move its diet rows down to match, so its
  # coverage is still complete: ages 0..nages-1, all present.
  shifted <- cleaned
  shifted$minage[1] <- 0L
  dd <- shifted$diet_data
  i <- dd$Pred == 1 & dd$Pred_age >= 1
  dd$Pred_age[i] <- dd$Pred_age[i] - 1L
  j <- dd$Prey == 1 & dd$Prey_age >= 1
  dd$Prey_age[j] <- dd$Prey_age[j] - 1L
  shifted$diet_data <- dd

  # minage:nages would demand a row at age `nages`, one past the plus group,
  # and report it as an uncovered age that exerts no predation.
  w <- testthat::capture_warnings(
    Rceattle:::.check_diet_age_coverage(shifted, shifted$diet_data))
  testthat::expect_length(grep("does not cover every age", w), 0)

  # And the gap that IS there is still found: drop the real plus group.
  holed <- shifted
  holed$diet_data <- dd[!(dd$Pred == 1 & dd$Pred_age == shifted$nages[1] - 1L), ]
  w2 <- testthat::capture_warnings(
    Rceattle:::.check_diet_age_coverage(holed, holed$diet_data))
  cov2 <- w2[grepl("does not cover every age", w2)]
  testthat::expect_length(cov2, 1)
  testthat::expect_match(cov2, paste0("as PREDATOR: no diet data at age ",
                                      shifted$nages[1] - 1L), fixed = TRUE)
})
