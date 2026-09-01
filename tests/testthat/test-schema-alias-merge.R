# A deprecated name present ALONGSIDE its canonical one.
#
# Provenance. `read_data()` canonicalises `fleet_control`, so a script that then
# assigns a deprecated name creates the old column beside the new one. Until
# 5.25.0 the upgrade kept the canonical column and deleted the other, silently
# discarding whatever the caller had just written behind a routine deprecation
# message. In `Rceattle-models` that hid a fishery selectivity deviation sd of
# 10.3 where the script asked for 0.35, a survey catchability prior that never
# took effect, `Proj_F_proportion = 0` on three GOA circulation-study
# assessments whose scripts set 1, and a "Dirichlet-multinomial" hake analysis
# that was running full multinomial.
#
# Nothing is merged, deliberately. `NA` is a real setting here -- `Sel_norm_bin`
# and `Sel_cap_bin` mean "do not normalize" / "no cap", `Proj_F_proportion`
# means "no F apportioned" -- so filling the canonical column's blanks from the
# deprecated one cannot express clearing a value, and an earlier draft of this
# change did exactly that, silently moving a GOA pollock bridging fleet from
# unnormalized to normalized at bin 3. Both names present must agree.
testthat::skip_on_cran()

testthat::test_that("the ordinary upgrade is unchanged", {
  testthat::skip_if_not_installed("Rceattle")
  data("GOAcod")

  d <- GOAcod
  d$fleet_control$Catchability_prior_sd <- NULL
  d$fleet_control$Q_sd_prior <- c(0.2, NA, NA, NA, NA)
  x <- suppressMessages(suppressWarnings(Rceattle:::switch_check(d)))
  testthat::expect_equal(x$fleet_control$Catchability_prior_sd[1], 0.2)
  testthat::expect_null(x$fleet_control$Q_sd_prior)
})

testthat::test_that("two spellings that agree are not a conflict", {
  testthat::skip_if_not_installed("Rceattle")
  data("GOAcod")

  # An integer code and its string name one setting.
  d <- GOAcod
  d$fleet_control$Comp_distribution <- "Multinomial"
  d$fleet_control$Comp_loglike      <- 0
  x <- suppressMessages(suppressWarnings(Rceattle:::switch_check(d)))
  testthat::expect_equal(unname(x$fleet_control$Comp_distribution[1]), "Multinomial")

  # Identical values under both names.
  d2 <- GOAcod
  d2$fleet_control$Catchability_prior_sd <- c(0.2, NA, NA, NA, NA)
  d2$fleet_control$Q_sd_prior            <- c(0.2, NA, NA, NA, NA)
  testthat::expect_equal(
    suppressMessages(suppressWarnings(
      Rceattle:::switch_check(d2)))$fleet_control$Catchability_prior_sd[1], 0.2)
})

testthat::test_that("a disagreement is an error, and nothing is merged", {
  testthat::skip_if_not_installed("Rceattle")
  data("GOAcod")

  # Differing values.
  d <- GOAcod
  d$fleet_control$Catchability_prior_sd <- c(0.9, NA, NA, NA, NA)
  d$fleet_control$Q_sd_prior            <- c(0.2, NA, NA, NA, NA)
  testthat::expect_error(
    suppressMessages(suppressWarnings(Rceattle:::switch_check(d))), "differ")

  # NA against a value. This is the case that must NOT be filled: NA is
  # "do not normalize", so taking the other column would rescale selectivity.
  d2 <- GOAcod
  d2$fleet_control$Sel_norm_bin     <- c(3, 10, 10, NA, NA)
  d2$fleet_control$Age_max_selected <- c(NA, NA, NA, NA, 3)
  testthat::expect_error(
    suppressMessages(suppressWarnings(Rceattle:::switch_check(d2))), "differ")

  # Differing lengths, which a fill would have collapsed to length 1 and sent
  # to a model built without bounds checking.
  d3 <- GOAcod
  d3$sigma_rec       <- c(NA, NA, NA)
  d3$sigma_rec_prior <- 0.7
  testthat::expect_error(
    suppressMessages(suppressWarnings(Rceattle:::switch_check(d3))), "are different lengths")

  # A value difference is caught even when the two columns are typed
  # differently; the type itself is not the disagreement.
  d4 <- GOAcod
  d4$fleet_control$N_sel_bins <- as.integer(c(1, 2, 3, 4, 5))
  d4$fleet_control$Nselages   <- c("1", "2", "3", "4", "9")
  testthat::expect_error(
    suppressMessages(suppressWarnings(Rceattle:::switch_check(d4))), "differ")

  # The same values written as integer and as character ARE the same setting,
  # and the canonical column keeps its own type.
  d5 <- GOAcod
  d5$fleet_control$N_sel_bins <- as.integer(c(1, 2, 3, 4, 5))
  d5$fleet_control$Nselages   <- c("1", "2", "3", "4", "5")
  x5 <- suppressMessages(suppressWarnings(Rceattle:::switch_check(d5)))
  testthat::expect_equal(x5$fleet_control$N_sel_bins, as.integer(1:5))
})

testthat::test_that("the error names the fleet, and the column the value came from", {
  testthat::skip_if_not_installed("Rceattle")
  data("GOA2018SS")

  # A wide table whose only disagreement is in the tail: the message has to name
  # the position, not print the first eight values of each and look identical.
  d <- GOA2018SS
  d$fleet_control$Catchability_prior_sd <- rep(0.2, 16)
  d$fleet_control$Q_sd_prior            <- c(rep(0.2, 15), 0.9)
  testthat::expect_error(
    suppressMessages(suppressWarnings(Rceattle:::switch_check(d))),
    d$fleet_control$Fleet_name[16], fixed = TRUE)

  # Where a first alias supplied the value, the message says so rather than
  # naming a canonical column the caller never wrote.
  data("GOAcod")
  d2 <- GOAcod
  d2$fleet_control$Sel_norm_bin     <- NULL
  d2$fleet_control$Age_max_selected <- rep(3, 5)
  d2$fleet_control$Sel_norm_bin1    <- rep(7, 5)
  testthat::expect_error(
    suppressMessages(suppressWarnings(Rceattle:::switch_check(d2))),
    "taken from 'Age_max_selected'", fixed = TRUE)
})

testthat::test_that("an unrecognized column that looks like a real one is reported", {
  testthat::skip_if_not_installed("Rceattle")
  data("GOAcod")

  # `fleet_control` is a data frame, so a typo creates a column and sets nothing.
  d <- suppressMessages(suppressWarnings(Rceattle:::switch_check(GOAcod)))
  d$fleet_control$Catchabilty_init <- 0.5
  # Collected rather than expect_warning()-ed: data_check() raises other,
  # unrelated warnings on this dataset, and suppressing them would suppress this.
  seen <- character(0)
  suppressMessages(withCallingHandlers(
    try(Rceattle:::data_check(d), silent = TRUE),
    warning = function(x) { seen <<- c(seen, conditionMessage(x))
                            invokeRestart("muffleWarning") }))
  testthat::expect_true(any(grepl("Catchability_init", seen, fixed = TRUE)))
  testthat::expect_true(any(grepl("did you mean", seen, fixed = TRUE)))

  # Assessment workbooks legitimately carry columns Rceattle does not read;
  # warning about those would be noise. 147 of the 183 fleet_control sheets in
  # the sibling repositories have this one.
  d2 <- suppressMessages(suppressWarnings(Rceattle:::switch_check(GOAcod)))
  d2$fleet_control$Accumatation_age_lower <- 3
  saw <- NULL
  suppressMessages(withCallingHandlers(
    try(suppressWarnings(Rceattle:::data_check(d2)), silent = TRUE),
    warning = function(x) {
      if (grepl("did you mean", conditionMessage(x))) saw <<- conditionMessage(x)
      invokeRestart("muffleWarning")
    }))
  testthat::expect_null(saw)
})

testthat::test_that("the pre-4.4 composition-accumulation columns are read", {
  testthat::skip_if_not_installed("Rceattle")
  data("GOAcod")

  # 183 of the fleet_control sheets in the sibling repositories carry these two
  # names; 79 carry a young-tail value that folds bins. They were absent from the
  # alias list, so every one of those settings was silently defaulted away.
  d <- GOAcod
  d$fleet_control$Comp_accum_young <- NULL
  d$fleet_control$Comp_accum_old   <- NULL
  d$fleet_control$Accumatation_age_lower <- c(3, NA, NA, NA, 2)
  d$fleet_control$Accumatation_age_upper <- c(10, NA, NA, NA, 10)
  x <- suppressMessages(suppressWarnings(Rceattle:::switch_check(d)))
  testthat::expect_equal(x$fleet_control$Comp_accum_young, c(3, NA, NA, NA, 2))
  testthat::expect_equal(x$fleet_control$Comp_accum_old,   c(10, NA, NA, NA, 10))
  testthat::expect_null(x$fleet_control$Accumatation_age_lower)

  # They are 1-based ordinals on the composition dimension, which equals the age
  # only when minage is 1 -- true of every workbook carrying the old names.
  testthat::expect_equal(unname(GOAcod$minage), 1)
})

testthat::test_that("Est_weights_mcallister stays an output, not an input", {
  testthat::skip_if_not_installed("Rceattle")

  # It is the pre-4.12 name of a diagnostic fit_mod() COMPUTES, mirrored on the
  # fitted object for back-compatibility. Treating a workbook copy of it as an
  # input would imply an effect on the fit that it does not have.
  sch <- Rceattle:::.rce_column_schema()
  aliases <- unlist(lapply(sch, function(r) r$aliases))
  testthat::expect_false("Est_weights_mcallister" %in% aliases)
  testthat::expect_false("Comp_weights_mcallister" %in%
                           vapply(sch, function(r) r$name, character(1)))
})
