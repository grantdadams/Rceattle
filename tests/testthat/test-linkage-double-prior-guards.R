# A prior counted twice (5.33.0): Catchability = "Estimated-with-prior" and M1_use_prior
# penalize the same log q / log M1 as a linkage intercept prior, as the Ricker alpha pair did.

.prior_row <- function(process, fleet = NA_integer_, species = NA_integer_,
                       prior_family = "lognormal", design_col = "(Intercept)") {
  data.frame(process = process, fleet = fleet, species = species,
             design_col = design_col, prior_family = prior_family,
             stringsAsFactors = FALSE)
}

testthat::test_that("a q intercept prior on an Estimated-with-prior fleet is refused", {
  fc  <- Rceattle::switch_check(Rceattle::clean_data(Rceattle::BS2017SS))$fleet_control
  flt <- which(fc$Fleet_name == "EIT_Pollock")
  fc$Catchability[flt] <- "Estimated-with-prior"
  testthat::expect_error(
    Rceattle:::.check_q_linkage_support(.prior_row("q", fleet = flt), fc),
    "counted twice")

  # The integer code is refused too.
  fc2 <- fc; fc2$Catchability[flt] <- 2
  testthat::expect_error(
    Rceattle:::.check_q_linkage_support(.prior_row("q", fleet = flt), fc2),
    "counted twice")

  # Not a double count: q estimated without a prior, a prior on a slope, or no prior.
  fc3 <- fc; fc3$Catchability[flt] <- "Estimated"
  testthat::expect_silent(
    Rceattle:::.check_q_linkage_support(.prior_row("q", fleet = flt), fc3))
  testthat::expect_silent(Rceattle:::.check_q_linkage_support(
    .prior_row("q", fleet = flt, design_col = "temp"), fc))
  testthat::expect_silent(Rceattle:::.check_q_linkage_support(
    .prior_row("q", fleet = flt, prior_family = "none"), fc))
})

testthat::test_that("the q guard reads a shared Catchability_index group's lead fleet", {
  # Fleets sharing a Catchability_index share one q; the template priors it on the
  # group's first estimated fleet, so that fleet's Catchability decides.
  fc   <- Rceattle::switch_check(Rceattle::clean_data(Rceattle::BS2017SS))$fleet_control
  lead <- which(fc$Fleet_name == "BT_Pollock")
  flw  <- which(fc$Fleet_name == "EIT_Pollock")
  fc$Catchability_index[flw] <- fc$Catchability_index[lead]

  # Lead carries the q prior, and the intercept prior sits on the follower.
  fc$Catchability[lead] <- "Estimated-with-prior"; fc$Catchability[flw] <- "Estimated"
  testthat::expect_error(
    Rceattle:::.check_q_linkage_support(.prior_row("q", fleet = flw), fc),
    "led by BT_Pollock.*counted twice")

  # The group has no q prior, whatever the follower's own Catchability says.
  fc$Catchability[lead] <- "Estimated"; fc$Catchability[flw] <- "Estimated-with-prior"
  testthat::expect_silent(
    Rceattle:::.check_q_linkage_support(.prior_row("q", fleet = flw), fc))
})

testthat::test_that("an M1 intercept prior with M1_use_prior is refused", {
  sp3 <- c("a", "b", "c")
  testthat::expect_error(Rceattle:::.check_M_linkage_prior(
    .prior_row("M", species = 2L), c(0, 1, 0), c(0, 0, 0), sp3), "species b:")

  # A row with no species targets species 1, as in the template.
  testthat::expect_error(Rceattle:::.check_M_linkage_prior(
    .prior_row("M"), c(1, 0, 0), c(0, 0, 0), sp3), "species a:")

  # Not a double count: the M prior on another species, on total M (M2_use_prior),
  # or no prior on the intercept.
  testthat::expect_silent(Rceattle:::.check_M_linkage_prior(
    .prior_row("M", species = 2L), c(1, 0, 1), c(0, 0, 0), sp3))
  testthat::expect_silent(Rceattle:::.check_M_linkage_prior(
    .prior_row("M", species = 2L), c(0, 1, 0), c(0, 1, 0), sp3))
  testthat::expect_silent(Rceattle:::.check_M_linkage_prior(
    .prior_row("M", species = 2L, prior_family = "none"), c(0, 1, 0), c(0, 0, 0), sp3))
})

testthat::test_that("fit_mod() runs both guards on the linkage table it builds", {
  d <- Rceattle::BS2017SS
  M_prior_link <- Rceattle::linkage_spec(~ 1, priors = list(
    `(Intercept)` = Rceattle::prior_lognormal(log(0.2), 0.1)))
  testthat::expect_error(suppressWarnings(suppressMessages(Rceattle::fit_mod(
    data_list = d, msmMode = 0, estimateMode = 3, verbose = 0,
    M1Fun = Rceattle::build_M1(M1_model = 1, M1_use_prior = TRUE, M_prior = 0.2,
                               M_prior_sd = 0.1, linkages = list(M1 = M_prior_link))))),
    "counted twice")

  flt <- which(d$fleet_control$Fleet_name == "EIT_Pollock")
  d$fleet_control$Catchability[flt]          <- 2
  d$fleet_control$Catchability_prior_sd[flt] <- 0.2
  q_prior_link <- Rceattle::linkage_spec(~ 1, fleet = flt, priors = list(
    `(Intercept)` = Rceattle::prior_lognormal(0, 0.2)))
  testthat::expect_error(suppressWarnings(suppressMessages(Rceattle::fit_mod(
    data_list = d, msmMode = 0, estimateMode = 3, verbose = 0,
    qFun = Rceattle::build_catchability(linkages = list(q = q_prior_link))))),
    "counted twice")
})
