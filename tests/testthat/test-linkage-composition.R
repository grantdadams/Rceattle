# =============================================================================
# Composition-weighting (comp) linkages: attach a prior to a Dirichlet-
# multinomial (DM) overdispersion weight via the linkage grammar.
#
# `build_composition()` exposes the DM weights as a "prior-only" process on the
# linkage table: an intercept row (`~ 1`) whose (Intercept) coefficient is
# mapped out at 0 and whose prior re-targets the underlying weight parameter
#   theta_comp  -> comp_weights        (age/length comps, Comp_distribution)
#   theta_caal  -> caal_weights        (conditional age-at-length, CAAL_distribution)
#   theta_diet  -> diet_comp_weights   (predator diet, Diet_distribution)
# The DM parameter is exp(weight), so the intercept prior is evaluated on the
# NATURAL scale: b_nat = exp(weight) (same log-link contract as the growth / M /
# q intercept priors, tested in test-linkage-intercept-prior-base-parameter.R).
#
# There is no accumulator for comp (theta is a scalar, not year-varying), so
# only intercept rows with a `priors` entry are meaningful -- covariate slopes,
# `init`, and `est_phase` are rejected (prior-only), and a linkage on a fleet /
# predator that is not fit with a DM likelihood is rejected (the weight is not
# estimated there, so the prior would have no effect).
# =============================================================================

# jnll row for the linkage-table priors, robust to row-name vs index.
comp_linkage_prior_nll <- function(fit) {
  jc <- fit$quantities$jnll_comp
  rn <- rownames(jc)
  if (!is.null(rn) && "Linkage-table priors" %in% rn) {
    return(sum(jc["Linkage-table priors", ]))
  }
  sum(jc[20, ])
}

# ---------------------------------------------------------------------------
# Fast, build-only (estimateMode = 3) checks: no optimization, so the prior is
# evaluated at the default comp_weights and can be pinned to machine precision.
# ---------------------------------------------------------------------------
testthat::test_that("comp DM gamma prior is evaluated on exp(comp_weights) per DM fleet", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  d <- Rceattle::BS2017SS
  comp_flts <- sort(unique(d$comp_data$Fleet_code[d$comp_data$Fleet_code > 0]))
  testthat::skip_if(length(comp_flts) == 0)
  d$fleet_control$Comp_distribution <- "MultinomialAFSC"
  d$fleet_control$Comp_distribution[comp_flts] <- "DirichletMultinomial"

  shape <- 2; rate <- 0.5
  cf <- Rceattle::build_composition(linkages = list(
    theta_comp = Rceattle::linkage_spec(
      ~ 1, by = ~ fleet, fleet = comp_flts,
      priors = list(`(Intercept)` = gamma(shape, rate)))))

  fit <- suppressMessages(Rceattle::fit_mod(
    data_list = d, estimateMode = 3, msmMode = 0, compFun = cf,
    fit_control = Rceattle::fit_control(phase = FALSE, verbose = 0)))

  # The DM weight is the base parameter; the intercept coefficient is mapped
  # out at 0 (carried entirely by comp_weights).
  testthat::expect_true(all(fit$estimated_params$beta_linkage == 0))

  # Prior on the NATURAL-scale DM parameter theta = exp(comp_weights), summed
  # over the DM fleets: dgamma(theta, shape, scale = 1/rate). R's gamma() prior
  # helper stores (shape, rate); the cpp uses dgamma(b_nat, p1, 1/p2).
  b_nat    <- exp(as.numeric(fit$estimated_params$comp_weights[comp_flts]))
  expected <- -sum(stats::dgamma(b_nat, shape = shape, scale = 1 / rate, log = TRUE))
  testthat::expect_equal(comp_linkage_prior_nll(fit), expected, tolerance = 1e-8)

  # Negative control: NOT evaluated on the mapped-out beta_linkage = 0
  # (which would give exp(0) = 1 for every fleet).
  broken <- -length(comp_flts) * stats::dgamma(1, shape = shape, scale = 1 / rate, log = TRUE)
  testthat::expect_false(isTRUE(all.equal(comp_linkage_prior_nll(fit), broken)))
})

testthat::test_that("theta_diet prior is evaluated on exp(diet_comp_weights) per predator", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  # Predator diet is keyed by species, not fleet, so theta_diet exercises the
  # per-predator (sp_idx) branch that theta_comp / theta_caal (fleet-keyed) do
  # not. make_msm_test_data() carries diet_data; switch its likelihood to DM.
  d <- make_msm_test_data()$data_list
  d$Diet_distribution <- rep(1L, d$nspp)   # 1 = DirichletMultinomial for diet

  shape <- 2; rate <- 0.5
  cf <- Rceattle::build_composition(linkages = list(
    theta_diet = Rceattle::linkage_spec(
      ~ 1, by = ~ species, species = seq_len(d$nspp),
      priors = list(`(Intercept)` = gamma(shape, rate)))))

  fit <- suppressMessages(Rceattle::fit_mod(
    data_list = d, estimateMode = 3, msmMode = 1, compFun = cf,
    fit_control = Rceattle::fit_control(phase = FALSE, verbose = 0)))

  b_nat    <- exp(as.numeric(fit$estimated_params$diet_comp_weights))
  expected <- -sum(stats::dgamma(b_nat, shape = shape, scale = 1 / rate, log = TRUE))
  testthat::expect_equal(comp_linkage_prior_nll(fit), expected, tolerance = 1e-8)
})

testthat::test_that("no comp linkage leaves the linkage-prior row at 0", {
  testthat::skip_if_not_installed("TMB")
  testthat::skip_if_not_installed("Rceattle")

  fit <- suppressMessages(Rceattle::fit_mod(
    data_list = Rceattle::BS2017SS, estimateMode = 3, msmMode = 0,
    fit_control = Rceattle::fit_control(phase = FALSE, verbose = 0)))
  testthat::expect_equal(comp_linkage_prior_nll(fit), 0)
})

# ---------------------------------------------------------------------------
# Prior-only guards: these must error before MakeADFun, so they need no fit and
# run fast (no skip_on_cran).
# ---------------------------------------------------------------------------
testthat::test_that("a comp linkage on a non-DM fleet is rejected", {
  testthat::skip_if_not_installed("Rceattle")
  d <- Rceattle::BS2017SS
  comp_flts <- sort(unique(d$comp_data$Fleet_code[d$comp_data$Fleet_code > 0]))
  testthat::skip_if(length(comp_flts) == 0)
  # leave Comp_distribution at its (non-DM) default
  testthat::expect_error(
    suppressMessages(Rceattle::fit_mod(
      data_list = d, estimateMode = 3, msmMode = 0,
      compFun = Rceattle::build_composition(linkages = list(
        theta_comp = Rceattle::linkage_spec(
          ~ 1, by = ~ fleet, fleet = comp_flts,
          priors = list(`(Intercept)` = gamma(2, 0.5))))),
      fit_control = Rceattle::fit_control(phase = FALSE, verbose = 0))),
    "DirichletMultinomial")
})

testthat::test_that("a comp covariate slope (non-intercept) is rejected as prior-only", {
  testthat::skip_if_not_installed("Rceattle")
  d <- Rceattle::BS2017SS
  comp_flts <- sort(unique(d$comp_data$Fleet_code[d$comp_data$Fleet_code > 0]))
  testthat::skip_if(length(comp_flts) == 0)
  d$fleet_control$Comp_distribution <- "MultinomialAFSC"
  d$fleet_control$Comp_distribution[comp_flts] <- "DirichletMultinomial"
  d$env_data <- data.frame(Year = d$styr:d$endyr, temp = 0)
  testthat::expect_error(
    suppressMessages(Rceattle::fit_mod(
      data_list = d, estimateMode = 3, msmMode = 0,
      compFun = Rceattle::build_composition(linkages = list(
        theta_comp = Rceattle::linkage_spec(
          ~ temp, by = ~ fleet, fleet = comp_flts))),
      fit_control = Rceattle::fit_control(phase = FALSE, verbose = 0))),
    "prior-only")
})

testthat::test_that("an explicit by = ~ species on theta_comp is rejected (fleet-indexed)", {
  testthat::skip_if_not_installed("Rceattle")
  # comp/caal weights are fleet-indexed. Omitting `by` now defaults to ~ fleet
  # (the base stratum), so it targets the fleets correctly. An *explicit*
  # by = ~ species would emit fleet = NA rows the cpp collapses to fleet 1 --
  # silently mis-targeting the prior -- and is still rejected.
  d <- Rceattle::BS2017SS
  comp_flts <- sort(unique(d$comp_data$Fleet_code[d$comp_data$Fleet_code > 0]))
  testthat::skip_if(length(comp_flts) == 0)
  d$fleet_control$Comp_distribution <- "MultinomialAFSC"
  d$fleet_control$Comp_distribution[comp_flts] <- "DirichletMultinomial"
  testthat::expect_error(
    suppressMessages(Rceattle::fit_mod(
      data_list = d, estimateMode = 3, msmMode = 0,
      compFun = Rceattle::build_composition(linkages = list(
        theta_comp = Rceattle::linkage_spec(
          ~ 1, by = ~ species, priors = list(`(Intercept)` = gamma(2, 0.5))))),
      fit_control = Rceattle::fit_control(phase = FALSE, verbose = 0))),
    "fleet-indexed")
})

testthat::test_that("`init` / `est_phase` on a comp spec are rejected as prior-only", {
  testthat::skip_if_not_installed("Rceattle")
  d <- Rceattle::BS2017SS
  comp_flts <- sort(unique(d$comp_data$Fleet_code[d$comp_data$Fleet_code > 0]))
  testthat::skip_if(length(comp_flts) == 0)
  d$fleet_control$Comp_distribution <- "MultinomialAFSC"
  d$fleet_control$Comp_distribution[comp_flts] <- "DirichletMultinomial"

  make <- function(...) Rceattle::build_composition(linkages = list(
    theta_comp = Rceattle::linkage_spec(~ 1, by = ~ fleet, fleet = comp_flts, ...)))

  testthat::expect_error(
    suppressMessages(Rceattle::fit_mod(
      data_list = d, estimateMode = 3, msmMode = 0,
      compFun = make(init = list(`(Intercept)` = 3)),
      fit_control = Rceattle::fit_control(phase = FALSE, verbose = 0))),
    "init")
  testthat::expect_error(
    suppressMessages(Rceattle::fit_mod(
      data_list = d, estimateMode = 3, msmMode = 0,
      compFun = make(est_phase = 0L),
      fit_control = Rceattle::fit_control(phase = FALSE, verbose = 0))),
    "est_phase")
})
